from paraview.util.vtkAlgorithm import *
import paraview
import paraview.simple

@smproxy.filter()
@smproperty.input(name="InputDataset", port_index=0)
@smdomain.datatype(dataTypes=["vtkUnstructuredGrid"], composite_data_supported=False)

class Clipper(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1, outputType="vtkUnstructuredGrid")
        self.save_selection_file = "none"
        self.selections = []

    @smproperty.stringvector(name="SavedSelection", default_values="none")
    def SetSavedSelection(self, value):
        self.saved_selection_file = value
        if self.selections == [] and self.saved_selection_file != None:
          import os
          if os.path.isfile(self.saved_selection_file):
            import numpy as np
            self.selections = np.loadtxt(self.saved_selection_file, delimiter = ",")
        else:
          if self.saved_selection_file == None:
            self.selections = []
          else:
            import numpy as np
            np.savetxt(self.saved_selection_file, self.selections, delimiter = ",")
        self.Modified()
        return

    @smproperty.intvector(name="ClearSelection", default_values=0)
    def SetClearSelection(self, value):
        print("CLEAR", value)
        if value > 0:
          self.selections = []
          print("clearing selections:", self.selections)
        self.Modified()
        return

    def FillInputPortInformation(self, port, info):
        info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid")
        return 1

    def RequestData(self, request, inInfoVec, outInfoVec):
        from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid
        from vtk import vtkAppendFilter, vtkCutter, vtkClipDataSet, vtkPlane
        from vtk.numpy_interface.dataset_adapter import numpy_support as ns
        from vtk.numpy_interface import dataset_adapter as dsa
        from math import sqrt
        import  numpy as np

        input = dsa.WrapDataObject(vtkUnstructuredGrid.GetData(inInfoVec[0], 0))
        output = vtkUnstructuredGrid.GetData(outInfoVec, 0)

        if 'selections' not in input.FieldData.keys():
          return 1

        selections = input.FieldData['selections']

        def cross(a,b):
          return np.array([a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]])

        def sub(a,b):
          return np.array([a[0] - b[0], a[1] - b[1], a[2] - b[2]])

        def nrm(a):
          d = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
          return np.array([a[0]/d, a[1]/d, a[2]/d])

        nsamples = 100
        disks = vtkAppendFilter()

        corner = selections[0]
        p1 = selections[1]
        p2 = selections[2]

        clip0 = vtkClipDataSet()
        clip0.SetInputData(input.VTKObject)
        cut0 = cross(corner, p1)
        plane0 = vtkPlane()
        plane0.SetOrigin(0.0, 0.0, 0.0)
        plane0.SetNormal(cut0) 
        clip0.SetClipFunction(plane0)

        clip1 = vtkClipDataSet()
        clip1.SetInputConnection(clip0.GetOutputPort())
        cut1 = cross(p2, corner)
        plane1 = vtkPlane()
        plane1.SetOrigin(0.0, 0.0, 0.0)
        plane1.SetNormal(cut1) 
        clip1.SetClipFunction(plane1)

        clip1.Update()
        output.ShallowCopy(clip1.GetOutput())

        del clip0
        del plane0
        del clip1
        del plane1

        return 1



