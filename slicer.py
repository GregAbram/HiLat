from paraview.util.vtkAlgorithm import *
import paraview
import paraview.simple

@smproxy.filter()
@smproperty.input(name="InputDataset", port_index=0)
@smdomain.datatype(dataTypes=["vtkUnstructuredGrid"], composite_data_supported=False)

class Slicer(VTKPythonAlgorithmBase):
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
        input = vtkUnstructuredGrid.GetData(inInfoVec[0], 0)
        output = vtkUnstructuredGrid.GetData(outInfoVec, 0)

        # If there are click points AND there are >2, they predominate

        print('hello', self.selections)

        nselections = -1
        proxy = paraview.simple.GetActiveSource()
        if proxy is not None:
          active_selection = proxy.GetSelectionInput(proxy.Port)
          if (active_selection is not None) and (len(active_selection.IDs) > 0):
            pids = active_selection.IDs[1::2]
            nselections = len(pids)

        if nselections > 0:
          if nselections < 2:
            print("need at least two selection points")
          else:
            self.selections = dsa.WrapDataObject(input).Points[pids,:]

        if len(self.selections) < 2:
          self.selections = []
          output.ShallowCopy(input)
        else:
          slices = vtkAppendFilter()
          p0 = self.selections[0]
          R = sqrt(p0[0]*p0[0] + p0[1]*p0[1] + p0[2]*p0[2])

          def cross(a,b):
            return[a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]
      
          for i,p1 in enumerate(self.selections[1:]):
            cutter = vtkCutter()
            cutter.SetInputData(input)
            cut_normal = cross(p0, p1)
            plane0 = vtkPlane()
            plane0.SetOrigin(0.0, 0.0, 0.0)
            plane0.SetNormal(cut_normal) 
            cutter.SetCutFunction(plane0)
            clipper = vtkClipDataSet()
            clipper.SetInputConnection(cutter.GetOutputPort())
            clip_normal = cross(cut_normal, p0)
            plane1 = vtkPlane()
            plane1.SetOrigin(0.0, 0.0, 0.0)
            plane1.SetNormal(clip_normal)
            clipper.SetClipFunction(plane1)
            clipper.Update()
            slice = clipper.GetOutput()
            ncells = slice.GetNumberOfCells()
            tArray = dsa.numpy_support.numpy_to_vtk((i*np.ones(ncells)).astype('u1'))
            tArray.SetName('tag')
            slice.GetCellData().AddArray(tArray)
            slices.AddInputData(slice)
            del cutter
            del plane0
            del clipper
            del plane1

          slices.Update()
          output.ShallowCopy(slices.GetOutput())
          selections = dsa.numpy_support.numpy_to_vtk(self.selections)
          selections.SetName("selections")
          output.GetFieldData().AddArray(selections)

          if self.saved_selection_file:
            try:
              np.savetxt(self.saved_selection_file, self.selections, delimiter = ",")
            except:
              print("unable to save selection", self.selections, " to %s" % self.saved_selection_file)

        return 1



