from paraview.util.vtkAlgorithm import *
import paraview
import paraview.simple
from os import getenv

@smproxy.filter()
@smproperty.input(name="InputDataset", port_index=0)
@smdomain.datatype(dataTypes=["vtkUnstructuredGrid"], composite_data_supported=False)

class Slicer(VTKPythonAlgorithmBase):
    def load_selections(self):
        import os
        if os.path.isfile(self.saved_selection_file):
          import numpy as np
          self.selections = np.loadtxt(self.saved_selection_file, delimiter = ",")
        else:
          print("no saved selections")
          self.selections = []
        self.Modified()

    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1, outputType="vtkUnstructuredGrid")
        self.saved_selection_file = getenv("HOME") + "/picks.csv"
        self.load_selections()

    @smproperty.stringvector(name="SavedSelection", default_values=getenv("HOME") + "/picks.csv")
    def SetSavedSelection(self, value):
        self.saved_selection_file = value
        self.load_selections()
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

        print('hello', self.selections)

        if len(self.selections) < 2:
          self.selections = []
          output.ShallowCopy(input)
          return 1

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

        return 1



