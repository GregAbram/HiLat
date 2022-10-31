from paraview.util.vtkAlgorithm import *
import paraview
import paraview.simple
from os import getenv

@smproxy.filter()
@smproperty.input(name="InputDataset", port_index=0)
@smdomain.datatype(dataTypes=["vtkUnstructuredGrid"], composite_data_supported=False)

class Tagger(VTKPythonAlgorithmBase):
    def load_selections(self):
        import os
        if os.path.isfile(self.saved_selection_file):
          import numpy as np
          self.selections = np.loadtxt(self.saved_selection_file, delimiter = ",").reshape(-1,3)
        else:
          print("no saved selections")
          self.selections = []
        self.reverse = 0
        self.Modified()

    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1, outputType="vtkUnstructuredGrid")
        self.saved_selection_file = getenv("HOME") + "/picks.csv"
        self.load_selections()

    @smproperty.xml("""
        <IntVectorProperty name="Reverse" number_of_elements="1" default_values="0" command="SetReverse">
            <Documentation>Reverse Single Slice</Documentation>
        </IntVectorProperty>""")
    def SetReverse(self, n):
        self.reverse = n
        self.Modified()

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

        print('aaa')
        input = dsa.WrapDataObject(vtkUnstructuredGrid.GetData(inInfoVec[0], 0))
        output = dsa.WrapDataObject(vtkUnstructuredGrid.GetData(outInfoVec, 0))

        output.VTKObject.ShallowCopy(input.VTKObject)

        if len(self.selections) == 0:
          self.selections = []
          return 1


        selections = self.selections.reshape(-1,3)

        if len(selections) == 1:
          selections = np.vstack(([0.0, 0.0, 1.0],selections))

        if len(selections) > 3:
          print("using first 3 pick points")
          selections = selections[:3]

        if len(selections) == 2:
          corner = False
        else:
          corner = True

        def cross(a,b):
          return np.array([a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]])

        def sub(a,b):
          return np.array([a[0] - b[0], a[1] - b[1], a[2] - b[2]])

        def nrm(a):
          d = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
          return np.array([a[0]/d, a[1]/d, a[2]/d])

        if corner:
          corner = selections[0]
          p1 = selections[1]
          p2 = selections[2]
          cut = cross(corner, p1)
          tag0 = np.sum(cut * input.Points, axis=1) > 0
          cut = cross(p2, corner)
          tag1 = np.sum(cut * input.Points, axis=1) > 0
          output.PointData.append(np.logical_or(tag0, tag1).astype('u1'), 'tag')
          del cut
        else:
          print("AAA")
          if self.reverse:
              p0 = nrm(selections[0])
              p1 = nrm(selections[1])
          else:
              p0 = nrm(selections[1])
              p1 = nrm(selections[0])
          cut = cross(p0, p1)
          a = cut * input.Points
          b = np.sum(a, axis=1)
          tags = (b > 0).astype('u1')
          output.PointData.append(tags, 'tag')
          del cut

        return 1
