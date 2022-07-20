from paraview.util.vtkAlgorithm import *
import paraview
import paraview.simple
from os import getenv

@smproxy.filter()
@smproperty.input(name="InputDataset", port_index=0)
@smdomain.datatype(dataTypes=["vtkUnstructuredGrid"], composite_data_supported=False)

class Disks(VTKPythonAlgorithmBase):
    def load_selections(self):
        import os
        if os.path.isfile(self.saved_selection_file):
          import numpy as np
          self.selections = np.loadtxt(self.saved_selection_file, delimiter = ",")
        else:
          print("no saved selections")
          self.selections = []

    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1, outputType="vtkUnstructuredGrid")
        self.saved_selection_file = getenv("HOME") + "/picks.csv"
        self.load_selections()
        self.offset = 0.001
        self.sides = 100
        self.Modified()

    @smproperty.xml("""
        <DoubleVectorProperty name="Offset" number_of_elements="1" default_values="0.001" command="SetOffset">
            <DoubleRangeDomain name="range" />
            <Documentation>Set offset of disks</Documentation>
        </DoubleVectorProperty>""")
    def SetOffset(self, o):
        self.offset = o
        self.Modified()

    @smproperty.xml("""
        <IntVectorProperty name="NSides" number_of_elements="1" default_values="100" command="SetNSides">
            <IntRangeDomain name="range" />
            <Documentation>Set offset of disks</Documentation>
        </IntVectorProperty>""")
    def SetNSides(self, n):
        self.sides = n
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
        from vtk import vtkUnstructuredGrid, VTK_TRIANGLE, vtkCellArray
        from vtk import vtkAppendFilter, vtkClipDataSet, vtkPlane
        from vtk.numpy_interface.dataset_adapter import numpy_support as ns
        from vtk.numpy_interface import dataset_adapter as dsa
        from math import sqrt
        import  numpy as np

        if len(self.selections) < 2:
          self.selections = []
          output.ShallowCopy(input)
          return 1

        if len(self.selections) > 3:
          print("using first 3 pick points")
          self.selections = self.selections[:3]

        if len(self.selections) == 2:
          corner = False
        else:
          corner = True

        slices = vtkAppendFilter()
        input = dsa.WrapDataObject(vtkUnstructuredGrid.GetData(inInfoVec[0], 0))
        output = vtkUnstructuredGrid.GetData(outInfoVec, 0)

        selections = self.selections

        def cross(a,b):
          return np.array([a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]])

        def sub(a,b):
          return np.array([a[0] - b[0], a[1] - b[1], a[2] - b[2]])

        def nrm(a):
          d = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
          return np.array([a[0]/d, a[1]/d, a[2]/d])

        ns = self.sides

        if corner:
          disks = vtkAppendFilter()
          corner = selections[0]
          R = sqrt(corner[0]*corner[0] + corner[1]*corner[1] + corner[2]*corner[2])
          corner = nrm(corner)
          first = True
          for p in selections[1:]:
            p0 = nrm(p)
            if first:
              slice_normal = nrm(cross(corner, p0))
            else:
              slice_normal = nrm(cross(p0, corner))
            v1 = nrm(cross(slice_normal, corner))
            A = 2 * np.pi * (np.arange(ns)/ns)
            xyz = (self.offset * slice_normal) + R*np.column_stack(corner[:,np.newaxis]*np.cos(A) + v1[:,np.newaxis]*np.sin(A)).astype('f4')
            xyz = dsa.numpy_support.numpy_to_vtk(xyz)
            ids = dsa.numpy_support.numpy_to_vtkIdTypeArray(np.column_stack(([3]*ns, [0]*ns, np.arange(ns), np.mod(np.arange(ns)+1, ns))))
            ca = vtkCellArray()
            ca.SetCells(ns, ids)
            co = dsa.numpy_support.numpy_to_vtkIdTypeArray(np.arange(ns).astype('i8')*4)
            ct =  dsa.numpyTovtkDataArray(np.array([VTK_TRIANGLE]*ns).astype('u1'))
            so = dsa.WrapDataObject(vtkUnstructuredGrid())
            so.Points = xyz
            so.VTKObject.SetCells(ct, co, ca)
            clipper = vtkClipDataSet()
            if first:
              clip_normal = nrm(cross(slice_normal, corner))
            else:
              clip_normal = nrm(cross(corner, slice_normal))
            plane0 = vtkPlane()
            plane0.SetOrigin(0.0, 0.0, 0.0)
            plane0.SetNormal(clip_normal)
            clipper.SetClipFunction(plane0)
            clipper.SetInputData(so.VTKObject)
            clipper.Update()
            disks.AddInputData(clipper.GetOutput())
            first = False
          disks.Update()
          output.ShallowCopy(disks.GetOutput())

        else:

          p0 = nrm(selections[0])
          p1 = nrm(selections[1])
          R = sqrt(p0[0]*p0[0] + p0[1]*p0[1] + p0[2]*p0[2])
          slice_normal = nrm(cross(p0, p1))
          v1 = nrm(cross(slice_normal, p0))
          A = 2 * np.pi * (np.arange(ns)/ns)
          xyz = (self.offset * slice_normal) + R*np.column_stack(p0[:,np.newaxis]*np.cos(A) + v1[:,np.newaxis]*np.sin(A)).astype('f4')
          xyz = dsa.numpy_support.numpy_to_vtk(xyz)
          ids = dsa.numpy_support.numpy_to_vtkIdTypeArray(np.column_stack(([3]*ns, [0]*ns, np.arange(ns), np.mod(np.arange(ns)+1, ns))))
          ca = vtkCellArray()
          ca.SetCells(ns, ids)
          co = dsa.numpy_support.numpy_to_vtkIdTypeArray(np.arange(ns).astype('i8')*4)
          ct =  dsa.numpyTovtkDataArray(np.array([VTK_TRIANGLE]*ns).astype('u1'))
          so = dsa.WrapDataObject(vtkUnstructuredGrid())
          so.Points = xyz
          so.VTKObject.SetCells(ct, co, ca)
          output.ShallowCopy(so.VTKObject)

        return 1
