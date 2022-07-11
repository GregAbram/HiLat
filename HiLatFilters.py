from paraview.util.vtkAlgorithm import *

#------------------------------------------------------------------------------
# A filter example.
#------------------------------------------------------------------------------
@smproxy.filter()
@smproperty.input(name="InputDataset", port_index=0)
@smdomain.datatype(dataTypes=["vtkUnstructuredGrid"], composite_data_supported=False)

class TruncatePathLine(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1, outputType="vtkUnstructuredGrid")
        self.max_age = 100

    def FillInputPortInformation(self, port, info):
        info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid")
        return 1

    @smproperty.intvector(name="MaxAge", label="MaxAge", default_values=100)
    @smdomain.intrange(min=10, max=1000)
    def SetMaxAge(self, x):
        self.max_age = x
        self.Modified()

    def RequestData(self, request, inInfoVec, outInfoVec):
        import numpy as np
        from vtk import vtkUnstructuredGrid, vtkCellArray, VTK_POLY_LINE
        from vtk.numpy_interface import dataset_adapter as dsa

        inpt = dsa.WrapDataObject(vtkUnstructuredGrid.GetData(inInfoVec[0], 0))
        outpt = vtkUnstructuredGrid.GetData(outInfoVec, 0)

        if 'age' in inpt.PointData.keys():
          age = inpt.PointData['age']
        elif 'timestamp' in inpt.PointData.keys():
            ts = inpt.PointData['timestamp']
            max_ts = np.max(ts)
            age = max_ts - ts
        else:
          print("neither an age or timestamp variable!")

        cl = [inpt.CellLocations[i] + i for i in range(inpt.GetNumberOfCells())]
        ck = inpt.Cells[cl]
        cells = [inpt.Cells[cl[i]+1:cl[i]+ck[i]+1] for i in range(len(cl))]

        new_cells = []
        new_cell_locations = []
        for cell in cells: 
          cell_age = age[cell]
          if cell_age[-1] < self.max_age:
              k = len(cell) - np.argmax(cell_age < self.max_age)
              if k > 1:
                new_cell_locations = new_cell_locations + [len(new_cells)]
                new_cells.append(list(cell[-k:]))
        new_cell_locations = np.array(new_cell_locations)

        n_new_cells = len(new_cells)
        print(age[cells[0]])
        print('aaaaaa')
        print(age[new_cells[0]])

        referenced_points = np.hstack(new_cells)
        referenced_points_mask = np.zeros(inpt.GetNumberOfPoints())
        referenced_points_mask[referenced_points] = 1
        referenced_points_map = (np.cumsum(referenced_points_mask) - 1).astype('i4')

        mapped_new_cells = [referenced_points_map[cell] for cell in new_cells]

        numbered_mapped_new_cells = [[len(cell)] + list(cell) for cell in mapped_new_cells]
        stacked_numbered_mapped_new_cells = np.hstack(numbered_mapped_new_cells)
        clipped = dsa.WrapDataObject(vtkUnstructuredGrid())

        clipped.Points = inpt.Points[referenced_points_mask == 1]

        ct = dsa.numpy_support.numpy_to_vtk(np.array([VTK_POLY_LINE]*len(new_cells)).astype('u1'))
        co = dsa.numpy_support.numpy_to_vtkIdTypeArray(new_cell_locations)
        nc = dsa.numpy_support.numpy_to_vtkIdTypeArray(stacked_numbered_mapped_new_cells)
        ca = vtkCellArray()
        ca.SetCells(n_new_cells, nc)
        clipped.VTKObject.SetCells(ct, co, ca)

        for array_name in inpt.PointData.keys():
            clipped.PointData.append(inpt.PointData[array_name][referenced_points_mask == 1], array_name)
            
        clipped.PointData.append(age[referenced_points_mask == 1], 'age')

        outpt.ShallowCopy(clipped.VTKObject)
        del clipped
        return 1
import paraview
import paraview.simple
from os import getenv

@smproxy.filter()
@smproperty.input(name="InputDataset", port_index=0)
@smdomain.datatype(dataTypes=["vtkUnstructuredGrid"], composite_data_supported=False)

class Clipper(VTKPythonAlgorithmBase):
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

        input = dsa.WrapDataObject(vtkUnstructuredGrid.GetData(inInfoVec[0], 0))
        output = vtkUnstructuredGrid.GetData(outInfoVec, 0)

        print('hello', self.selections)

        if len(self.selections) < 2:
          self.selections = []
          output.ShallowCopy(input.VTKObject)
          return 1

        selections = self.selections

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
        from vtk import vtkUnstructuredGrid, VTK_TRIANGLE, vtkCellArray
        from vtk import vtkAppendFilter, vtkClipDataSet, vtkPlane
        from vtk.numpy_interface.dataset_adapter import numpy_support as ns
        from vtk.numpy_interface import dataset_adapter as dsa
        from math import sqrt
        import  numpy as np

        print('hello', self.selections)

        if len(self.selections) < 2:
          self.selections = []
          output.ShallowCopy(input)
          return 1

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

        nsamples = 100
        disks = vtkAppendFilter()

        corner = selections[0]
        R = sqrt(corner[0]*corner[0] + corner[1]*corner[1] + corner[2]*corner[2])
        print("R=",R)
        corner = nrm(corner)
        first = True
        for p in selections[1:]:
          p0 = nrm(p)
          if first:
            slice_normal = nrm(cross(corner, p0))
          else:
            slice_normal = nrm(cross(p0, corner))
          v1 = nrm(cross(slice_normal, corner))
          A = 2 * np.pi * (np.arange(nsamples)/nsamples)
          xyz = (0.001 * slice_normal) + R*np.column_stack(corner[:,np.newaxis]*np.cos(A) + v1[:,np.newaxis]*np.sin(A)).astype('f4')
          xyz = dsa.numpy_support.numpy_to_vtk(xyz)
          ids = dsa.numpy_support.numpy_to_vtkIdTypeArray(np.column_stack(([3]*nsamples, [0]*nsamples, np.arange(nsamples), np.mod(np.arange(nsamples)+1, nsamples))))
          ca = vtkCellArray()
          ca.SetCells(nsamples, ids)
          co = dsa.numpy_support.numpy_to_vtkIdTypeArray(np.arange(nsamples).astype('i8')*4)
          ct =  dsa.numpyTovtkDataArray(np.array([VTK_TRIANGLE]*nsamples).astype('u1'))
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
        return 1
import paraview
import paraview.simple

@smproxy.filter()
@smproperty.input(name="InputDataset", port_index=0)
@smdomain.datatype(dataTypes=["vtkUnstructuredGrid"], composite_data_supported=False)

class DScale(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1, outputType="vtkUnstructuredGrid")
        self.dscale = 3

    @smproperty.doublevector(name="Scale", default_values=[100])
    @smdomain.doublerange()
    def SetScale(self, d):
        self.dscale = d
        self.Modified()

    def FillInputPortInformation(self, port, info):
        info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid")
        return 1

    def RequestData(self, request, inInfoVec, outInfoVec):
        from vtk import vtkUnstructuredGrid, vtkPoints
        from vtk.numpy_interface import dataset_adapter as dsa
        import numpy as np

        input = vtkUnstructuredGrid.GetData(inInfoVec[0], 0)
        output = vtkUnstructuredGrid.GetData(outInfoVec, 0)

        # print("AAA", input.GetNumberOfPoints(), input.GetNumberOfCells())

        P = dsa.numpy_support.vtk_to_numpy(input.GetPoints().GetData())
        # print("P")
        # for p in P:
          # print(p)

        # Distance each point is from the center of the Earth
        W = np.linalg.norm(P, axis=1)
        np.where(W == 0, 1, W)

        # Normalized vector from each point to the center of the earth
        N = P / W[:,np.newaxis]

        # Wmax is the max radius of the data set
        Wmax = np.max(W)
        Wmin = np.min(W)
        print("MM", Wmax, Wmin)

        # proportional depth 
        D = (Wmax - W) / Wmax

        # scale it
        D = self.dscale*D
        # print("D2")
        # for d in D:
          # print(d)
  
        output.ShallowCopy(input)

        O = Wmax * (1.0 - D[:,np.newaxis]) * N
        # print("O")
        # for o in O:
          # print(o)
  
        pts = vtkPoints()
        pts.SetData(dsa.numpy_support.numpy_to_vtk(O))
        output.SetPoints(pts)

        # print("BBB", output.GetNumberOfPoints(), output.GetNumberOfCells())
        return 1

import paraview
import paraview.simple
import os

@smproxy.filter()
@smproperty.input(name="InputDataset", port_index=0)
@smdomain.datatype(dataTypes=["vtkUnstructuredGrid"], composite_data_supported=False)

class Picker(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1, outputType="vtkUnstructuredGrid")
        self.saved_selection_file = os.getenv('HOME') + '/picks.csv'

    @smproperty.stringvector(name="SavedSelection", default_values=os.getenv('HOME')+"/picks.csv")
    def SetSavedSelection(self, value):
        self.saved_selection_file = value
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

        output.ShallowCopy(input)

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
            try:
              selections = dsa.WrapDataObject(input).Points[pids,:]
              np.savetxt(self.saved_selection_file, selections, delimiter = ",")
              print("selections saved")
            except:
              print("unable to save selection", selections, " to %s" % self.saved_selection_file)

        return 1



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



