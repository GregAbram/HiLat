from paraview.util.vtkAlgorithm import *
import paraview
import paraview.simple

@smproxy.filter()
@smproperty.input(name="InputDataset", port_index=0)
@smdomain.datatype(dataTypes=["vtkUnstructuredGrid"], composite_data_supported=False)

class Disks(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1, outputType="vtkUnstructuredGrid")

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

#         nselections = -1
#         proxy = paraview.simple.GetActiveSource()
#         if proxy is not None:
#           active_selection = proxy.GetSelectionInput(proxy.Port)
#           if (active_selection is not None) and (len(active_selection.IDs) > 0):
#             pids = active_selection.IDs[1::2]
#             nselections = len(pids)
#         if nselections < 2:
#           print("need at least two selection points")
#           output.ShallowCopy(input)
#         else:
#           selections = dsa.WrapDataObject(input).Points[pids,:]
#           slices = vtkAppendFilter()
#           p0 = selections[0]
#           R = sqrt(p0[0]*p0[0] + p0[1]*p0[1] + p0[2]*p0[2])
# 
#           def cross(a,b):
#             return[a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]
#       
#           for i,p1 in enumerate(selections[1:]):
#             cutter = vtkCutter()
#             cutter.SetInputData(input)
#             cut_normal = cross(p0, p1)
#             plane0 = vtkPlane()
#             plane0.SetOrigin(0.0, 0.0, 0.0)
#             plane0.SetNormal(cut_normal) 
#             cutter.SetCutFunction(plane0)
#             clipper = vtkClipDataSet()
#             clipper.SetInputConnection(cutter.GetOutputPort())
#             clip_normal = cross(cut_normal, p0)
#             plane1 = vtkPlane()
#             plane1.SetOrigin(0.0, 0.0, 0.0)
#             plane1.SetNormal(clip_normal)
#             clipper.SetClipFunction(plane1)
#             clipper.Update()
#             slice = clipper.GetOutput()
#             ncells = slice.GetNumberOfCells()
#             tArray = dsa.numpy_support.numpy_to_vtk((i*np.ones(ncells)).astype('u1'))
#             tArray.SetName('tag')
#             slice.GetCellData().AddArray(tArray)
#             slices.AddInputData(slice)
#             del cutter
#             del plane0
#             del clipper
#             del plane1
# 
#           slices.Update()
#           output.ShallowCopy(slices.GetOutput())
# 
#           selections = dsa.numpy_support.numpy_to_vtk(selections)
#           selections.SetName("selections")
#           output.GetFieldData().AddArray(selections)




