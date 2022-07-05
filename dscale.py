from paraview.util.vtkAlgorithm import *
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

