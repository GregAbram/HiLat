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
