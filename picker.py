from paraview.util.vtkAlgorithm import *
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
            print('PICKER: as', active_selection.IDs, 'pids', pids)
            nselections = len(pids)

        if nselections > 0:
            try:
              selections = dsa.WrapDataObject(input).Points[pids,:]
              print('Picker selections', selections)
              np.savetxt(self.saved_selection_file, selections, delimiter = ",")
              print("selections saved")
            except:
              print("unable to save selection", selections, " to %s" % self.saved_selection_file)

        return 1
