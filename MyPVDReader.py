from paraview.util.vtkAlgorithm import *

def createModifiedCallback(anobject):
    import weakref
    weakref_obj = weakref.ref(anobject)
    anobject = None
    def _markmodified(*args, **kwars):
        o = weakref_obj()
        if o is not None:
            o.Modified()
    return _markmodified

@smproxy.source(name="MyPVDReader", label="My PVD Reader")

class MyPVDReader(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=0, nOutputPorts=1, outputType='vtkUnstructuredGrid')
        from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader
        self._filename = ""
        self._timestep = 0
        self._vtus = []
        self._ndata = None
        self._rdr = vtkXMLUnstructuredGridReader()
        from vtkmodules.vtkCommonCore import vtkDataArraySelection
        self._arrayselection = vtkDataArraySelection()
        self._arrayselection.AddObserver("ModifiedEvent", createModifiedCallback(self))

    def _get_array_selection(self):
        return self._arrayselection

    @smproperty.stringvector(name="FileName")
    @smdomain.filelist()
    @smhint.filechooser(extensions="pvd", file_description="PVD files")
    def SetFileName(self, name):
        if name != 'None' and self._filename != name:
            from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid
            self._filename = name
            directory = self._filename.rsplit('/', 1)[0]
            f = open(self._filename)
            self._vtus = [directory + '/' + l.strip().split('"')[5] for l in f if 'DataSet' in l]
            self._rdr.SetFileName(self._vtus[self._timestep])
            self._rdr.UpdateInformation()
            self._arrayselection.RemoveAllArrays()
            for i in range(self._rdr.GetNumberOfPointArrays()):
              r = self._rdr.GetPointArrayName(i)
              self._arrayselection.AddArray(r)
            self.Modified()

    @smproperty.intvector(name="timestep", default_values=0)
    @smdomain.intrange()
    def SetTimestep(self, t):
        self._timestep = t
        self.Modified()

    @smproperty.dataarrayselection(name="Arrays")
    def GetDataArraySelection(self):
        return self._get_array_selection()

    def RequestData(self, request, inInfoVec, outInfoVec):
        from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid
        from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader
        if self._timestep < 0 or self._timestep >= len(self._vtus):
          print("XXXXXXXXXXXXXX error accessing", self._timestep, 'of', self._filename)
          return False
        print("MYPVD reader: ", self._timestep, self._vtus[self._timestep])
        self._rdr.SetFileName(self._vtus[self._timestep])
        for i in range(self._arrayselection.GetNumberOfArrays()):
          a = self._arrayselection.GetArrayName(i)
          self._rdr.SetPointArrayStatus(a, self._arrayselection.ArrayIsEnabled(a))
        self._rdr.Update()
        output = vtkUnstructuredGrid.GetData(outInfoVec, 0)
        output.ShallowCopy(self._rdr.GetOutput())

        return 1

