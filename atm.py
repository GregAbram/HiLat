from mpi4py import MPI
import os, sys
from netCDF4 import Dataset
import numpy as np
import vtk
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.numpy_interface.dataset_adapter import numpy_support as ns
from glob import glob

def ll_to_xyz(la, lo):
  dla = (la / 90.0) * (np.pi / 2)
  dlo = (lo / 180.0) * np.pi
  x = (np.cos(dla) * np.cos(dlo))
  y = (np.cos(dla) * np.sin(dlo))
  z = np.sin(dla)
  return np.column_stack((x,y,z)).astype('f4')

def ll_to_ll0(la, lo):
  return np.column_stack((lo,la,[0]*len(la))).astype('f4')

def pull_mesh(meshfile, radius):
  vtu = vtk.vtkUnstructuredGrid()
  ds = Dataset(meshfile)
  ncells = ds.dimensions['grid_size'].size
  cola = np.array(ds['grid_corner_lat']).flatten()
  colo = np.array(ds['grid_corner_lon']).flatten()
  num_corners = colo.shape[0]
  corners = np.column_stack((cola, colo))   # The complete corner point list
  corners = (1000000*corners).astype('i4')  # Rounded
  uniq = np.unique(corners, axis=0, return_inverse=True)
  opts = uniq[0] / 1000000.0                # output unique points
  pmap = uniq[1]                            # pmap[i] is the index of the output point for the i'th input point
  quads = np.array([pmap[i] for i in np.arange(num_corners)]).reshape(-1,4)
  xyz = ll_to_xyz(opts[:,0], opts[:,1]) * radius
  pts = vtk.vtkPoints()
  vxyz = ns.numpy_to_vtk(xyz)
  pts.SetData(ns.numpy_to_vtk(xyz))
  vtu.SetPoints(pts)
  vtx = ns.numpy_to_vtkIdTypeArray(np.column_stack(([4]*ncells, quads)).flatten())
  carray = vtk.vtkCellArray()
  carray.SetCells(ncells, vtx)
  ctypes = ns.numpy_to_vtk(vtk.VTK_QUAD * np.ones(ncells).astype('u1'))
  coff = ns.numpy_to_vtkIdTypeArray(5 * np.arange(ncells))
  vtu.SetCells(ctypes, coff, carray)
  return vtu

def pull_data(vtu, variables, nc, template):
  ds = Dataset(nc)
  times = np.array(ds['time'])
  for i,time in enumerate(times): 
    hour = int(time * 24)
    fname = template % hour
    if not os.path.exists(fname):
      dvtu = vtk.vtkUnstructuredGrid()
      dvtu.DeepCopy(vtu)
      for variable in variables:
        array = np.array(ds[variable])[i,:].astype('f4')
        vtk_array = ns.numpy_to_vtk(array)
        vtk_array.SetName(variable)
        dvtu.GetCellData().AddArray(vtk_array)
      wrtr = vtk.vtkXMLUnstructuredGridWriter()
      wrtr.SetInputData(dvtu)
      wrtr.SetFileName(fname)
      wrtr.Write()
      del wrtr
      del dvtu
      print(fname, 'done')
    else:
      print(fname, 'skipped')

if __name__ == '__main__':

  variables = []
  meshfile = ""
  data_regex = ""
  output_template = "atm-%064d.vtu"
  radius = 6372000

  def syntax():
    if rank == 0:
      print("args:")
      print("  -meshfile filename   .nc file containing grid information (required)")
      print("  -datafile regex      regex for datafiles (required)")
      print("  -vars gridvars      comma separated array of data variables (all variables)")
      print("  -o template          template for output file name - one %s for timestep, .vtu extension ('atm-%04d-%02d.vtu')")
      print("  -R radius            scale value for radius (6372000)")
      sys.exit(1)

  if len(sys.argv) == 1:
    syntax()

  args = sys.argv[1:]
  while len(args) > 0:
    if args[0] == '-rs':
      rank = int(args[1])
      size = int(args[2])
      args = args[3:]
    elif args[0] == '-meshfile':
      meshfile = args[1]
      args = args[2:]
    elif args[0] == '-datafile':
      data_regex = args[1]
      args = args[2:]
    elif args[0] == '-vars':
      variables = args[1].split(',')
      print(variables)
      args = args[2:]
    elif args[0] == '-o':
      output_template = args[1]
      args = args[2:]
    elif args[0] == '-R':
      radius = float(args[1])
      args = args[2:]
    else:
      print('unknown arg:', args[0])
      syntax()

  if meshfile == "" or data_regex == "":
    syntax()

  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  size = comm.Get_size()

  vtu = pull_mesh(meshfile, radius)
  print("rank %d mesh ready" % rank)
  if rank == 0:
    wrtr = vtk.vtkXMLUnstructuredGridWriter()
    wrtr.SetFileName('grid.vtu')
    wrtr.SetInputData(vtu)
    wrtr.Write()
    print("grid written")

  dfiles = glob(data_regex)[rank::size]
  for dfile in dfiles:
    print("rank %d working on %s" % (rank, dfile))
    pull_data(vtu, variables, dfile, output_template)
    print("rank %d done with %s" % (rank, dfile))

