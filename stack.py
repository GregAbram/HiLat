import sys, pdb, os
from mpi4py import MPI
from netCDF4 import Dataset
from vtk import *
import numpy as np
from glob import glob
from vtk.numpy_interface.dataset_adapter import numpy_support as ns

month_offset = [0, 31,  59,  90, 120, 151, 181, 212, 243, 273, 304, 334]

def load_grid(gridfile, max_layers):
  gds = Dataset(gridfile)

  xyz = np.column_stack((np.array(gds['xCell']), np.array(gds['yCell']), np.array(gds['zCell']))).astype('f4')

  thickness = gds['layerThickness'][0]
  if max_layers > 0:
    thickness = thickness[:,:max_layers]

  bot = np.array(gds['bottomDepth'])

  # Thickness has one row for each cell center (eventual vertex of the tri-grid).
  # The i'th row of thickness are the vertical intervals of the vertices under point i,
  # stating at the topmost interval.  However, they need to be accumulated up from the
  # bottom  (I believe).

  # The following  reverses the *rows* of layerThickness, runs a row-wise cumsum, and
  # reverses them again. The i'th row now are the offsets of the *tops* of the intervals
  # from the bottom, starting from the topmost, which is thus more-or-less the surface.

  layer_offsets = np.fliplr(np.cumsum(np.fliplr(thickness), axis=1))

  # We need to add one more column: a zero at the end of each row so the lowest interval
  # starts at the bottom.

  layer_offsets = np.column_stack((layer_offsets, np.zeros(thickness.shape[0])))

  # Now we get the points on the bottom by subtracting the corresponding bottom depth
  # times the vector toward the center of the Earth from the surface point.  The bottom
  # depth is given in positive meters.

  R = xyz / np.linalg.norm(xyz, axis=1)[:,np.newaxis]
  bottom = xyz - bot[:,np.newaxis]*R

  # Now, for a given point on the bottom, we create a i'th layer point by adding the
  # i'th layer_offset of that point times the R vector to the bottom point

  layers = [bottom + layer_offsets[:,i][:,np.newaxis]*R for i in range(layer_offsets.shape[1])]

  # Now layers is a *list* of layers, with the topmost first.  vstack them to get the vertex
  # set for the volumetric data

  prism_points = np.vstack(layers)

  grid = vtkUnstructuredGrid()

  gridinfo = vtkIntArray()
  gridinfo.SetNumberOfComponents(1);
  gridinfo.InsertNextValue(thickness.shape[-1] + 2) # Number of levels - original number of layers + top and bottom
  gridinfo.InsertNextValue(thickness.shape[0])      # Points in each level of the original grid
  gridinfo.SetName('gridinfo')
  grid.GetFieldData().AddArray(gridinfo)

  pts = vtkPoints()
  pts.SetData(ns.numpy_to_vtk(prism_points.astype('f4')))
  grid.SetPoints(pts)

  cov = np.array(gds.variables['cellsOnVertex'])
  boundary_vertices = np.unique(np.where(cov == 0))
  tris = np.delete(cov, boundary_vertices, axis=0) - 1

  cells = [np.column_stack(([6]*tris.shape[0], tris+(i*xyz.shape[0]), tris+((i+1)*xyz.shape[0]))) for i in range(thickness.shape[-1])]
  cells = np.vstack(cells)
  ctype = VTK_WEDGE


  ncells = cells.shape[0]
  cverts = cells.shape[1]

  carray = vtkCellArray()
  carray.SetCells(ncells, ns.numpy_to_vtkIdTypeArray(cells.flatten()))
  ctypes = ns.numpy_to_vtk((ctype * np.ones(ncells)).astype('u1'))
  coff = ns.numpy_to_vtkIdTypeArray(cverts * np.arange(ncells))
  grid.SetCells(ctypes, coff, carray)

  return (grid, thickness)

def load_data(datafile, grid, thickness, variables, time_dim, output_template, max_layers):

  # Variables are cell-centered - now centered at the center of each interval under each
  # surface point.  To get the values at each *internal interface i* we compute a weighted
  # sum of the interval-centered values above and below.   The weights vary because the
  # distance between the interface points and the interval centers differs due to the varying
  # layer thicknesses

  # So for a given *internal* interface, the linear weight for the *above* value is
  # belowThickness / (aboveThickness + belowThickness) and the *below* value is
  # aboveThickness / (aboveThickness + belowThickness)

  # We then use the top inverval-center value for the surface vertex and the bottom interval-
  # center for the bottom.  Avoid divide by zero; won't matter since the numerators will be
  # zero

  delta = (thickness[:,1:] + thickness[:,:-1])
  delta = np.where(delta == 0, 1, delta)
  
  t0  = thickness[:,:-1] / delta
  t1 = thickness[:,1:] / delta

  y,m,d = [int(i) for i in datafile.split('_')[-2].split('.')[-1].split('-')]
  day = month_offset[(m - 1)] + (d - 1)

  print("DAY ", day)

  dds = Dataset(datafile)
  timesteps = dds.dimensions[time_dim].size

  wrtr = vtkXMLUnstructuredGridWriter()

  for tstep in range(timesteps):
    day_of_the_year = (day + tstep)*time_int
    fname = output_template % day_of_the_year

    if os.path.isfile(fname):
      print(fname, 'is already done')
    else:
      nvtu = vtkUnstructuredGrid()
      nvtu.ShallowCopy(grid)

      for vname in variables:
        var = dds[vname][tstep]
        if max_layers > 0:
          var = var[:,:max_layers]
        a = [var[:,0]]
        b = [t0[:,i]*var[:,i+1] + t1[:,i]*var[:,i] for i in range(var.shape[1]-1)]
        c = [var[:,-1]]
        tvar = np.vstack((a + b + c))
        print(vname, tvar.shape, tstep)
        tvar = ns.numpy_to_vtk(tvar.flatten().astype('f4'))
        tvar.SetName(vname)
        nvtu.GetPointData().AddArray(tvar)

      wrtr.SetInputData(nvtu)
      wrtr.SetFileName(fname)
      wrtr.Write()
      print(fname, 'done!')


if __name__ == '__main__':

  rank = 0
  size = 1
  datavars = []
  meshfile = ""
  data_regex = ""
  time_int = 24
  time_dim = 'Time'
  output_template = "timestep_%08d.vtu"
  max_layers = -1

  def syntax():
    if rank == 0:
      print("args:")
      print("  -meshfile filename   .nc file containing grid information (required)")
      print("  -datafile regex      regex for datafiles (required)")
      print("  -vars datavars       comma separated array of data variables (all variables)")
      print("  -o template          template for output file name - one %s for timestep, .vtu extension ('timestep_%08d.vtu')")
      print("  -tint  hours         interval between timesteps (24)")
      print("  -tdim  dimname       name of dimension containing the number of timesteps ('Time')")
      print("  -ml max_layers       limit the number of layers (all)")
      sys.exit(1)

  if len(sys.argv) == 1:
    syntax()

  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  size = comm.Get_size()

  print("rank:", rank)

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
      datavars = args[1].split(',')
      args = args[2:]
    elif args[0] == '-o':
      output_template = args[1]
      args = args[2:]
    elif args[0] == '-tint':
      time_int = int(args[1])
      args = args[2:]
    elif args[0] == '-tdim':
      time_dim = args[1]
      args = args[2:]
    elif args[0] == '-ml':
      max_layers = int(args[1])
      args = args[2:]
    else:
      print('unknown arg:', args[0])
      syntax()
      
  if meshfile == "" or data_regex == "":
    syntax()
        
  datafiles = glob(data_regex)[rank::size]

  if len(datafiles) > 0:

    grid, thickness = load_grid(meshfile, max_layers)
    print("mesh loaded")

    for datafile in datafiles:
      print("working on ", datafile)
      load_data(datafile, grid, thickness, datavars, time_dim, output_template, max_layers)
      print(datafile, 'done.')
