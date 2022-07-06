import os, sys
from  mpi4py import MPI
from glob import glob
from vtk import *
import numpy as np
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.numpy_interface.dataset_adapter import numpy_support as ns
from netCDF4 import Dataset

month_offset = [0, 31,  59,  90, 120, 151, 181, 212, 243, 273, 304, 334]

# This will convert MPAS netCDF files to VTK unstructured grids and write them
# in the .vtu format.   
#
# For the purposes of the comments herein, the term 'MPAS grid' will refer to the
# mostly-hex grid that MPAS operates on, and the 'VTK grid' will refer to the 
# all-triangle VTK-format grid
#
# Every non-edge vertex of the MPAS grid (eg. with 3 incident MPAS cells) will 
# correspond to a triangle in the VTK grid, linking the centers of the three
# MPAS cells that are incident on the vertex. Similarly, every cell of the MPAS

# grid corresponds to a vertex of the VTK grid, where the vertex is at the
# center of the MPAS cell.
#
# Variables on the MPAS grid that are defined on the MPAS cells are directly
# transported to the vertices of the VTK grid.   Variables that are defined 
# on the vertices of the MPAS grid are migrated to the cell centers (by averaging
# the vertex values for the cell) and then transported to the vertices of the 
# VTK grid.
#
# NOTE: this will autodetect results and not re-create them!

def MPAS_to_VTU(rank, size, meshfile, data_regex, datavars, gridvars, layer, time_int, time_dim, output_template, uvw):

    print(data_regex)
    ncfiles = sorted(glob(data_regex))

    un_done_ncfiles = []
    done_ncfiles = []

    # Go through datafiles to check which have already been created.   While at it,
    # do several checks on data consistency.  Here we trim to the set of nc files that
    # this rank will consider. 

    print(ncfiles)
    for ncfile in ncfiles[rank::size]:
      ds = Dataset(ncfile)

      if time_dim not in ds.dimensions.keys():
        print("time dimension ", timedim, "not found")
        sys.exit(1)

      timesteps = ds.dimensions[time_dim].size

      # What are the timesteps in the file?   If the corresponding output
      # files are all present, we don't 'need' the nc file

      y,m,d = [int(i) for i in ncfile.split('.')[-2].split('_')[-2].split('-')]
      print(y,m,d)

      day = month_offset[(m - 1)] + (d - 1)
      hour = day * 24

      for tindx in range(timesteps):

        tstep = hour + tindx*time_int
        fname = output_template % tstep
        print(fname)

        if not os.path.isfile(fname):
          need = 1
          break

      if need:
        un_done_ncfiles.append(ncfile)
      else:
        done_ncfiles.append(ncfile)

      tmp = []
      for name in datavars:
        if name not in ds.variables.keys():
            print("var ", name, "not in dataset",  ds.variables.keys())
        else:
          variable = ds.variables[name]
          if len(variable.dimensions) == 2 and (variable.dimensions[0] != 'Time' or (variable.dimensions[1] != 'nCells' and variable.dimensions[1] != 'nVertices')):
            print("requested var ", name, "not correct shape")
          elif len(variable.dimensions) == 3 and (variable.dimensions[0] != 'Time' or (variable.dimensions[1] != 'nCells' and variable.dimensions[1] != 'nVertices') or (variable.dimensions[2][:11] != 'nVertLevels')):
            print("requested var ", name, "not correct shape")
          else:
            tmp.append(name)

      datavars = tmp

    print(rank, "doing", un_done_ncfiles)

    if len(un_done_ncfiles) == 0:
      print("nothing to do")
      sys.exit(1)

    # Now we load the mesh, creating an unstructured grid with triangles 
    # joining the hex-cell centers.  In this grid, vertices correspond to
    # cells in the hex grid and triangles correspond to vertices of the 
    # hex grid that have 3 incident edges - eg. are internal to the hex grid

    ds = Dataset(meshfile)

    # points of VTK grid are the cell centers of the MPAS grid, which
    # are included as variables in the netCDF.

    xyz = np.column_stack((np.array(ds['xCell']), np.array(ds['yCell']), np.array(ds['zCell'])))

    vxyz = ns.numpy_to_vtk(xyz)
    npts = len(xyz)

    grid = vtkUnstructuredGrid()
    pts = vtkPoints()
    pts.SetData(vxyz)
    grid.SetPoints(pts)

    # Note - we'll need to use the verticesOnCell as indices if/when
    # we need to merge the MPAS cell vertex data elements to the MPAS
    # cell center, where we'll use them on the VTK grid

    voc = np.array(ds.variables['verticesOnCell'])
    neoc = np.array(ds.variables['nEdgesOnCell'])

    # We'll use the 'cellsOnVertex' array to create the triangles.

    cov = np.array(ds.variables['cellsOnVertex'])

    # Which elements of the cov array have a zero, indicating a
    # missing adjacent cell?   These are boundary vertices (of the
    # MPAS grid) which won't correspond to an output triangle.

    boundary_vertices = np.unique(np.where(cov == 0))

    # The VTK triangle k corresponds to the k'th non-boundary 
    # vertex of the MPAS grid, and will join the VTK vertices 
    # corresponding to the 3 MPAS cells incident on that MPAS 
    # vertex.   Subtract 1 due to FORTRAN indexing in the MPAS
    # grid

    tris = np.delete(cov, boundary_vertices, axis=0) - 1

    # Create a vtk triangle list: 3 p0 p1 p2 3 q0 q1 q2...
    # then add the triangles to the vtu

    tris = np.column_stack(([3]*tris.shape[0], tris))
    ntris = tris.shape[0]

    vtris = ns.numpy_to_vtkIdTypeArray(tris.flatten())
    carray = vtkCellArray()
    carray.SetCells(ntris, vtris)
    ctypes = ns.numpy_to_vtk(VTK_TRIANGLE * np.ones(ntris).astype('u1'))
    coff = ns.numpy_to_vtkIdTypeArray(4 * np.arange(ntris))
    grid.SetCells(ctypes, coff, carray)

    if npts != grid.GetNumberOfPoints() or ntris != grid.GetNumberOfCells():
      print('grid error')
      sys.exit(1)

    if uvw:
      # Add UVW so we can later build vectors

      U = np.column_stack((-xyz[:,1], xyz[:,0], [0]*npts))
      U = U / np.linalg.norm(U, axis=1)[:,np.newaxis]
      W = xyz / np.linalg.norm(xyz, axis=1)[:,np.newaxis]
      V = np.cross(W, U)

      vU = ns.numpy_to_vtk(U.astype('f4'))
      vU.SetName('U')
      grid.GetPointData().AddArray(vU)

      vV = ns.numpy_to_vtk(V.astype('f4'))
      vV.SetName('V')
      grid.GetPointData().AddArray(vV)

      vW = ns.numpy_to_vtk(W.astype('f4'))
      vW.SetName('W')
      grid.GetPointData().AddArray(vW)

    for vname in gridvars:

      var = ds.variables[vname]
      dimensions = var.dimensions

      nvar = np.array(var)

      if dimensions[0] == 'Time':  # got to be 0... the mesh oughtn't contain time varying data
        nvar = nvar[0]
        dimensions = dimensions[1:]

      if len(dimensions) == 2 and dimensions[-1][:11] == 'nVertLevels':
        nvar = nvar[:,layer]
        dimensions = dimensions[:-1]

      if len(dimensions) != 1 or (dimensions[0] != 'nVertices' and dimensions[0] != 'nCells'):
        print("gvar shape error: ", vname, var.dimensions)
        sys.exit(0)

      if var.dimensions[0] == 'nVertices':

        # Then the variable is vertex-associated on the MPAS grid and we need
        # to average each cell's vertex values to get a MPAS cell-associated
        # value. The np.where creates a new variable [nCells, nEdgesOnCells] with 0
        # where voc[cell,vertex] == 0 (eg. this cell has fewer than nEdgesOnCells
        # vertices) or the vertex value otherwise.   Remember FORTRAN.  Then sum
        # along rows and divide by the actual number of vertices of each cell to
        # get the result that we can associate with the vertices of the VTK grid

        nvar = np.sum(np.where(voc == 0, 0, np.array(var)[voc-1]), axis=1) / neoc

      if len(nvar.shape) > 2:
        print("shape error:", vname, ds.variables[vname].dimensions, nvar.shape, dimensions)
      else:
        vvar = ns.numpy_to_vtk(nvar.astype('f4'))
        vvar.SetName(vname)
        grid.GetPointData().AddArray(vvar)

    print('grid ready,', npts, 'vertices', ntris, "cells")

    wrtr = vtkXMLUnstructuredGridWriter()

    # now for every nc file we need to process...

    for ncfile in  un_done_ncfiles:

      ds = Dataset(ncfile)

      timesteps = ds.dimensions[time_dim].size

      y,m,d = [int(i) for i in ncfile.split('.')[-2].split('_')[-2].split('-')]
      print(y,m,d)

      day = month_offset[(m - 1)] + (d - 1)
      hour = day * 24

      for tindx in range(timesteps):

        tstep = hour + tindx*time_int
        fname = output_template % tstep

        # While we know from above that *some* timestep in this nc file
        # hasn't been created, some may have been..

        if not os.path.isfile(fname):
          print("rank", rank, "working on", fname)

          # Create output vtu and copy grid and non-time-varying info in...

          out = vtk.vtkUnstructuredGrid()
          out.ShallowCopy(grid)

          # Use python data object wrapping for convenience:

          vout = dsa.WrapDataObject(out)
          p = vout.PointData
          c = vout.CellData
          vc = c.VTKObject

          for vname in datavars:


            var = ds.variables[vname]
            nvar = np.array(var)
            dims = var.dimensions

            if dims[0] == 'Time':
              nvar = nvar[tindx]
              dims = dims[1:]

            if dims[-1][:11] == 'nVertLevels':
              nvar = nvar[:,layer]
              dims = dims[:-1]

            if len(dims) != 1 or (dims[0] != 'nCells' and dims[0] != 'nVertices'):
              print("shape error:", vname, ds.variables[vname].dimensions)
            else:
              if dims[0] == 'nVertices':  # Then sum vertex values and divide
                nvar = np.sum(np.where(voc == 0, 0, nvar[voc-1]), axis=1) / neoc

              vvar = ns.numpy_to_vtk(nvar.astype('f4'))
              vvar.SetName(vname)
              out.GetPointData().AddArray(vvar)

          print('output:', vout.GetNumberOfPoints(), vout.GetNumberOfCells())
          wrtr.SetFileName(fname)
          wrtr.SetInputData(vout.VTKObject)
          wrtr.Write()
          print(fname, "done")

if __name__ == '__main__':

  rank = 0
  size = 1
  gridvars = []
  datavars = []
  meshfile = ""
  data_regex = ""
  time_int = 24
  time_dim = 'Time'
  output_template = "timestep_%08d.vtu"
  uvw = False
  layer=0

  def syntax():
    if rank == 0:
      print("args:")
      print("  -meshfile filename   .nc file containing grid information (required)")
      print("  -datafile regex      regex for datafiles (required)")
      print("  -gvars gridvars      comma separated array of data variables (all variables)")
      print("  -vars datavars       comma separated array of data variables (all variables)")
      print("  -o template          template for output file name - one %s for timestep, .vtu extension ('timestep_%08d.vtu')")
      print("  -tint  hours         interval between timesteps (24)")
      print("  -tdim  dimname       name of dimension containing the number of timesteps ('Time')")
      print("  -uvw                 add UVW tangent/radial vectors to output (no)")
      print("  -l layer             for multi-layer data, whch layer (0)")
      sys.exit(1)

  if len(sys.argv) == 1:
    syntax()

  # comm = MPI.COMM_WORLD
  # rank = comm.Get_rank()
  # size = comm.Get_size()
  rank = 0
  size = 1

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
    elif args[0] == '-gvars':
      gridvars = args[1].split(',')
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
    elif args[0] == '-uvw':
      uvw = True
      args = args[1:]
    elif args[0] == '-l':
      layer = int(args[1])
      args = args[2:]
    else:
      print('unknown arg:', args[0])
      syntax()
      
  if meshfile == "" or data_regex == "":
    syntax()
        

  MPAS_to_VTU(rank, size, meshfile, data_regex, datavars, gridvars, layer, time_int, time_dim, output_template, uvw)
