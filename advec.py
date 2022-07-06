from glob import glob
from vtk import *
import sys, os
import numpy as np
from vtk.util import numpy_support as ns
from vtk.numpy_interface import dataset_adapter as dsa
import pdb

N_SUBSTEPS = 10
MINLAT = 0.7
MAX_NUMBER_OF_POINTS_ALONG_LINES = 100
TERMINATION_STEPS = 200
NSEEDS = 500
FILENAME_TEMPLATE = 'tt_%05d.vtu'
WRITE_EVERY = 1
WIND_ARRAY='WBOT'
X = -1
template = ""

def sample(vtu, minlat, nsamp):
  a = vtu.Points
  b = a[a[:,2] > minlat]
  k = float(nsamp + 10) / b.shape[0] 
  s = np.array(b[np.random.random(b.shape[0]) < k])
  return s[:nsamp]

def xyz_to_latlon(xyz):
  x = xyz[:,0]
  y = xyz[:,1]
  z = xyz[:,2]
  d = np.sqrt(x*x + y*y + z*z)
  nxyz = xyz / d[:,np.newaxis]
  x = nxyz[:,0]
  y = nxyz[:,1]
  z = nxyz[:,2]
  lat = np.arcsin(z)
  lon = np.arctan2(-y, -x)+ np.pi
  return np.column_stack([lon, lat, [0]*lat.shape[0]])

def latlon(vtu):
  return xyz_to_latlon(vtu.Points)

def toLatLon(v):
  import numpy as np
  ll = latlon(v)
  lon = ll[:,0]
  cells = np.array(v.Cells).reshape(-1, 5)[:,1:]
  clon = lon[cells]
  cwid = np.max(clon, axis=1) - np.min(clon, axis=1)
  wrappers = cells[cwid > 1]
  non_wrappers = cells[cwid <= 1]

  wflat = wrappers.flatten()

  if 1 == 1:
    wmask = lon[wflat] < 4
    u0 = np.unique(wflat[wmask])
    x = ll.shape[0] + np.arange(u0.shape[0])
    map = np.arange(ll.shape[0])
    map[u0] = x
    wflat0 = map[wflat]
    wrappers0 = wflat0.reshape(-1,4)
    wrapped_ll0 = ll[u0] + [2*np.pi, 0, 0]

  if 1 == 1:
    wmask = lon[wflat] > 4
    u1 = np.unique(wflat[wmask])
    x = ll.shape[0] + wrapped_ll0.shape[0] + np.arange(u1.shape[0])
    map = np.arange(ll.shape[0])
    map[u1] = x
    wflat1 = map[wflat]
    wrappers1 = wflat1.reshape(-1,4)
    wrapped_ll1 = ll[u1] - [2*np.pi, 0, 0]

  cells = np.vstack((non_wrappers, wrappers0, wrappers1))
  points = np.vstack((ll, wrapped_ll0, wrapped_ll1))

  n_cells = cells.shape[0]
  cells = np.column_stack(([4]*n_cells, cells)).flatten()
  ctype = [9]*n_cells
  coff = 5*np.arange(cells.shape[0])

  v1 = vtkUnstructuredGrid()
  p = vtkPoints()
  p.SetData(ns.numpy_to_vtk(points))
  v1.SetPoints(p)

  carray = vtkCellArray()
  carray.SetCells(n_cells, ns.numpy_to_vtkIdTypeArray(cells))
  ctype = ns.numpy_to_vtk(np.array([VTK_QUAD]*n_cells).astype('u1'))
  coff = ns.numpy_to_vtkIdTypeArray(5*np.arange(n_cells))

  v1.SetCells(ctype, coff, carray)

  wv1 = dsa.WrapDataObject(v1)
  for name in v.PointData.keys():
    array = v.PointData[name]
    if len(array.shape) == 1:
      add = np.hstack((array, array[u0], array[u1]))
    else:
      add = np.vstack((array, array[u0], array[u1]))
    wv1.PointData.append(add, name)

  return wv1.VTKObject

def read_vtu(n):
  r = vtkXMLUnstructuredGridReader()
  r.SetFileName(vtus[n])
  r.Update()
  return dsa.WrapDataObject(r.GetOutput())

def vtu_timestamp(n):
  return float(vtus[n].rsplit('/')[-1].rsplit('_')[-1].split('.')[0])
  
def advec():
  v_next = read_vtu(0)
  vll_next = toLatLon(v_next)
  t_next = vtu_timestamp(0)

  samples = sample(v_next, MINLAT, NSEEDS)
  if X != -1:
    samples = [samples[X]]

  d = np.linalg.norm(samples, axis=1)
  samples = samples / d[:,np.newaxis]

  n_samples = len(samples)

  # vary initial termination - at least 2

  termination = 2+np.array((0.5 + np.random.random(n_samples)) * (TERMINATION_STEPS-2)).astype('i4')
  terminated = []

  number_of_points_along_each_line = 0
  first = True

  terminated_lines = []
  terminated_wind_magnitude = []
  terminated_timestamps = []

  for tstep in range(len(vtus))[1:]:
  
    vll_last = vll_next
    v_next = read_vtu(tstep)

    t_last = t_next
    t_next = vtu_timestamp(tstep)

    vll_next = toLatLon(v_next)
  
    for sub_timestep in range(0 if first else 1, N_SUBSTEPS+1):

      d_sub_timestep = sub_timestep / float(N_SUBSTEPS)

      this_timestep = int(t_last + d_sub_timestep*(t_next - t_last))

      if first:
        endpoints = samples
        first = False
        lines = [np.array([]).reshape(-1,3) for i in range(n_samples)]
        wind_magnitude = [np.array([]) for i in range(n_samples)]
        timestamps = [np.array([]) for i in range(n_samples)]
      else:
        termination = termination - 1
        tflags = termination < 0
        if np.sum(np.where(tflags, 1, 0)) > 0:
          # Transfer terminated lines and re-init the corresponding slots
          terminated_lines = terminated_lines + [lines[i] for i in range(len(lines)) if tflags[i]]
          terminated_wind_magnitude = terminated_wind_magnitude + [wind_magnitude[i] for i in range(len(lines)) if tflags[i]]
          terminated_timestamps = terminated_timestamps + [timestamps[i] for i in range(len(lines)) if tflags[i]]
          lines = [lines[i] if not tflags[i] else np.array([]).reshape(-1,3) for i in range(len(lines))]
          wind_magnitude = [wind_magnitude[i] if not tflags[i] else np.array([]) for i in range(len(lines))]
          timestamps = [timestamps[i] if not tflags[i] else np.array([]) for i in range(len(lines))]
          endpoints = np.array([s if t else e for t,s,e in zip(tflags, samples, endpoints)])
          termination = np.where(tflags, TERMINATION_STEPS, termination)

      ll_endpoints = xyz_to_latlon(endpoints)
  
      # Create flat VTU with current samples
  
      vs = vtkUnstructuredGrid()
      p = vtkPoints()
      p.SetData(ns.numpy_to_vtk(ll_endpoints))
      vs.SetPoints(p)
      carray = vtkCellArray()
      cells = np.column_stack(([1]*ll_endpoints.shape[0], np.arange(ll_endpoints.shape[0]))).flatten()
      carray.SetCells(ll_endpoints.shape[0], ns.numpy_to_vtkIdTypeArray(cells.astype('i8')))
      ctype = ns.numpy_to_vtk(np.array([VTK_VERTEX]*ll_endpoints.shape[0]).astype('u1'))
      coff = ns.numpy_to_vtkIdTypeArray(2*np.arange(ll_endpoints.shape[0]))
      vs.SetCells(ctype, coff, carray)

      sampler = vtkResampleWithDataSet()
      sampler.SetSourceData(vll_last)
      sampler.SetInputData(vs)
      sampler.Update()
      r = dsa.WrapDataObject(sampler.GetOutput())
      wind0 = r.PointData[WIND_ARRAY]

      sampler.SetSourceData(vll_next)
      sampler.Update()
      r = dsa.WrapDataObject(sampler.GetOutput())
      wind1 = r.PointData[WIND_ARRAY]
  
      wind = d_sub_timestep * wind1 + (1.0 - d_sub_timestep) * wind0
      wmag = np.sqrt(wind[:,0]*wind[:,0] + wind[:,1]*wind[:,1])

      ### Timesteps are 6 hourly, wind vector is m/s on Earth scale; we are working in meters
      wind = (wind * 21600) / N_SUBSTEPS     

      # new_endpoints = endpoints + wind[:,0] + wind[:,1]
      new_endpoints = endpoints + wind
      d = np.linalg.norm(new_endpoints, axis=1)
      new_endpoints = new_endpoints / d[:,np.newaxis]

      lines = [np.vstack((l, e)) for l,e in zip(lines, new_endpoints)]
      wind_magnitude = [np.concatenate((wm, [w])) for wm,w in zip(wind_magnitude, wmag)]
      timestamps = [np.concatenate((ts, [this_timestep])) for ts in timestamps]

      number_of_points_along_each_line = number_of_points_along_each_line + 1
      endpoints = new_endpoints

      terminated_lines = [t[1:] for t in terminated_lines]
      terminated_wind_magnitude = [t[1:] for t in terminated_wind_magnitude]
      terminated_timestamps = [t[1:] for t in terminated_timestamps]
  
    if tstep % WRITE_EVERY == 0:
  
      # Trim lines

      terminated_lines = [t for t in terminated_lines if t.shape[0] > 1]
      terminated_wind_magnitude = [t for t in terminated_wind_magnitude if t.shape[0] > 1]
      terminated_timestamps = [t for t in terminated_timestamps if t.shape[0] > 1]

      lines = [l[-MAX_NUMBER_OF_POINTS_ALONG_LINES:] for l in lines]
      wind_magnitude = [wm[-MAX_NUMBER_OF_POINTS_ALONG_LINES:] for wm in wind_magnitude]
      timestamps = [ts[-MAX_NUMBER_OF_POINTS_ALONG_LINES:] for ts in timestamps]

      pts = np.vstack(([l for l in lines if l.shape[0] > 1] + terminated_lines))
      wmags = np.concatenate(([w for w in wind_magnitude if w.shape[0] > 1] + terminated_wind_magnitude))
      tstamps = np.concatenate(([ts for ts in timestamps if ts.shape[0] > 1] + terminated_timestamps))

      line_offsets = []
      line_indices = []
      line_offset = 0
      index_offset = 0
  
      for l in [l for l in lines if l.shape[0] > 1] + terminated_lines:
        line_offsets.append(line_offset)
        line_offset = line_offset + (l.shape[0] + 1)
        line_indices.append([l.shape[0]] + list(range(index_offset, index_offset+l.shape[0])))
        index_offset = index_offset + l.shape[0]

      line_count = len(line_indices)
      line_indices = np.concatenate(line_indices)

      p = vtkPoints()
      p.SetData(ns.numpy_to_vtk(pts.astype('f4')))
      vs.SetPoints(p)
      carray = vtkCellArray()
      carray.SetCells(line_count, ns.numpy_to_vtkIdTypeArray(line_indices))
      ctype = ns.numpy_to_vtk(np.array([VTK_POLY_LINE]*line_count).astype('u1'))
      coff = ns.numpy_to_vtkIdTypeArray(np.array(line_offsets))
      vs.SetCells(ctype, coff, carray)
      cnum = ns.numpy_to_vtk(np.arange(line_count).astype('i4'))
      cnum.SetName("cellID")
      vs.GetCellData().AddArray(cnum)
      wmags = ns.numpy_to_vtk(wmags.astype('f4'))
      wmags.SetName("velocity")
      vs.GetPointData().AddArray(wmags)
      tstamps = ns.numpy_to_vtk(tstamps.astype('f4'))
      tstamps.SetName("timestamp")
      vs.GetPointData().AddArray(tstamps)
      wrtr = vtkXMLUnstructuredGridWriter()
      wrtr.SetFileName(FILENAME_TEMPLATE % (this_timestep))
      wrtr.SetInputData(vs)
      wrtr.Write()

      print('timestep', tstep-1, 'written')
  
    print("done")

def syntax():
  print("syntax: vtkpython advec.py [args] template")
  print("args:")
  print("    -N substeps        substeps in each timestep interval(10)")
  print("    -P points          max number of points in any streamline (100)")
  print("    -T term            max number of samples from any one seed (200)")
  print("    -S nseeds          appx number of seeds (500)")
  print("    -O template        output file template ('tt_%05d.vtu')")
  print("    -A variable        wind vector (WBOT)")
  print("    -X n               (testing seed n)")

  sys.exit()

if __name__ == '__main__':

  args = sys.argv[1:]
  while len(args) > 0:
    if args[0] == '-N':
      N_SUBSTEPS = int(args[1])
      args = args[2:]
    elif args[0] == '-P': 
      MAX_NUMBER_OF_POINTS_ALONG_LINES = int(args[1])
      args = args[2:]
    elif args[0] == '-T': 
      TERMINATION_STEPS = int(args[1])
      args = args[2:]
    elif args[0] == '-S': 
      NSEEDS = int(args[1])
      args = args[2:]
    elif args[0] == '-O': 
      FILENAME_TEMPLATE = args[1]
      args = args[2:]
    elif args[0] == '-A':
      WIND_ARRAY = args[1]
      args = args[2:]
    elif args[0] == '-R':
      rs = np.random.seed(int(args[1]))
      args = args[2:]
    elif args[0] == '-X':
      X = int(args[1])
      args = args[2:]
    elif template == "":
      template = args[0]
      args = args[1:]
    else:
      syntax()

  if template == '':
    syntax()

  vtus = sorted(glob(template))

  if '/' in FILENAME_TEMPLATE:
    dname = FILENAME_TEMPLATE.rsplit('/', 1)[0]
    if not os.path.isdir(dname):
      os.mkdir(dname)

  advec()
