import vtk, os, sys, glob
from mpi4py import MPI
from vtk.numpy_interface import dataset_adapter as dsa
import numpy as np

vectors = {
  'W300': ['U300', 'V300'],
  'W500': ['U500', 'V500'],
  'W850': ['U850', 'V850'],
  'WBOT': ['UBOT', 'VBOT']
}

def norm(v):
  x = v[:,0]
  y = v[:,1]
  z = v[:,2]
  d = 1.0 / np.sqrt(x*x + y*y + z*z)
  return v * d[:,np.newaxis]

def uv(pyvtu):
  W = norm(pyvtu.Points)
  U = np.column_stack((-W[:,1], W[:,0], [0]*W.shape[0]))
  U = norm(U)
  V = np.cross(W, U)
  return(U,V)

def wind(vtu):
  c2p = vtk.vtkCellDataToPointData()
  c2p.SetInputData(vtu)
  c2p.Update()
  pyvtu = dsa.WrapDataObject(c2p.GetOutput())
  u,v = uv(pyvtu)
  for w in vectors:
    pyvtu.GetPointData().append(u*pyvtu.PointData[vectors[w][0]] + v*pyvtu.PointData[vectors[w][1]], w)
  return pyvtu

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

rdr = vtk.vtkXMLUnstructuredGridReader()
wrtr = vtk.vtkXMLUnstructuredGridWriter()

for i in glob.glob('atm/atm_*vtu')[rank::size]:
  oname = 'wind/%s\n' % i.split('/',1)[1]
  rdr.SetFileName(i)
  rdr.Update()
  w = wind(rdr.GetOutput())
  wrtr.SetInputData(w.VTKObject)
  wrtr.SetFileName(oname)
  wrtr.Write()
  print(oname, 'done')
