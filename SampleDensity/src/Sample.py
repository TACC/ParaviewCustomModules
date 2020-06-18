import sys, os
import vtk
from vtk.numpy_interface import dataset_adapter as dsa

import ctypes as C
import numpy as np

if len(sys.argv) > 1:
  nSamples = int(sys.argv[1])
else:
  nSamples = 100;

ccode_so = os.environ['HOME'] + '/python/_sample.so'
if not os.path.isfile(ccode_so):
  print('can\'t find so (', ccode_so, ')')
  sys.exit()

print('loading ', ccode_so)

ccode = C.CDLL(ccode_so)
ccode.SampleTetrahedra.argtypes = [C.c_int, C.c_int, C.c_int, 
                         np.ctypeslib.ndpointer(C.c_float, flags="C_CONTIGUOUS"),
                         np.ctypeslib.ndpointer(C.c_float, flags="C_CONTIGUOUS"),
                         np.ctypeslib.ndpointer(C.c_int, flags="C_CONTIGUOUS")]

ccode.SampleVTI.argtypes = [C.c_int, 
                            np.ctypeslib.ndpointer(C.c_int, flags="C_CONTIGUOUS"),
                            np.ctypeslib.ndpointer(C.c_float, flags="C_CONTIGUOUS"),
                            np.ctypeslib.ndpointer(C.c_float, flags="C_CONTIGUOUS"),
                            np.ctypeslib.ndpointer(C.c_float, flags="C_CONTIGUOUS")]


ccode.GetNumberOfSamples.restype = C.c_int
ccode.GetSamples.restype = C.c_void_p

def SampleTetrahedra(vtu, pdf, n):
  points = np.ascontiguousarray(vtu.Points).astype('f4')
  pdf    = np.ascontiguousarray(vtu.CellData[pdf]).astype('f4')
  tets   = np.ascontiguousarray(vtu.Cells).astype('i4')
  ccode.SampleTetrahedra(n, vtu.GetNumberOfPoints(), vtu.GetNumberOfCells(), points, pdf, tets)
  n = ccode.GetNumberOfSamples()
  s = np.ctypeslib.as_array(C.cast(ccode.GetSamples(), C.POINTER(C.c_float)),shape=(n, 3))
  samples = dsa.WrapDataObject(vtk.vtkUnstructuredGrid())
  co = dsa.numpy_support.numpy_to_vtkIdTypeArray(np.arange(n).astype('i8')*2)
  ca = vtk.vtkCellArray()
  ca.SetCells(n, dsa.numpy_support.numpy_to_vtkIdTypeArray(np.column_stack(([1]*n, range(n))).reshape((2*n,))))
  ct = dsa.numpyTovtkDataArray(np.array([vtk.VTK_VERTEX]*n).astype('u1'))
  samples.VTKObject.SetCells(ct, co, ca)
  samples.Points = dsa.numpy_support.numpy_to_vtk(s)
  return samples

def SampleVTI(vti, pdf, n):
  e = vti.VTKObject.GetExtent()
  o = vti.VTKObject.GetOrigin()
  s = vti.VTKObject.GetSpacing()
  dimensions = np.array([(e[2*i+1] - e[2*i])+1 for i in range(3)]).astype('i4')
  origin     = np.array([o[i] - e[2*i]*s[i] for i in range(3)]).astype('f4')
  spacing    = np.array(s).astype('f4')
  pdf        = np.ascontiguousarray(vti.PointData[pdf]).astype('f4')
  ccode.SampleVTI(n, dimensions, origin, spacing, pdf);
  n = ccode.GetNumberOfSamples()
  s = np.ctypeslib.as_array(C.cast(ccode.GetSamples(), C.POINTER(C.c_float)),shape=(n, 3))
  samples = dsa.WrapDataObject(vtk.vtkUnstructuredGrid())
  co = dsa.numpy_support.numpy_to_vtkIdTypeArray(np.arange(n).astype('i8')*2)
  ca = vtk.vtkCellArray()
  ca.SetCells(n, dsa.numpy_support.numpy_to_vtkIdTypeArray(np.column_stack(([1]*n, range(n))).reshape((2*n,)).astype('i8')))
  ct = dsa.numpyTovtkDataArray(np.array([vtk.VTK_VERTEX]*n).astype('u1'))
  samples.VTKObject.SetCells(ct, co, ca)
  samples.Points = dsa.numpy_support.numpy_to_vtk(s)
  return samples
  
rdr = vtk.vtkXMLImageDataReader()
rdr.SetFileName('../noise.vti')
rdr.Update()

igrid = rdr.GetOutput()

samples = SampleVTI(dsa.WrapDataObject(igrid), 'noise', nSamples)

wrtr = vtk.vtkXMLUnstructuredGridWriter()
wrtr.SetFileName('vti-samples.vtu')
wrtr.SetInputData(samples.VTKObject)
wrtr.Write()

e = igrid.GetExtent()
o = igrid.GetOrigin()
s = igrid.GetSpacing()
n = igrid.GetDimensions()

x = dsa.numpy_support.numpy_to_vtk(np.array([o[0] - e[0]*s[0] + i*s[0] for i in range(e[1]+1)]).astype('f8'))
y = dsa.numpy_support.numpy_to_vtk(np.array([o[1] - e[2]*s[1] + i*s[1] for i in range(e[3]+1)]).astype('f8'))
z = dsa.numpy_support.numpy_to_vtk(np.array([o[2] - e[4]*s[2] + i*s[2] for i in range(e[5]+1)]).astype('f8'))

rgrid = vtk.vtkRectilinearGrid()
rgrid.SetDimensions(n)
rgrid.SetXCoordinates(x)
rgrid.SetYCoordinates(y)
rgrid.SetZCoordinates(z)
rgrid.GetPointData().ShallowCopy(igrid.GetPointData())

to_tets = vtk.vtkRectilinearGridToTetrahedra()
to_tets.SetInputData(rgrid)
to_tets.Update()

tets = to_tets.GetOutput()
tets.GetPointData().ShallowCopy(rgrid.GetPointData())

w = dsa.WrapDataObject(tets)

p_to_c = vtk.vtkPointDataToCellData()
p_to_c.SetInputConnection(to_tets.GetOutputPort())
p_to_c.Update()
ctets = p_to_c.GetOutput()

vtu = dsa.WrapDataObject(p_to_c.GetOutput())

samples = SampleTetrahedra(vtu, 'noise', nSamples)

wrtr = vtk.vtkXMLUnstructuredGridWriter()
wrtr.SetFileName('tet-samples.vtu')
wrtr.SetInputData(samples.VTKObject)
wrtr.Write()
