Name = 'OrientedGlyph'
Label = 'Oriented Glyph'
Help = 'Place a glyph at each point using two vectors for orientation'

NumberOfInputs = 2
InputDataType = 'vtkUnstructuredGrid'
OutputDataType = 'vtkUnstructuredGrid'
ExtraXml = ''

Properties = dict(
forward = 'velocity',
up = 'Normals',
scale = 1.0,
xscale = 1.0,
yscale = 1.0,
zscale = 1.0
)

def RequestData():
  print 'here'

  import numpy as np
  from vtk.numpy_interface import dataset_adapter as dsa
  from vtk.util import numpy_support as ns

  ipoints = inputs[0]
  number_of_glyphs = ipoints.GetNumberOfPoints()

  print 'nG:', number_of_glyphs

  glyph = inputs[1]
  glyph_points = [scale*xscale, scale*yscale, scale*zscale] * glyph.Points
  glyph_polygons = glyph.Cells.reshape((-1,4))

  print 'nP:', inputs[1].GetNumberOfPoints()

  points_per_glyph = glyph_points.shape[0]
  triangles_per_glyph = glyph_polygons.shape[0]

  if forward not in ipoints.PointData.keys():
    print 'can\'t find forward array'
    return

  U = ipoints.PointData[forward]

  if up not in ipoints.PointData.keys():
    print 'can\'t find up array'
    return

  V = ipoints.PointData[up]

  W = dsa.VTKArray(np.cross(U, V))

  l = np.linalg.norm(U, axis=1)
  U = U / np.where(l == 0, 1.0, l)
  l = np.linalg.norm(U, axis=1)
  V = V / np.where(l == 0, 1.0, l)
  l = np.linalg.norm(W, axis=1)
  W = W / np.where(l == 0, 1.0, l)

  P = ipoints.Points

  p = P[0]
  u = U[0]
  v = V[0]
  w = W[0]

  opoints = []
  opolys = []

  iDataArrays = []
  oDataArrays = []
  print '111'
  for i,k in enumerate(ipoints.PointData.keys()):
    print i, k
    iDataArrays.append(ipoints.PointData[k])
    oDataArrays.append([])
  print '222'

  for i,p,u,v,w in zip(range(len(P)), P, U, V, W):
    new_points = p + glyph_points[:,0][:,np.newaxis]*u + glyph_points[:,1][:,np.newaxis]*v + glyph_points[:,2][:,np.newaxis]*w 
    new_polys  = glyph_polygons + [0, i*points_per_glyph, i*points_per_glyph, i*points_per_glyph]
    opoints.append(new_points)
    opolys.append(new_polys)
    for id,od in zip(iDataArrays, oDataArrays):
       od.append([id[i]]*points_per_glyph)

  opoints = dsa.numpyTovtkDataArray(np.vstack(opoints))
  opolys = ns.numpy_to_vtkIdTypeArray(np.ascontiguousarray(np.vstack(opolys).flatten()))

  ids = [np.array([i]*points_per_glyph) for i in range(len(P))]
  ids = dsa.numpyTovtkDataArray(np.vstack(ids).flatten(), name='ID')

  print 'AAAAAAA'
  oug = self.GetUnstructuredGridOutput()

  pts = vtk.vtkPoints()
  pts.SetData(opoints)
  oug.SetPoints(pts)

  ct = dsa.numpyTovtkDataArray(np.array([vtk.VTK_TRIANGLE]*triangles_per_glyph*number_of_glyphs).astype('u1'))
  co = dsa.numpy_support.numpy_to_vtkIdTypeArray(np.arange(0, 4*triangles_per_glyph*number_of_glyphs, 4))

  ca = vtk.vtkCellArray()
  ca.SetCells(number_of_glyphs*triangles_per_glyph, opolys)

  print 'BBBBBBBB'
  oug.SetCells(ct, co, ca)
  oug.GetPointData().AddArray(ids)
  print 'CCCCCCC'

  for i,n in enumerate(ipoints.PointData.keys()):
    oug.GetPointData().AddArray(dsa.numpyTovtkDataArray(np.vstack(oDataArrays[i]).flatten(), name=n))

  print 'DDDDDDD'
  return


