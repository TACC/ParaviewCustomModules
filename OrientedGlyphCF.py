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
  import numpy as np
  from vtk.numpy_interface import dataset_adapter as dsa
  from vtk.util import numpy_support as ns

  ipoints = inputs[0]
  number_of_glyphs = ipoints.GetNumberOfPoints()

  glyph = inputs[1]

  glyph_points     = [scale*xscale, scale*yscale, scale*zscale] * glyph.Points
  points_per_glyph = glyph_points.shape[0]

  cells_per_glyph  = len(glyph.CellTypes)

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
  for i,p,u,v,w in zip(range(len(P)), P, U, V, W):
    opoints.append(p + glyph_points[:,0][:,np.newaxis]*u + glyph_points[:,1][:,np.newaxis]*v + glyph_points[:,2][:,np.newaxis]*w)

  print "Cells", glyph.Cells

  opolys = [glyph.Cells]
  for i in range(1, len(P)):
    o = np.zeros(len(glyph.Cells))
    k = 0
    for j in range(len(glyph.Cells)):
      if k == 0:
        k = glyph.Cells[j]
        o[j] = k
      else:
        k = k - 1
        o[j] = glyph.Cells[j] + i*points_per_glyph
    opolys.append(o)

  opoints = dsa.numpyTovtkDataArray(np.vstack(opoints))

  print '44', opolys

  ids = [np.array([i]*points_per_glyph) for i in range(len(P))]
  ids = dsa.numpyTovtkDataArray(np.vstack(ids).flatten(), name='ID')

  print 'AAAAAAA'

  oug = self.GetUnstructuredGridOutput()

  pts = vtk.vtkPoints()
  pts.SetData(opoints)
  oug.SetPoints(pts)

  ct = np.hstack([glyph.CellTypes for i in range(number_of_glyphs)])
  co = np.hstack([glyph.CellLocations + i*len(glyph.Cells) for i in range(number_of_glyphs)])
  opolys = np.hstack(opolys).astype('i8')

  print 'CT', ct
  print 'CO', co
  print 'OP', opolys

  ct = dsa.numpyTovtkDataArray(ct)
  co = dsa.numpy_support.numpy_to_vtkIdTypeArray(co)
  opolys = ns.numpy_to_vtkIdTypeArray(opolys)

  ca = vtk.vtkCellArray()
  ca.SetCells(number_of_glyphs*cells_per_glyph, opolys)

  print 'BBBBBBBB'
  oug.SetCells(ct, co, ca)
  oug.GetPointData().AddArray(ids)
  print 'CCCCCCC'

  for i,n in enumerate(ipoints.PointData.keys()):
    oug.GetPointData().AddArray(dsa.numpyTovtkDataArray(np.vstack([ipoints.PointData[n] for i in range(number_of_glyphs)]), name=n))

  print 'DDDDDDD'
  return


