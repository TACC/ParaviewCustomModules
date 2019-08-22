## ========================================================================== ##
## Copyright (c) 2019 The University of Texas at Austin.                      ##
## All rights reserved.                                                       ##
##                                                                            ##
## Licensed under the Apache License, Version 2.0 (the "License");            ##
## you may not use this file except in compliance with the License.           ##
## A copy of the License is included with this software in the file LICENSE.  ##
## If your copy does not contain the License, you may obtain a copy of the    ##
## License at:                                                                ##
##                                                                            ##
##     https://www.apache.org/licenses/LICENSE-2.0                            ##
##                                                                            ##
## Unless required by applicable law or agreed to in writing, software        ##
## distributed under the License is distributed on an "AS IS" BASIS, WITHOUT  ##
## WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.           ##
## See the License for the specific language governing permissions and        ##
## limitations under the License.                                             ##
##                                                                            ##
## ========================================================================== ##

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
zscale = 1.0,
dbg = 999
)

def RequestData():
  import numpy as np
  from time import time
  from vtk.numpy_interface import dataset_adapter as dsa
  from vtk.util import numpy_support as ns

  ipoints = inputs[0]
  number_of_glyphs = ipoints.GetNumberOfPoints()

  glyph = inputs[1]

  glyph_points     = [scale*xscale, scale*yscale, scale*zscale] * glyph.Points
  points_per_glyph = glyph_points.shape[0]

  cells_per_glyph  = len(glyph.CellTypes)

  if forward not in ipoints.PointData.keys():
    print('can\'t find forward array')
    return

  U = ipoints.PointData[forward]

  if up not in ipoints.PointData.keys():
    print('can\'t find up array')
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

  xpts = glyph_points[:,0][:,np.newaxis]
  ypts = glyph_points[:,1][:,np.newaxis]
  zpts = glyph_points[:,2][:,np.newaxis]

  opoints = []
  for i,p,u,v,w in zip(range(len(P)), P, U, V, W):
    opoints.append(p + xpts*u + ypts*v + zpts*w)

  opolys = [glyph.Cells.reshape(-1,4)]

  ijk = glyph.Cells.reshape((-1,4))[:,1:4]

  for i in range(1, len(P)):
    nijk = np.column_stack(([3]*ijk.shape[0], ijk + i*points_per_glyph))
    opolys.append(nijk)

  opoints = dsa.numpyTovtkDataArray(np.vstack(opoints))

  ids = [np.array([i]*points_per_glyph) for i in range(len(P))]
  ids = dsa.numpyTovtkDataArray(np.vstack(ids).flatten(), name='ID')

  oug = vtk.vtkUnstructuredGrid()

  pts = vtk.vtkPoints()
  pts.SetData(opoints)
  oug.SetPoints(pts)

  ct = np.hstack([glyph.CellTypes for i in range(number_of_glyphs)])
  co = np.hstack([glyph.CellLocations + i*len(glyph.Cells) for i in range(number_of_glyphs)])
  opolys = np.hstack(opolys).astype('i8')

  ct = dsa.numpyTovtkDataArray(ct)
  co = dsa.numpy_support.numpy_to_vtkIdTypeArray(co)
  opolys = ns.numpy_to_vtkIdTypeArray(opolys)

  ca = vtk.vtkCellArray()
  ca.SetCells(number_of_glyphs*cells_per_glyph, opolys)

  oug.SetCells(ct, co, ca)
  oug.GetPointData().AddArray(ids)

  oug.GetPointData().AddArray(dsa.numpyTovtkDataArray(np.vstack([glyph.PointData['Normals'] for i in range(number_of_glyphs)]), name='Normals'))

  if 'Texture Coordinates' in glyph.PointData.keys():  
    a = np.vstack([glyph.PointData['Texture Coordinates'] for	 i in range(len(ipoints.Points))])
    oug.GetPointData().SetTCoords(dsa.numpyTovtkDataArray(a))
	
  for n in ipoints.PointData.keys():
    if n != 'Normals':
      a = [[ipoints.PointData[n][i]]*points_per_glyph for i in range(number_of_glyphs)]
      oug.GetPointData().AddArray(dsa.numpyTovtkDataArray(np.concatenate(a), name=n))

  self.GetUnstructuredGridOutput().Initialize()
  self.GetUnstructuredGridOutput().ShallowCopy(oug)

  return

if __name__ == '__main__':

  class s:
    def __init__(self, o):
      from vtk import vtkUnstructuredGrid
      self.outpt = vtkUnstructuredGrid()
      self.outpt.DeepCopy(o.VTKObject)
    def GetUnstructuredGridOutput(self):
      return self.outpt

  from vtk import *
  from vtk.numpy_interface import dataset_adapter as dsa

  inputs = []

  prdr = vtkXMLUnstructuredGridReader()
  prdr.SetFileName('/Users/gda/Glyphs/points.vtu')
  prdr.Update()
  inputs.append(dsa.WrapDataObject(prdr.GetOutput()))

  grdr = vtkXMLUnstructuredGridReader()
  grdr.SetFileName('/Users/gda/Glyphs/bird2-small.vtu')
  grdr.Update()
  inputs.append(dsa.WrapDataObject(grdr.GetOutput()))

  self = s(inputs[1])
  locals().update(Properties)
  eval('RequestData()', globals(), locals())

