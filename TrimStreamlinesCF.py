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

Name = 'TrimStreamlines'
Label = 'Trim Streamlines'
Help = 'Trim streamlines based to a range of IntegrationTime'

NumberOfInputs = 1
InputDataType = 'vtkUnstructuredGrid'
OutputDataType = 'vtkUnstructuredGrid'
ExtraXml = ''

Properties = dict(
  start = -1,
  end   = -1
)

def RequestData():
  import numpy as np
  from vtk.numpy_interface import dataset_adapter as dsa

  sl = self.GetPolyDataInput()
  nsl = dsa.WrapDataObject(sl)

  itime = nsl.PointData['IntegrationTime']

  nv = nsl.Cells[nsl.CellLocations]    # number of verts in each line
  ns = nsl.CellLocations + 1           # index of first vertex in each line
  ne = ns + nv                         # index one past the last vertex
  olines = [nsl.Cells[i:j] for i,j in zip(ns, ne)]
  nlines = []

  iarrays = {'points': nsl.Points}    # initialize source arrays with input points
  oarrays = {'points': []}
  for n in nsl.PointData.keys():
    iarrays[n] = nsl.PointData[n]     # add input point data arrays to source arrays
    oarrays[n] = []                   # add empty destination arrays

  knt = 0
  for line in olines:
    if start != -1: line = [l for l in line if itime[l] > start]
    if end   != -1: line = [l for l in line if itime[l] < end]
    for n in iarrays.keys(): oarrays[n].append(iarrays[n][line])
    nlines.append(range(knt, knt+len(line)))
    knt = knt + len(line)

  tsl = vtk.vtkUnstructuredGrid()

  line_lengths = [len(l)+1 for l in nlines]  # number of vertices + 1 for count

  ct = dsa.numpyTovtkDataArray(np.array([vtk.VTK_POLY_LINE]*len(nlines)).astype('u1'))
  co = dsa.numpy_support.numpy_to_vtkIdTypeArray(np.hstack(([0], np.cumsum(line_lengths)[:-1])))

  ca = vtk.vtkCellArray()
  for l in nlines:
    ca.InsertNextCell(len(l), l)

  tsl.SetCells(ct, co, ca)

  cid = dsa.numpyTovtkDataArray(np.arange(len(nlines)).astype('i4'))
  cid.SetName('cell id')
  tsl.GetCellData().AddArray(cid)

  ptsa = np.concatenate(oarrays['points']).reshape((-1, 3)).astype('f4')

  op = vtk.vtkPoints()
  op.SetNumberOfPoints(ptsa.shape[0])

  for i, p in enumerate(ptsa):
    op.InsertPoint(i, p[0], p[1], p[2])

  tsl.SetPoints(op)

  for n in oarrays:
    if n != 'points':
      a = dsa.numpyTovtkDataArray(np.concatenate(oarrays[n]))
      a.SetName(n)
      tsl.GetPointData().AddArray(a)

  self.GetUnstructuredGridOutput().ShallowCopy(tsl)


