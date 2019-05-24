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

Name = 'CurvatureStreamlines'
Label = 'Curvature Streamlines'
Help = 'Create curvature vector for streamlines'

NumberOfInputs = 1
InputDataType = 'vtkUnstructuredGrid'
OutputDataType = 'vtkUnstructuredGrid'
ExtraXml = ''

Properties = dict(
  stepsize = 100000
)

def RequestData():
  import numpy as np
  from vtk.numpy_interface import dataset_adapter as dsa
  from math import ceil

  isl = self.GetUnstructuredGridInput()
  osl = self.GetUnstructuredGridOutput()

  nisl = dsa.WrapDataObject(isl)
  nosl = dsa.WrapDataObject(osl)

  osl.SetPoints(isl.GetPoints())
  for i in range(isl.GetPointData().GetNumberOfArrays()):
    osl.GetPointData().AddArray(isl.GetPointData().GetArray(i))

  cell_locations = nisl.GetCellLocations()
  
  vPerL  = nisl.Cells[cell_locations]
  starts = cell_locations + 1
  ends   = starts + vPerL

  lines = []
  for s,e in zip(starts, ends):
    lines.append(nisl.Cells[s:e])

  if 'arclen' in nisl.PointData.keys():
    arclen = nisl.PointData['arclen']
  else:
    arclen = np.zeros(nisl.GetNumberOfPoints())
    nosl.PointData.append(arclen, 'arclen')

    for line in lines:
      points = nisl.Points[line]
      seglen = np.linalg.norm(points[1:] - points[:-1], axis=1)
      arclen[line] = np.concatenate(([0], np.cumsum(seglen)))

  curvature = np.zeros(len(nisl.Points)*3).reshape((-1, 3))
  binormal = np.zeros(len(nisl.Points)*3).reshape((-1, 3))

  for indx,line in enumerate(lines):
    line_points = nisl.Points[line]
    line_arclen = arclen[line]
    llen = line_arclen[-1]
    nsteps = ceil(llen / stepsize)
    if nsteps < 3:
      nsteps = 3
    samples = np.linspace(0.0, llen, nsteps)
    sampled_points = np.column_stack([np.interp(samples, line_arclen, line_points[:,i]) for i in range(3)])
    vectors = sampled_points[1:] - sampled_points[:-1]
    d = np.linalg.norm(vectors, axis=1)
    d = np.where(d == 0.0, 1.0, d)
    vectors = vectors / d[:,np.newaxis]
    crosses = np.cross(vectors[1:], vectors[:-1], axis=1)
    crosses = np.concatenate(([crosses[0]], crosses, [crosses[-1]]))
    crosses = np.column_stack([np.interp(line_arclen, samples, crosses[:,i]) for i in range(3)])
    # d = np.linalg.norm(crosses, axis=1)
    # d = np.where(d == 0.0, 1.0, d)
    # crosses = crosses / d[:,np.newaxis]
    binorm = 0.5 * (vectors[1:] - vectors[:-1])
    print indx, len(vectors), vectors[:10]
    # d = np.linalg.norm(binorm, axis=1)
    # d = np.where(d == 0.0, 1.0, d)
    # binorm = binorm / d[:,np.newaxis]
    print indx, len(binorm), binorm[:10]
    binorm = np.concatenate(([binorm[0]], binorm, [binorm[-1]]))
    binorm = np.column_stack([np.interp(line_arclen, samples, binorm[:,i]) for i in range(3)])
    binormal[line] = binorm
    curvature[line] = crosses

  nosl.PointData.append(curvature, 'curvature')
  nosl.PointData.append(binormal, 'binormal')

