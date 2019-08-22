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

Name = 'SampleStreamlinesInTime'
Label = 'Sample Streamlines In Time'
Help = 'One sample at of each pathline at t (0 to 1) of the way along'

NumberOfInputs = 1
InputDataType = 'vtkUnstructuredGrid'
OutputDataType = 'vtkUnstructuredGrid'
ExtraXml = ''

Properties = dict(
t = 0.5,
nt = 1
)

def RequestData():
  import numpy as np
  from vtk.numpy_interface import dataset_adapter as dsa
  from math import ceil, floor

  sl = self.GetUnstructuredGridInput()    # The streamlines
  nsl = dsa.WrapDataObject(sl)            # wrap with a Python interface

  itime = nsl.PointData['IntegrationTime']  # the time component of the streamlines

  nv = nsl.Cells[nsl.CellLocations]    # number of verts in each line
  ns = nsl.CellLocations + 1           # index of first vertex in each line
  ne = ns + nv                         # index one past the last vertex

  lines = [nsl.Cells[i:j] for i,j in zip(ns, ne)] # divide into distinct lines
                                                  # lines[i] is a list of the ids of the
                                                  # vertices that comprise each streamline

  # Get length (in integration time) of longest line.   Note - this assumes that the 
  # forward- and backward- integrations are combined into one streamline (see JoinStreamlines)

  mint = itime[lines[0][0]]
  maxt = itime[lines[0][-1]]
  maxlen = itime[lines[0][-1]] - itime[lines[0][0]]
  if len(lines) > 1:
    for i in range(1, len(lines)):
      mt = itime[lines[i][-1]] - itime[lines[i][0]]
      if mt > maxlen:
        maxlen = mt;
      if mint > itime[lines[i][0]]: mint = itime[lines[i][0]]
      if maxt < itime[lines[i][-1]]: maxt = itime[lines[i][-1]]

  # dt is the distance between samples in integration time - nt samples distributed along longest line

  dt = (maxt - mint) / nt

  # destination arrays for streamline points and any point-dependent data - eg. orientation data

  iarrays = {'points': nsl.Points}    # initialize source arrays with input points
  oarrays = {'points': []}            # initialize destination arrays with (empty) points

  for n in nsl.PointData.keys():      
    iarrays[n] = nsl.PointData[n]     # add input point data arrays to source arrays
    oarrays[n] = []                   # add empty destination arrays 

  # for each sample time...

  for it in range(nt):

    sample_t = mint + (it + t) * dt     # the point in time to interpolate at

    for i,line in enumerate(lines):     # for each input line...

      # if this sample time is in the range for the current line...

      if sample_t >= itime[line[0]] and sample_t <= itime[line[-1]]:

        # index of first elt greater than sample_x (or 0, in which case we use the last)

        interval_end = np.argmax(itime[line] > sample_t)      # linear search?
        if interval_end == 0: interval_end = len(line) - 1

        # get indices of points and point-dependent data at either end of the interval
        endi = line[interval_end]
        starti = line[interval_end - 1]

        # interpolant value in interval
        d = (sample_t - itime[starti]) / (itime[endi] - itime[starti])  # interpolant in interval

        for n in iarrays:                 #   for each array we are interpolating...
          ia = iarrays[n]                 #     input array
          sv = ia[starti]                 #     start values
          ev = ia[endi]                   #     end values
          v  = sv + d*(ev - sv)           #     interpolation
          oarrays[n].append(v)

  # create an output vtkUnstructured data with the interpolated points and data

  ptsa = np.concatenate(oarrays['points']).reshape((-1, 3)).astype('f4')
  oug = vtk.vtkUnstructuredGrid()

  op = vtk.vtkPoints()
  op.SetNumberOfPoints(ptsa.shape[0])

  for i, p in enumerate(ptsa):
    op.InsertPoint(i, p[0], p[1], p[2])


  oug.SetPoints(op)

  for n in oarrays:
    if n != 'points':
      if oarrays[n][0].__class__ == dsa.VTKArray:
        ncomp = len(oarrays[n][0])
        a = dsa.numpyTovtkDataArray(np.concatenate(oarrays[n]).reshape((-1, ncomp)))
      else:
        a = dsa.numpyTovtkDataArray(oarrays[n])
      a.SetName(n)
      oug.GetPointData().AddArray(a)

  ct = dsa.numpyTovtkDataArray(np.array([vtk.VTK_VERTEX]*oug.GetNumberOfPoints()).astype('u1'))
  co = dsa.numpy_support.numpy_to_vtkIdTypeArray(np.array(range(0, 2*oug.GetNumberOfPoints(), 2)))

  ca = vtk.vtkCellArray()
  for i in range(oug.GetNumberOfPoints()):
    ca.InsertNextCell(1, [i])

  oug.SetCells(ct, co, ca)

  self.GetUnstructuredGridOutput().ShallowCopy(oug)
