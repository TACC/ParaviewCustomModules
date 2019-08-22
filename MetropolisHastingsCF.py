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

Name = 'MetropolisHastingsSampler'
Label = 'Metropolis Hastings Sampler'
Help = 'The MH Sampler uses a Markov Chain algorithm to produces samples of a volumetric dataset in the vicinity of a particular chosen data value.   The probability density is derived from an underlying scalar field as the difference of the scalar field from the target value, scaled to range from 0 to 1 in an interval given by the \'kernelwidth\' parameter, raised to a power given by the \'power\' parameter.  \'allabove\' and \'allbelow\' can be set to assign a probability of 1 in the data space above (or below) the target value.'

NumberOfInputs = 2
InputDataType = ['vtkPolyData', 'vtkUnstructuredGrid', 'vtkImageData']
OutputDataType = 'vtkUnstructuredGrid'
ExtraXml = ''

Properties = dict(
  arrayName = "RTData",
  starts = 100,
  allabove = False,
  allbelow = False,
  kernelwidth = 10.0,
  nmisses = 1000,
  nsamples = 1000,
  power = 10,
  target = 99999.0,
  sscale = 100,
  asscale = -1
)


def RequestData():
  import vtk
  import random
  import numpy as np
  from vtk.numpy_interface import dataset_adapter as dsa

  class Interpolator:
  
    def __init__(self, dset):
      self.dset = dset.VTKObject
      self.xyz = [-10000, -20000, -30000]
      self.pids = vtk.reference([0]*10)
      self.nverts = -1
      self.pc = [0]*3
      self.wts = [0]*10
      self.gc = vtk.vtkGenericCell()
      self.sid = 2
      if self.dset.IsA('vtkUnstructuredGrid'):
        self.locator = vtk.vtkCellTreeLocator()
        self.locator.SetDataSet(dset.VTKObject)
        self.locator.BuildLocator()
        self.is_vtu = True
      else:
        self.is_vtu = False
  
    def Locate(self, xyz):
      if self.is_vtu:
        cid = self.locator.FindCell(xyz, 0.0, self.gc, self.pc, self.wts)
        if cid < 0 or min(self.wts[:4]) < 0 or max(self.wts[:4]) > 1:
          self.xyz = []
          return False
        idl = vtk.vtkIdList()
        self.dset.GetCellPoints(cid, idl)
        self.ids = [idl.GetId(i) for i in range(idl.GetNumberOfIds())]
      else:
        vox = self.dset.FindAndGetCell(xyz, None, 0, 0.0, vtk.reference(self.sid), self.pc, self.wts)
        if vox == None:
          self.xyz = []
          return None
        self.ids = [vox.GetPointId(i) for i in range(vox.GetNumberOfPoints())]
      self.xyz = xyz
      return True
  
    def Interpolate(self, xyz, a):
      if list(xyz) != list(self.xyz):
        if not self.Locate(xyz):
          return None
      if len(a.shape) == 1:
        return sum(self.wts[i]*a[self.ids[i]] for i in range(len(self.ids)))
      else:
        return [sum(self.wts[i]*a[self.ids[i]][j] for i in range(len(self.ids))) for j in range(a.shape[1])]
  
  class Samples:
  
    def __init__(self, dset):
      self.points = []
      self.vars = []
      self.I = []
      self.Q = []
      for i in dset.PointData.keys():
        self.vars.append([i, dset.PointData[i], []])
  
    def num(self):
      return len(self.points)
  
    def add(self, I, p, q, i):
      self.points.append(p)
      self.Q.append(q)
      self.I.append(q)
      for i in self.vars:
        value = I.Interpolate(p, i[1])
        if value == None:
          print 'oops'
        i[2].append(value)
  
    def stuff_vtu(self, outpt):
      outpt.SetPoints(dsa.VTKArray(np.array(self.points).astype('f4')))
      outpt.PointData.append(dsa.VTKArray(np.array(self.Q).astype('f4')), 'Q')
      outpt.PointData.append(dsa.VTKArray(np.array(self.I).astype('f4')), 'I')
      ct = dsa.numpyTovtkDataArray(np.array([vtk.VTK_VERTEX]*outpt.GetNumberOfPoints()).astype('u1'))
      co = dsa.numpy_support.numpy_to_vtkIdTypeArray(np.array(range(0, 2*outpt.GetNumberOfPoints(), 2)))
      ca = vtk.vtkCellArray()
      for i in range(outpt.GetNumberOfPoints()):
        ca.InsertNextCell(1, [i])
      outpt.VTKObject.SetCells(ct, co, ca)
      for v in self.vars:
        outpt.PointData.append(dsa.VTKArray(np.array(v[2]).astype('f4')), v[0])

  print 'arrayName', arrayName
  print 'target', target
  print 'foo', kernelwidth, asscale

  print 'inputs[0]', inputs[0].PointData.keys()
  print 'inputs[1]', inputs[1].PointData.keys()

  volume = inputs[0]
  
  print 'start'

  mx,MX,my,MY,mz,MZ = volume.VTKObject.GetBounds()
  ll = np.array([mx,my,mz])
  ur = np.array([MX,MY,MZ])
  diagonal = ur - ll

  print 'aa', diagonal

  if len(inputs) > 1:
    points = inputs[1].GetPoints()
  else:
    i = int(starts)
    points = [[ll[i] + random.random()*diagonal[i] for i in range(3)] for j in range(starts)]

  if asscale == -1:
    asscale = np.linalg.norm(diagonal) / sscale

  volume = inputs[0]
  if volume.VTKObject.IsA('vtkImageData'):
    is_vtu = False
  elif volume.VTKObject.IsA('vtkUnstructuredGrid'):
    is_vtu = True
  else:
    print 'wha?'
    return

  print 'bb'

  if arrayName == "" and len(volume.PointData.keys()) > 1:
    print 'need to know what variable to operate on'
    return 

  if arrayName == -1:
    array = volume.GetPointData().GetArray(0)
  elif arrayName in volume.PointData.keys():
    array = volume.PointData[arrayName]
  else:
    print 'can\'t find requested data array:', arrayName
    print volume.PointData.keys()
    return 

  if target == 99999.0:
    mv = np.min(array)
    MV = np.max(array)
    target = (mv + MV) / 2.0

  if allabove:
    def map(a):
      if a == None: return 0.0
      if a > target: return 1.0
      d = target - a
      if d > kernelwidth: return 0.0
      else: return pow(1.0 - (d / kernelwidth), power)
  elif allbelow:
    def map(a):
      if a == None: return 0.0
      if a < target: return 1.0
      d = a - target
      if d > kernelwidth: return 0.0
      else: return pow(1.0 - (d / kernelwidth), power)
  else:
    def map(a):
      if a == None: 
        return 0.0
      d = target - a
      if d < 0: d = -d
      if d > kernelwidth: 
        return 0.0
      else:
        return pow(1.0 - (d / kernelwidth), power)

  samples = Samples(volume)
  interp  = Interpolator(volume)

  initial_points = []
  initial_pqs = []

  current_points = []
  current_pqs = []

  print 'dd'

  indx = 0
  for i, p in enumerate(points):
    if interp.Locate(p):
      v = interp.Interpolate(p, array)
      q = map(v)
      if q > 0.0:
        initial_points.append(p)
        initial_pqs.append(q)
        samples.add(interp, p, q, indx)
        indx = indx + 1

  current_points = list(initial_points)
  current_pqs = list(initial_pqs)

  print 'ee'

  misses = [0]*len(initial_points)

  done = False
  indx = 0

  print 'ff'

  while not done and samples.num() < nsamples:
    print samples.num(), indx
    if misses[indx] >= nmisses:
      print 'g'
      misses[indx] = 0
      current_points[indx] = initial_points[indx]
      current_pqs[indx] = initial_pqs[indx]
    print 'h'
    cpoint = current_points[indx] + np.random.normal(loc=0.0, scale=asscale, size=3)
    cv = interp.Interpolate(cpoint, array)
    cq = map(cv)
    print 'i'
    if cq > 0.0:
      if cq >= current_pqs[indx]:
        samples.add(interp, cpoint, cq, indx)
        print indx, samples.num()
        misses[indx] = 0
      else:
        u = np.random.rand()
        if u < cq/current_pqs[indx]:
          samples.add(interp, cpoint, cq, indx)
          print indx, samples.num()
          misses[indx] = 0
    else:
        misses[indx] = misses[indx] + 1
    print 'j'
    current_points[indx] = list(cpoint)
    current_pqs[indx] = cq
    indx = indx + 1
    if indx >= len(misses):
      indx = 0

  samples.stuff_vtu(output)
  return
