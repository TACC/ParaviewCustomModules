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

Name = 'MH4'
Label = 'MH4'
Help = ''

NumberOfInputs = 2
InputDataType = ['vtkPolyData', 'vtkUnstructuredGrid', 'vtkImageData']
OutputDataType = 'vtkUnstructuredGrid'
ExtraXml = ''

Properties = dict(
  arrayName = 'timeMonthly_avg_ecosysTracers_DON',
  starts = 25,
  power = 1.0,
  sscale = 10000,
  loop_count = 1,
  target = 999999.0,
  spread = -1.0,
  maxpoints = -1,
)

def RequestData():
  import vtk
  import random
  import numpy as np
  from vtk.numpy_interface import dataset_adapter as dsa
  import paraview.vtk.util.numpy_support as vnp
  
  def P(v):
    if v < 0.0:
      return 0
    if power != 1.0:
      try:
        pv = pow(v, power)
        return pv
      except:
        print('E', v, power)
        return 0
    else:
      return v
  
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
        if cid < 0:
          self.xyz = []
          return False
        idl = vtk.vtkIdList()
        self.dset.GetCellPoints(cid, idl)
        self.ids = [idl.GetId(i) for i in range(idl.GetNumberOfIds())]
        #print("LOCATE cid", cid)
        #print("vids", self.ids)
        #print('wts', self.wts)
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
      self.V = []
      self.PV = []
      self.I = []
      for i in dset.PointData.keys():
        self.vars.append([i, dset.PointData[i], []])
  
    def num(self):
      return len(self.points)
  
    def add(self, I, p, v, pv, i):
      vals = [] 
      for var in self.vars:
        value = I.Interpolate(p, var[1])
        if value == None:
          value = -99999
        vals.append(value)
      self.points.append(p)
      self.V.append(v)
      self.PV.append(pv)
      self.I.append(i)
      for j,var in enumerate(self.vars):
        var[2].append(vals[j])
  
    def stuff_vtu(self, outpt):
      print("stuff", len(self.points), len(self.V), len(self.PV), len(self.I))
      outpt.SetPoints(dsa.VTKArray(np.array(self.points).astype('f4')))
      outpt.PointData.append(dsa.VTKArray(np.array(self.V).astype('f4')), 'V')
      outpt.PointData.append(dsa.VTKArray(np.array(self.PV).astype('f4')), 'PV')  
      outpt.PointData.append(dsa.VTKArray(np.array(self.I).astype('f4')), 'I')
      ct = dsa.numpyTovtkDataArray(np.array([vtk.VTK_VERTEX]*outpt.GetNumberOfPoints()).astype('u1'))
      co = dsa.numpy_support.numpy_to_vtkIdTypeArray(np.array(range(0, 2*outpt.GetNumberOfPoints(), 2)).astype('i8'))
      ca = vtk.vtkCellArray()
      for i in range(outpt.GetNumberOfPoints()):
        ca.InsertNextCell(1, [i])
      outpt.VTKObject.SetCells(ct, co, ca)
      for v in self.vars:
        outpt.PointData.append(dsa.VTKArray(np.array(v[2]).astype('f4')), v[0])
  
  np.random.seed(12346)
  
  volume = inputs[0]
  
  if volume.VTKObject.IsA('vtkImageData'):
    is_vtu = False
  elif volume.VTKObject.IsA('vtkUnstructuredGrid'):
    is_vtu = True
  else:
    print('wha?')
    return
  
  components = volume
  
  samples = Samples(volume)
  interp  = Interpolator(volume)
  
  if target == 999999.0:
    target = np.max(volume.PointData[arrayName])
  
  if spread < 0:
    a = target - np.min(volume.PointData[arrayName])
    b = np.max(components.PointData[arrayName]) - target
    if a > b: 
      spread = 0.1*a
    else:
      spread = 0.1*b

  if sscale == -1:
    from math import sqrt
    mx,my,mz,Mx,My,Mz = volume.VTKObject.GetBounds()
    dx = Mx - mx
    dy = My - my
    dz = Mz - mz
    sscale = 0.01 * sqrt(dx*dx + dy*dy + dz*dz)
  
  print('target', target, 'spread', spread, 'sscale', sscale)
  
  # Array is the PDF
  array = 1.0 - np.minimum(np.abs(volume.PointData[arrayName] - target) / spread, 1.0)
  
  initial_points = []
  initial_pqs = []
  
  r = np.random.rand(len(array))
  selections = array > r

  # if there are too many...
  if np.sum(selections) > starts:
    selections = selections & (((np.sum(selections)/float(starts)) * np.random.random(len(selections)) < 1))

  if is_vtu:
    pts = volume.Points[selections]
  else:
    indices = np.arange(len(array)).astype('i4')[selections]
    pts = []
    for i in indices:
      pts.append(volume.VTKObject.GetPoint(i))
      
  vs = array[selections]
  
  for p,v in zip(pts, vs):
    initial_points.append(p)
    initial_pqs.append(v) 
  
  
  print('target', target, 'spread', spread, 'seeds', len(initial_points))
  current_points = list(initial_points)
  current_pqs    = list(initial_pqs)
  misses         = [0]*len(initial_points)
  steps          = [0]*len(initial_points)
  
  done = False
  indx = 0
  
  accept_count = 0
  
  permute = np.arange(len(initial_points))
  np.random.shuffle(permute)
  
  for l in range(loop_count):
    if maxpoints > 0 and samples.num() >= maxpoints:
      break
  
    for indx in permute:
      if maxpoints > 0 and samples.num() >= maxpoints:
        break
  
      p0 = initial_points[indx]
      v0 = initial_pqs[indx]
  
      p1 = p0 + np.random.normal(loc=0.0, scale=sscale, size=3)
      v1 = interp.Interpolate(p1, array)
      if not v1: 
        continue
  
      accept = 0
      if v1 >= v0:
        accept = 1
      else:
        u = np.random.rand()
        if u < v1/v0:
          accept = 1
  
      if accept:
        samples.add(interp, p1, v1, v1, 0)
  
  samples.stuff_vtu(output)
  return
  
