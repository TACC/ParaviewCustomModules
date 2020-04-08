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

Name = 'MH'
Label = 'MH'
Help = ''

NumberOfInputs = 2
InputDataType = ['vtkPolyData', 'vtkUnstructuredGrid', 'vtkImageData']
OutputDataType = 'vtkUnstructuredGrid'
ExtraXml = ''

Properties = dict(
  arrayName = 'timeMonthly_avg_ecosysTracers_DON',
  nmisses = 10,
  nsamples = 10000,
  nstarts  = 10,
  power = 1.0,
  sscale = 10000,
  max_loop_count = 100000,
  minvalue = -1.0,
  every = 10,
)

def RequestData():
  import vtk
  import random
  import numpy as np
  from vtk.numpy_interface import dataset_adapter as dsa
  import paraview.vtk.util.numpy_support as vnp

  def P(v):
    if power != 1.0:
      return pow(v, power)
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
      self.V = []
      self.PV = []
      self.I = []
      for i in dset.PointData.keys():
        self.vars.append([i, dset.PointData[i], []])
  
    def num(self):
      return len(self.points)
  
    def add(self, I, p, v, pv, i):
      err = 0
      vals = []
      for var in self.vars:
        value = I.Interpolate(p, var[1])
        if value == None:
          err = 1
          print('oops', var[0])
          break
        vals.append(value)
      if err == 0:
        self.points.append(p)
        self.V.append(v)
        self.PV.append(pv)
        self.I.append(i)
        for j,var in enumerate(self.vars):
          # print("XXX", var[0], j)
          var[2].append(vals[j])

    def stuff_vtu(self, outpt):
      outpt.SetPoints(dsa.VTKArray(np.array(self.points).astype('f4')))
      outpt.PointData.append(dsa.VTKArray(np.array(self.V).astype('f4')), 'V')
      outpt.PointData.append(dsa.VTKArray(np.array(self.PV).astype('f4')), 'PV')
      outpt.PointData.append(dsa.VTKArray(np.array(self.I).astype('f4')), 'I')
      ct = dsa.numpyTovtkDataArray(np.array([vtk.VTK_VERTEX]*outpt.GetNumberOfPoints()).astype('u1'))
      co = dsa.numpy_support.numpy_to_vtkIdTypeArray(np.array(range(0, 2*outpt.GetNumberOfPoints(), 2)))
      ca = vtk.vtkCellArray()
      for i in range(outpt.GetNumberOfPoints()):
        ca.InsertNextCell(1, [i])
      outpt.VTKObject.SetCells(ct, co, ca)
      for v in self.vars:
        outpt.PointData.append(dsa.VTKArray(np.array(v[2]).astype('f4')), v[0])

  np.random.seed(12346)
  
  volume = inputs[0]

  array = volume.PointData[arrayName]
  if volume.VTKObject.IsA('vtkImageData'):
    is_vtu = False
  elif volume.VTKObject.IsA('vtkUnstructuredGrid'):
    is_vtu = True
  else:
    print('wha?')
    return

  samples = Samples(volume)
  interp  = Interpolator(volume)

  # This stuff thresholds out the part of the incoming dataset that at or above the selected
  # minvalue, then selects one initial point in each connected component of the result

  tf = vtk.vtkThreshold()  
  tf.SetInputData(volume.VTKObject)
  tf.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, arrayName)
  tf.ThresholdByUpper(minvalue)

  cf = vtk.vtkConnectivityFilter()      
  cf.SetInputConnection(tf.GetOutputPort())
  cf.SetExtractionModeToAllRegions();
  cf.ColorRegionsOn();

  cf.Update();
  components = dsa.WrapDataObject(cf.GetOutput())

  rids = components.PointData['RegionId']
  n,i = np.unique(rids, return_index=True)

  initial_points = components.Points[i]
  initial_pqs = []
  for p in initial_points:
    v = interp.Interpolate(p, array)
    pv = P(v)
    samples.add(interp, p, v, pv, 0)
    initial_pqs.append(P(interp.Interpolate(p, array)))

  current_points = list(initial_points)
  current_pqs    = list(initial_pqs)
  misses         = [0]*len(initial_points)
  steps          = [0]*len(initial_points)

  done = False
  indx = 0

  accept_count = 0
  loop_count = 0
  while not done and samples.num() < nsamples:

    loop_count = loop_count + 1
    if  loop_count % 1000 == 0:
      print(loop_count)

    if loop_count > max_loop_count:
      print("broke on total loop count")
      done = True

    if misses[indx] >= nmisses:
      misses[indx] = 0
      current_points[indx] = initial_points[indx]
      current_pqs[indx] = initial_pqs[indx]

    cpoint = current_points[indx] + np.random.normal(loc=0.0, scale=sscale, size=3)
    cv = interp.Interpolate(cpoint, array)
    if not cv:
      continue

    cq = P(cv)

    accept = 0
    if cq >= current_pqs[indx]:
      accept = 1
      misses[indx] = 0
    else:
      u = np.random.rand()
      if u < cq/current_pqs[indx]:
        accept = 1
        misses[indx] = 0
      else:
        accept = 0
        misses[indx] = misses[indx] + 1

    if accept:
      if accept_count % every == 0:
        samples.add(interp, cpoint, cv, cq, steps[indx])
      misses[indx] = 0
      steps[indx] = steps[indx] + 1
      current_points[indx] = list(cpoint)
      current_pqs[indx] = cq
      accept_count = accept_count + 1

    indx = indx + 1
    if indx >= len(misses):
      indx = 0

  samples.stuff_vtu(output)
  return
