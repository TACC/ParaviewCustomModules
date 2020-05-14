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

Name = 'DecimateStreamlines'
Label = 'DecimateStreamlines'
Help = ''

NumberOfInputs = 1
InputDataType = ['vtkPolyData', 'vtkUnstructuredGrid', 'vtkImageData']
OutputDataType = 'vtkUnstructuredGrid'
ExtraXml = ''

Properties = dict(
  factor = 5,
)

def RequestData():
  import vtk
  import numpy as np
  from vtk.numpy_interface import dataset_adapter as dsa

  inpt = inputs[0]

  cl = inpt.GetCellLocations()
  ct = inpt.GetCellTypes()
  cells = inpt.GetCells()

  ncl = []
  ncells = []
  remap = []

  while len(cells) > 0:
    onids = cells[0]
    oids = cells[1:onids+1]
    nids = [oids[0]] + list(oids[factor:-1:factor]) + [oids[-1]]
    nnids = len(nids)
    ncell = [nnids] + list(range(len(remap), len(remap)+len(nids)))
    remap = remap + nids
    ncells = ncells + ncell
    cells = cells[onids+1:]

  new_points = inpt.GetPoints()[remap]
  ncells = dsa.VTKArray(ncells).astype('i8')
  ncl = dsa.VTKArray(np.array(ncl)).astype('i8')

  o = dsa.WrapDataObject(vtk.vtkUnstructuredGrid())
  o.SetPoints(dsa.VTKArray(np.array(new_points).astype('f4')))
  o.SetCells(ct, ncl, ncells)

  ipd = inpt.GetPointData()
  opd = o.GetPointData()
  for i in ipd.keys():
    opd.append(ipd[i][remap], i)

  icd = inpt.GetCellData()
  ocd = o.GetCellData()
  for i in icd.keys():
    ocd.append(icd[i], i)

  if 0 == 1:
    wrtr = vtk.vtkXMLUnstructuredGridWriter()
    wrtr.SetFileName('reduced.vtu')
    wrtr.SetInputData(o.VTKObject)
    wrtr.Write()
  else:
    output.VTKObject.ShallowCopy(o.VTKObject)

  
