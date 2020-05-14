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

Name = 'PDF'
Label = 'PDF'
Help = ''

NumberOfInputs = 1
InputDataType = ['vtkPolyData', 'vtkUnstructuredGrid', 'vtkImageData']
OutputDataType = 'vtkUnstructuredGrid'
ExtraXml = ''

Properties = dict(
  arrayName = 'RTData',
  one_side = True,
  power = 150.0,
  target = 0.0
)

def RequestData():
  import numpy as np
  import vtk
  from vtk.numpy_interface import dataset_adapter as dsa
  import paraview.vtk.util.numpy_support as vnp
  
  print("here 2")
  print(inputs[0].Points)
  print("here 2a")
  
  raw = inputs[0].PointData['RTData']
  if one_side:
    a = np.abs(target - np.where(raw > target, target, raw))
  else:  
    a = np.abs(target - raw)
  
  m = np.max(a)
  b = 1.0 - (a / np.max(a))
  
  c = np.power(b, power)

  print('here 3', len(inputs[0].Points))
  output.VTKObject.ShallowCopy(inputs[0].VTKObject)
  print('here 4', len(output.Points))
  output.PointData.append(c, 'PDF')
  
  print("done 2", output.VTKObject)
  return
