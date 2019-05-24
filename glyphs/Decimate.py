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

import sys

if len(sys.argv) != 4:
  print 'syntax: pvpython', sys.argv[0], 'input.vtu decimation output.vtu'
  print 'where decimation is the reduction ratio'
  sys.exit(1)

from paraview.simple import *

ifile = sys.argv[1]
decimation = float(sys.argv[2])
ofile = sys.argv[3]

vtu = XMLUnstructuredGridReader(FileName=[ifile])
extractSurface1 = ExtractSurface(Input=vtu)
clean1 = Clean(Input=extractSurface1)
decimate1 = Decimate(Input=clean1)
decimate1.TargetReduction = decimation
clean2 = Clean(Input=decimate1)
generateSurfaceNormals1 = GenerateSurfaceNormals(Input=clean2)
generateSurfaceNormals1.FeatureAngle = 90
appendDatasets1 = AppendDatasets(Input=generateSurfaceNormals1)
SaveData(ofile, proxy=appendDatasets1)
