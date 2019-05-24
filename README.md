# ParaviewCustomModules

## Overview

This repo contains a set of custom modules for use with ParaView[1].  These include:

**SampleSurface**: Given a PolyData surface (think Contour output) sample randomly and independently of triangle size.  I use it to seed streamlines in a contour just below sea level.  Remember this data is not flat; I contour the magnitude of the coordinates to get distance to the center of the world, which I contour at 0.999 to get a surface plane.

**JoinStreamlines**:  Join the forward and backward streamlines from each seed point into a single polyline.  

**SampleArcLength**: Sample a set of poly lines based on arc length.   I think I can use integration time also to get samples equally spaced in time rather than space.

**OrientedGlyph**: Takes two inputs: a set of samples with at least two direction vectors and a polydata glyph.  Places the glyph at each sample point in the coordinate space defined by the two given directional vectors and their cross product. Using SampleArcLength on JoinStreamline’d VTK streamlines produces tangent and normal.

**CurvatureStreamlines**: Calculates local curvature in polylines.  Creates a curvature map using a downsampled version of the input polylines to avoid noise in almost-straight or oversampled sections.    Dot the result with the local up vector to get the direction that eddies are spinning.

**TrimStreamlines**:  Trim a streamline to an interval of integration time, given by two properties.

To use these, run 'Make' in the download directory.   This will invoke Kitware's 	python_filter_generator.py to create .xml files for each.   These can then be added to Paraview using Tools->Manage Plugins.

## License

ParaviewCustomModules have been developed primarily by Greg Abram with funding from the US Department of Energy's Office of Science ASCR program (Dr. Laura Biven, program manager) through sub-award from Los Alamos National Laboratory; and from the Texas Advanced Computing Center at the University of Texas at Austin.

ParaviewCustomModules is licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License. A copy of the License is included with this software in the file `LICENSE`. If your copy does not contain the License, you may obtain a copy of the License at: [Apache License Version 2.0][2]. 
Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.  

Copyright (c) 2019 The University of Texas at Austin. All rights reserved.

