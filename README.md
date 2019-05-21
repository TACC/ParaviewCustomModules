# ParaviewCustomModules

This repo contains a set of custom modules for use with Paraview.  These include:

SampleSurface: Given a PolyData surface (think Contour output) sample randomly and independently of triangle size.  I use it to seed streamlines in a contour just below sea level.  Remember this data is not flat; I contour the magnitude of the coordinates to get distance to the center of the world, which I contour at 0.999 to get a surface plane.

JoinStreamlines:  Join the forward and backward streamlines from each seed point into a single polyline.  

SampleArcLength: Sample a set of poly lines based on arc length.   I think I can use integration time also to get samples equally spaced in time rather than space.

OrientedGlyph: Takes two inputs: a set of samples with at least two direction vectors and a polydata glyph.  Places the glyph at each sample point in the coordinate space defined by the two given directional vectors and their cross product. Using SampleArcLength on JoinStreamline’d VTK streamlines produces tangent and normal.

CurvatureStreamlines: Calculates local curvature in polylines.  Creates a curvature map using a downsampled version of the input polylines to avoid noise in almost-straight or oversampled sections.    Dot the result with the local up vector to get the direction that eddies are spinning.

TrimStreamlines:  Trim a streamline to an interval of integration time, given by two properties.

To use these, run 'Make' in the download directory.   This will invoke Kitware's 	python_filter_generator.py to create .xml files for each.   These can then be added to Paraview using Tools->Manage Plugins.

