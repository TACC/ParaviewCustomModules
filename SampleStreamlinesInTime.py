Name = 'SampleStreamlinesInTime'
Label = 'Sample Streamlines In Time'
Help = 'One sample at of each pathline at t (0 to 1) of the way along'

NumberOfInputs = 1
InputDataType = 'vtkUnstructuredGrid'
OutputDataType = 'vtkUnstructuredGrid'
ExtraXml = ''

Properties = dict(
t = 0.5
)

def RequestData():
  import numpy as np
  from vtk.numpy_interface import dataset_adapter as dsa
  from math import ceil, floor

  sl = self.GetUnstructuredGridInput()
  nsl = dsa.WrapDataObject(sl)

  print 'hello'
  arclen = nsl.PointData['arclen']

  nv = nsl.Cells[nsl.CellLocations]    # number of verts in each line
  ns = nsl.CellLocations + 1           # index of first vertex in each line
  ne = ns + nv                         # index one past the last vertex

  lines = [nsl.Cells[i:j] for i,j in zip(ns, ne)] # divide into distinct lines

  iarrays = {'points': nsl.Points}    # initialize source arrays with input points
  oarrays = {'points': []}            # initialize destination arrays with (empty) points

  for n in nsl.PointData.keys():      
    iarrays[n] = nsl.PointData[n]     # add input point data arrays to source arrays
    oarrays[n] = []                   # add empty destination arrays

  for i,line in enumerate(lines):     # for each input line...

    llen = arclen[line[-1]]
    sample_x = t * llen               #   sample distance along line

    # print 'arclen', arclen[line]
    # print 'sample_x', sample_x

    # index of first elt greater than sample_x (or 0, in which case we use the last)

    interval_end = np.argmax(arclen[line] > sample_x) 
    if interval_end == 0: interval_end = len(line) - 1

    # print 'interval_end', interval_end

    # get indices of points and point-dependent data at either end of the interval
    starti = line[interval_end]
    endi   = line[interval_end - 1]

    # print 'starti', starti, 'endi', endi

    # interpolant value in interval
    d = (sample_x - arclen[starti] / (arclen[endi] - arclen[starti]))  # interpolant in interval

    for n in iarrays:                 #   for each array we are interpolating...
      ia = iarrays[n]                 #     input array
      sv = ia[starti]                 #     start values
      ev = ia[endi]                   #     end values
      v  = sv + d*(ev - sv)           #     interpolation
      if n == 'points': 
        print v
      oarrays[n].append(v)

  print 'AAA'

  ptsa = np.concatenate(oarrays['points']).reshape((-1, 3)).astype('f4')

  print 'BBB', oarrays['points'], ptsa, ptsa.shape

  oug = self.GetUnstructuredGridOutput()

  print 'CCCC'

  op = vtk.vtkPoints()
  op.SetNumberOfPoints(ptsa.shape[0])

  for i, p in enumerate(ptsa):
    op.InsertPoint(i, p[0], p[1], p[2])

  print 'DDDD'

  oug.SetPoints(op)

  for n in oarrays:
    if n != 'points':
      a = dsa.numpyTovtkDataArray(np.concatenate(oarrays[n]))
      a.SetName(n)
      oug.GetPointData().AddArray(a)

  print 'EEE'

  ct = dsa.numpyTovtkDataArray(np.array([vtk.VTK_VERTEX]*oug.GetNumberOfPoints()).astype('u1'))
  co = dsa.numpy_support.numpy_to_vtkIdTypeArray(np.array(range(0, 2*oug.GetNumberOfPoints(), 2)))

  print 'FFFG'

  ca = vtk.vtkCellArray()
  for i in range(oug.GetNumberOfPoints()):
    ca.InsertNextCell(1, [i])

  print 'GGG'
  oug.SetCells(ct, co, ca)
