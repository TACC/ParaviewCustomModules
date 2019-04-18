Name = 'SampleStreamlines'
Label = 'Sample Streamlines'
Help = 'Sample streamlines evenly in arc len'

NumberOfInputs = 1
InputDataType = 'vtkUnstructuredGrid'
OutputDataType = 'vtkUnstructuredGrid'
ExtraXml = ''

Properties = dict(
nsamples = 100
)

def RequestData():
  import numpy as np
  from vtk.numpy_interface import dataset_adapter as dsa
  from math import ceil, floor

  sl = self.GetUnstructuredGridInput()
  nsl = dsa.WrapDataObject(sl)

  arclen = nsl.PointData['arclen']

  nv = nsl.Cells[nsl.CellLocations]    # number of verts in each line
  ns = nsl.CellLocations + 1           # index of first vertex in each line
  ne = ns + nv                         # index one past the last vertex

  lines = [nsl.Cells[i:j] for i,j in zip(ns, ne)] # divide into distinct lines

  llen = arclen[[l[-1] for l in lines]]   # length of each line
  totlen = sum(llen)                      # total length
  sdist = totlen / nsamples               # appx distance between samples

  nPerLine = np.array([int(ceil(l / sdist)) for l in llen])   # samples in each line
  nPerLine = np.where(nPerLine < 3, 3, nPerLine)

  iarrays = {'points': nsl.Points}    # initialize source arrays with input points
  oarrays = {'points': []}            # initialize destination arrays with (empty) points

  for n in nsl.PointData.keys():      
    iarrays[n] = nsl.PointData[n]     # add input point data arrays to source arrays
    oarrays[n] = []                   # add empty destination arrays

  for i,line in enumerate(lines):     # for each input line...
    ns = nPerLine[i]                  #   number in this line
    sdist1 = float(llen[i] / (ns-1))  #   inter-sample distance in this line
    x = arclen[line]                  #   X axis is arc len along line
    y = range(len(line))              #   Y is index along line
    s = [i*sdist1 for i in range(ns)] #   s's are the sample points along the line in arclen
    d = np.interp(s, x, y)            #   d's are the interpolant values
    ds = np.floor(d).astype('i4')     #   index of interval start
    dd = d - ds                       #   delta in interval
    de = np.where(dd == 0, ds, ds+1)  #   index of interval end (unless dd is zero)
    si = line[ds]                     #   offset of starting value in arrays
    se = line[de]                     #   offset of ending value in arrays
    for n in iarrays:                 #   for each array we are interpolating...
      ia = iarrays[n]                 #     input array
      sv = ia[si]                     #     start values
      ev = ia[si]                     #     end values
      v  = sv + dd*(ev - sv)          #     interpolation
      oarrays[n].append(v)

  ptsa = np.concatenate(oarrays['points']).astype('f4')

  oug = self.GetUnstructuredGridOutput()

  op = vtk.vtkPoints()
  op.SetNumberOfPoints(ptsa.shape[0])

  for i, p in enumerate(ptsa):
    op.InsertPoint(i, p[0], p[1], p[2])

  oug.SetPoints(op)

  for n in oarrays:
    if n != 'points':
      a = dsa.numpyTovtkDataArray(np.concatenate(oarrays[n]))
      a.SetName(n)
      oug.GetPointData().AddArray(a)

  ct = dsa.numpyTovtkDataArray(np.array([vtk.VTK_VERTEX]*oug.GetNumberOfPoints()).astype('u1'))
  co = dsa.numpy_support.numpy_to_vtkIdTypeArray(np.array(range(0, 2*oug.GetNumberOfPoints(), 2)))

  ca = vtk.vtkCellArray()
  for i in range(oug.GetNumberOfPoints()):
    ca.InsertNextCell(1, [i])

  oug.SetCells(ct, co, ca)
