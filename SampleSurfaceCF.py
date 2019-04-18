Name = 'RandomSampleSurface'
Label = 'Sample Surface'
Help = 'Sample a surface randomly'

NumberOfInputs = 1
InputDataType = 'vtkPolyData'
OutputDataType = 'vtkUnstructuredGrid'
ExtraXml = ''

Properties = dict(
nsamples = 100
)

def RequestData():
  import random
  from vtk.numpy_interface import dataset_adapter as dsa
  from numpy import cross, ascontiguousarray, array, column_stack, arange
  from numpy.linalg import norm

  inp = inputs[0]

  tf = vtk.vtkTriangleFilter()
  tf.SetInputData(inp.VTKObject)
  tf.Update()

  dobj = dsa.WrapDataObject(tf.GetOutput())
  
  nCells = dobj.GetNumberOfCells()

  tris = dobj.GetPolygons().reshape((-1, 4))[:,1:]

  points = dobj.GetPoints()

  p0 = points[tris[:,0]]
  p1 = points[tris[:,1]]
  p2 = points[tris[:,2]]

  areas = norm(cross(p1-p0, p2-p0), axis=1)/2.0
  area_per_sample = sum(areas) / nsamples
  samples_per_triangle = [(a / area_per_sample) for a in areas]
  samples_per_triangle = [(int(i) + 1) if (random.random() < (i - int(i))) else int(i) for i in samples_per_triangle]
  sampled_triangles = [[i]*j for i,j in enumerate(samples_per_triangle)]
  sampled_triangles = [j for i in sampled_triangles for j in i]
  selected_triangles = tris[sampled_triangles,:]

  p = points[selected_triangles[:,0]]
  q = points[selected_triangles[:,1]]
  r = points[selected_triangles[:,2]]

  qp = q - p
  rp = r - p

  unit = [[random.random(), random.random()] for i in range(len(selected_triangles))]
  unit = [u if (u[0]+u[1]) < 1 else [1-u[0], 1-u[1]] for u in unit]
  
  v = [u[0]*q + u[1]*r for u,q,r in zip(unit, qp, rp)]
  samples = p + v

  output.SetPoints(dsa.VTKArray(samples))

  # w0 = array([random.random()/3.0 for t in selected_triangles])
  # w1 = array([random.random()/3.0 for t in selected_triangles])
  # w2 = 1.0 - (w0 + w1)
  # starting_points = points[selected_triangles[:,0]]*w0 + points[selected_triangles[:,1]]*w1 + points[selected_triangles[:,2]]*w2
  # npts = len(starting_points)
  # output.SetPoints(dsa.VTKArray(starting_points))
