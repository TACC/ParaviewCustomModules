Name = 'JoinStreamlines'
Label = 'Join Streamlines'
Help = 'Join forward and backward fragments of streamlines'

NumberOfInputs = 1
InputDataType = 'vtkPolyData'
OutputDataType = 'vtkUnstructuredGrid'
ExtraXml = ''

Properties = dict(
)

def RequestData():
  import numpy as np
  from vtk.numpy_interface import dataset_adapter as dsa

  sl = self.GetPolyDataInput()
  nsl = dsa.WrapDataObject(sl)

  osl = self.GetUnstructuredGridOutput()

  osl.SetPoints(sl.GetPoints())
  for i in range(sl.GetPointData().GetNumberOfArrays()):
    osl.GetPointData().AddArray(sl.GetPointData().GetArray(i))

  arclen = np.zeros(len(nsl.Points))

  sids = nsl.CellData['SeedIds']
  unique_sids = np.unique(sids)

  polylines = []
  idlists = [vtk.vtkIdList() for i in range(2)]

  if 'IntegrationTime' not in nsl.PointData.keys() or 'Normals' not in nsl.PointData.keys():
    print "need integration time and Normals variables"
    return

  for indx, seed in enumerate(unique_sids):
    print 'aa'
    parts = np.where(sids == seed)[0]
    sl.GetCellPoints(parts[0], idlists[0])
    w0 = [idlists[0].GetId(i) for i in range(idlists[0].GetNumberOfIds())]
    if len(parts) == 2:
      sl.GetCellPoints(parts[1], idlists[1])
      w1 = [idlists[1].GetId(i) for i in range(1, idlists[1].GetNumberOfIds())]
      nsl.PointData['Normals'][w1] = 0 - nsl.PointData['Normals'][w1]
      t0 = nsl.PointData['IntegrationTime'][w0[-1]]
      t1 = nsl.PointData['IntegrationTime'][w1[-1]]
      if t0 > t1:
        tmin = t1
        ids = w1[::-1] + w0
      else:
        tmin = t0
        nsl.PointData['Normals'][w0] = 0 - nsl.PointData['Normals'][w0]
        ids = w0[::-1] + w1
      nsl.PointData['IntegrationTime'][ids]  = nsl.PointData['IntegrationTime'][ids] - tmin
    else:
      ids = w0
    polylines.append(ids)
    p = nsl.Points[ids]
    pv = p[1:] - p[:-1]
    l = np.linalg.norm(pv, axis=1)
    arclen[ids] = np.hstack(([0], np.cumsum(l)))

  ct = dsa.numpyTovtkDataArray(np.array([vtk.VTK_POLY_LINE]*len(polylines)).astype('u1'))
  line_lengths = [len(p)+1 for p in polylines]  # number of vertices + 1 for count
  co = dsa.numpy_support.numpy_to_vtkIdTypeArray(np.hstack(([0], np.cumsum(line_lengths)[:-1])))

  ca = vtk.vtkCellArray()
  for pl in polylines:
    ca.InsertNextCell(len(pl), pl)

  osl.SetCells(ct, co, ca)

  cid = dsa.numpyTovtkDataArray(np.arange(len(polylines)).astype('i4'))
  cid.SetName('cell id')
  osl.GetCellData().AddArray(cid)

  al = dsa.numpyTovtkDataArray(arclen.astype('f4'))
  al.SetName('arclen')
  osl.GetPointData().AddArray(al)
