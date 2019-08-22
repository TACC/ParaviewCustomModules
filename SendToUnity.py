Name = 'SendToUnity'
Label = 'SendToUnity'
Help = 'Send unstructured grid to Unity'

NumberOfInputs = 1
InputDataType = 'vtkUnstructuredGrid'
OutputDataType = 'vtkUnstructuredGrid'
ExtraXml = ''

Properties = dict(
  label = 'label',
  host = 'localhost',
  port = 1900
)

def RequestData():
  import socket
  from struct import pack, unpack
  from vtk import vtkUnstructuredGridWriter
  try:
    s = socket.socket()
    s.connect((host, port))
  except:
    print 'unable to connect to server', host, port
  else:
    print 'connected'

    for i in range(self.GetUnstructuredGridInput().GetPointData().GetNumberOfArrays()):
      self.GetUnstructuredGridOutput().GetPointData().AddArray(self.GetUnstructuredGridInput().GetPointData().GetArray(i))

    n = self.GetUnstructuredGridInput().GetPointData().GetNormals()
    if n:
      self.GetUnstructuredGridOutput().GetPointData().SetNormals(n)

    wrtr = vtkUnstructuredGridWriter()
    wrtr.SetInputData(self.GetUnstructuredGridInput())

    wrtr = vtkUnstructuredGridWriter()
    wrtr.WriteToOutputStringOn()
    wrtr.SetInputData(self.GetUnstructuredGridInput())
    wrtr.Update()

    lab = label.encode()
    s.send(pack('!i', len(lab)))
    s.send(lab)

    data = wrtr.GetOutputString().encode()
    s.send(pack('!i', len(data)))
    s.send(data)

    ack = s.recv(2).decode()
    print 'got', ack
