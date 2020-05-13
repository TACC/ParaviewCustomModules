from paraview.util.vtkAlgorithm import *

#------------------------------------------------------------------------------
# A filter example.
#------------------------------------------------------------------------------
@smproxy.filter()
@smproperty.input(name="InputDataset", port_index=0)
@smdomain.datatype(dataTypes=["vtkUnstructuredGrid"], composite_data_supported=False)

class SampleDensity(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1, outputType="vtkUnstructuredGrid")
        self.vars = ['aaa', 'bbb']

    def FillInputPortInformation(self, port, info):
        info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid")
        return 1

    @smproperty.stringvector(name="StringInfo", information_only="1")
    def GetStrings(self):
        return self.vars

    @smproperty.stringvector(name="String", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="StringInfo" function="StringInfo"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def SetString(self, value):
        self.inputArray = value
        self.Modified()
        return 1

    @smproperty.intvector(name="NSamples", label="Number of samples", default_values=1000)
    @smdomain.intrange(min=0, max=1000000)
    def SetNSamples(self, x):
        self.nSamples = x
        self.Modified()

    def RequestInformation(self, request, inInfoVec, outInfoVec):
        from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid
        executive = self.GetExecutive()
        input0 = vtkUnstructuredGrid.GetData(inInfoVec[0], 0)
        self.vars = []
        for i in range(input0.GetCellData().GetNumberOfArrays()):
          self.vars.append(input0.GetCellData().GetArray(i).GetName())
        print("here", self.vars)
        return 1

    def RequestData(self, request, inInfoVec, outInfoVec):
        from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid
        import sys, os
        import vtk
        from vtk.numpy_interface import dataset_adapter as dsa
        import ctypes as C
        import numpy as np
        
        if 'HOME' in os.environ:
          ccode_so = os.environ['HOME'] + '/python/_sample.so'
        elif 'HOMEPATH' in os.environ:
          ccode_so = os.environ['HOMEPATH'] + '/python/_sample.so'
        else:
          code_so =  '_sample.so'


        if not os.path.isfile(ccode_so):
          print('can\'t find so:', ccode_so)
          return

        print("found", ccode_so)
        
        ccode = C.CDLL(ccode_so)
        if not ccode:
          print('failed to load ', ccode.so)


        ccode.Sample.argtypes = [C.c_int, C.c_int, C.c_int,
                                 np.ctypeslib.ndpointer(C.c_float, flags="C_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(C.c_float, flags="C_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(C.c_int, flags="C_CONTIGUOUS")]
        
        ccode.GetNumberOfSamples.restype = C.c_int
        ccode.GetSamples.restype = C.c_void_p
        
        inpt = dsa.WrapDataObject(vtkUnstructuredGrid.GetData(inInfoVec[0], 0))

        pdf = self.inputArray
        nSamples = self.nSamples

        points = np.ascontiguousarray(inpt.Points).astype('f4')
        pdf    = np.ascontiguousarray(inpt.CellData[pdf]).astype('f4')
        tets   = np.ascontiguousarray(inpt.Cells).astype('i4')
        ccode.Sample(nSamples, inpt.GetNumberOfPoints(), inpt.GetNumberOfCells(), points, pdf, tets)
        n = ccode.GetNumberOfSamples()
        s = np.ctypeslib.as_array(C.cast(ccode.GetSamples(), C.POINTER(C.c_float)),shape=(n, 3))
        samples = dsa.WrapDataObject(vtk.vtkUnstructuredGrid())
        co = dsa.numpy_support.numpy_to_vtkIdTypeArray(np.arange(n).astype('i8')*2)
        ca = vtk.vtkCellArray()
        ca.SetCells(n, dsa.numpy_support.numpy_to_vtkIdTypeArray(np.column_stack(([1]*n, range(n))).reshape((2*n,))))
        ct = dsa.numpyTovtkDataArray(np.array([vtk.VTK_VERTEX]*n).astype('u1'))
        samples.VTKObject.SetCells(ct, co, ca)
        samples.Points = dsa.numpy_support.numpy_to_vtk(s, deep=1)
        
        outpt = vtkUnstructuredGrid.GetData(outInfoVec, 0)
        outpt.ShallowCopy(samples.VTKObject)

        ccode.Cleanup()
        return 1
        

