import ctypes
import numpy as np

class ThisStruct(ctypes.Structure):
    _fields_ = [("nx", ctypes.c_int),
                ("ny", ctypes.c_int),
                ("x", ctypes.POINTER(ctypes.c_double)),
                ("y", ctypes.POINTER(ctypes.c_double))
    ]

nx = 10;
ny = 10;
    
nstruct = 2
TheseStructs = ThisStruct * nstruct

xs, ys = [None]*nstruct, [None]*nstruct
for i in range(nstruct):
    xs[i] = np.array([i*10 + ii for ii in range(nx)], dtype=np.float64)
    ys[i] = np.array([-i*10 + -ii for ii in range(ny)], dtype=np.float64)

thesestructs = TheseStructs()
for i in range(nstruct):
    thesestructs[i].nx = nx
    thesestructs[i].ny = ny
    thesestructs[i].x = xs[i].ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    thesestructs[i].y = ys[i].ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    

lib = ctypes.cdll.LoadLibrary("./clib.so")
lib.thisfunction.restype = ctypes.c_int
lib.thisfunction.argtypes = [ctypes.c_int, ctypes.POINTER(TheseStructs)]

lib.thisfunction(nstruct, ctypes.pointer(thesestructs))
                             
                             


