import numpy as np

import os

import ctypes as ct
from numpy.ctypeslib import ndpointer

import typing

# Kudos to KubaO:
# https://github.com/python/mypy/issues/7540
class POINTER[T](ct._Pointer):
    class _GenericAlias(typing.GenericAlias): # type: ignore
        def __repr__(self):
            val = super().__repr__()
            ibra = val.find('[')
            idot = val.rfind('.', 0, ibra)
            return f"{val[:idot+1]}POINTER{val[ibra:]}"

    def __class_getitem__(cls, *args):
        ptrtype = ct.POINTER(*args)
        alias = POINTER._GenericAlias(ptrtype, *args)
        return alias

class Path(ct.Structure):
    _fields_ = [
        ("points", ct.POINTER(ct.c_int32)),
        ("plen", ct.c_int32),
        ("dlen", ct.c_int32)
    ]

class Point(ct.Structure):
    _fields_ = [
        ("data", ct.POINTER(ct.c_int32)),
        ("seg_i", ct.c_int32)
    ]

class Segmented(ct.Structure):
    _fields_ = [
        ("points", ct.POINTER(Point)),
        ("rsize", ct.c_size_t),
        ("n_seg", ct.c_int32),
        ("paths", ct.POINTER(Path))
    ]

class Mesh(ct.Structure):
    _fields_ = [
        ("verts", ndpointer(ct.c_double, flags="C_CONTIGUOUS")),
        ("vsize", ct.c_size_t),
        ("faces", ndpointer(ct.c_int, flags="C_CONTIGUOUS")),
        ("fsize", ct.c_size_t)
    ]

lib_path = os.path.abspath("./libremove.dll")

lib = ct.cdll.LoadLibrary(lib_path)

nn = lib.c_nearest_neighbors
nn.restype = None
nn.argtypes = [
    ct.POINTER(Segmented),
    ct.POINTER(Mesh),
    ct.POINTER(Segmented)
]

def nearest_neighbors(ref: Segmented, mesh: Mesh, 
                      seg: Segmented):
    nn(ct.byref(ref), ct.byref(mesh), ct.byref(seg))

nns = lib.c_neighbors_from_segmented
nns.restype = None
nns.argtypes = [
    ct.POINTER(Segmented),
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ndpointer(ct.c_int32, flags="C_CONTIGUOUS"),
    ct.c_int32,
    ct.POINTER(Segmented)
]

def neighbors_from_segmented(ref: Segmented, pts: np.ndarray,
                             skeleton: np.ndarray
                             ) -> Segmented:
    tseg = Segmented(
        ct.POINTER(Point)(),
        0,
        0,
        ct.POINTER(Path)()
    )

    nns(ct.byref(ref), pts, skeleton, skeleton.size, ct.byref(tseg))

    return tseg

spc = lib.c_seg_pv2c
spc.restype = None
spc.argtypes = [
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ndpointer(ct.c_int32, flags="C_CONTIGUOUS"),
    ct.c_int32,
    ct.POINTER(Segmented)
]

def seg_pv2c(pts: np.ndarray, 
             lines: np.ndarray) -> Segmented:
    seg: Segmented = Segmented(
        ct.POINTER(Point)(),
        0,
        0,
        ct.POINTER(Path)()
    )

    spc(pts, lines, lines.size, ct.byref(seg))

    return seg

scp = lib.seg_c2pv
scp.restype = np.ndarray
scp.argtypes = [
    ct.POINTER(Segmented)
]

def seg_c2pv(seg: Segmented) -> np.ndarray:
    return scp(ct.byref(seg))