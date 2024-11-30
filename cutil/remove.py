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
        ("points", ct.POINTER(ct.c_int)),
        ("plen", ct.c_int),
        ("dlen", ct.c_int)
    ]

class Point(ct.Structure):
    _fields_ = [
        ("data", ct.POINTER(ct.c_int)),
        ("seg_i", ct.c_int)
    ]

class Segmented(ct.Structure):
    _fields_ = [
        ("points", ct.POINTER(Point)),
        ("rsize", ct.c_size_t),
        ("n_seg", ct.c_int),
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

def nearest_neighbors(ref: POINTER[Segmented], mesh: POINTER[Mesh], 
                      seg: POINTER[Segmented]):
    nn(ref, mesh, seg)

nns = lib.c_neighbors_from_segmented
nns.restype = None
nns.argtypes = [
    ct.POINTER(Segmented),
    ndpointer(ct.c_int, flags="C_CONTIGUOUS"),
    ndpointer(ct.c_int, flags="C_CONTIGUOUS"),
    ct.c_int,
    ct.POINTER(Segmented)
]

def neighbors_from_segmented(ref: POINTER[Segmented], pts: np.ndarray,
                             skeleton: np.ndarray
                             ) -> POINTER[Segmented] | None:
    tseg: POINTER[Segmented] | None = None

    nns(ref, pts, skeleton, skeleton.size, tseg)

    return tseg

spc = lib.c_seg_pv2c
spc.restype = None
spc.argtypes = [
    ndpointer(ct.c_double, flags="C_CONTIGUOUS"),
    ndpointer(ct.c_int, flags="C_CONTIGUOUS"),
    ct.c_int,
    ct.POINTER(Segmented)
]

something = ct.pointer(ct.c_int(10))

def seg_pv2c(pts: np.ndarray, 
             lines: np.ndarray) -> POINTER[Segmented] | None:
    seg: POINTER[Segmented] | None = None

    spc(pts, lines, lines.size, seg)

    return seg

scp = lib.seg_c2pv
scp.restype = None
scp.argtypes = [
    ct.POINTER(Segmented),
    ndpointer(ct.c_int, flags="C_CONTIGUOUS")
]

def seg_c2pv(seg: POINTER[Segmented]) -> np.ndarray:
    lines = np.ndarray(4)

    scp(seg, lines)

    return lines