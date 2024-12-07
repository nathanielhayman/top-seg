import pyvista as pv
import numpy as np

import ctypes as ct

from cutil import *

verts = np.array([
    # target curve
    [0, 0, 0],
    [0, 0, 1],
    [0, 0, 2],
    [0, 0, 3],
    [0, 0.5, 4],
    [0, 1, 5],
    [0, 1, 6],

    # anomaly
    [0.25, 0, 1],
    [0.5, 0, 2],
    [0.25, 0, 3],
])

ref_verts = np.array([
    # reference curve
    [1, 0, 0],
    [1, 0, 1],
    [1, 0, 2],
    [1, 0, 3],
    [1, 0.5, 4],
    [1, 1, 5],
    [1, 1, 6],
])

conn = np.array([
    # target line
    7, 0, 1, 2, 3, 4, 5, 6,

    # anomaly line
    5, 0, 7, 8, 9, 4,

    # null terminator
    0
], dtype=np.int32)

ref_conn: np.ndarray = np.array([
    7, 0, 1, 2, 3, 4, 5, 6, 0
], dtype=np.int32)

target   = pv.PolyData(verts, lines=conn)

reference = pv.PolyData(ref_verts, lines=ref_conn)

plotter = pv.Plotter()

_ = plotter.add_mesh(target, color='blue', render_points_as_spheres=True,
                     point_size=1, show_vertices=True)

_ = plotter.add_mesh(reference, color='red', render_points_as_spheres=True,
                     point_size=1, show_vertices=True)

plotter.show(cpos='xy')