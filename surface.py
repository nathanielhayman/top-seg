import pyvista as pv

from skimage import measure

import numpy as np

FILENAME = "./data/topcow_mr_002_0000.nii.gz"

REF_FILENAME = "./data/topcow_mr_002.nii.gz"

PLOT_BOX = False

CONTOUR = False

reader = pv.get_reader(REF_FILENAME)

mesh = reader.read()

reference = mesh.threshold(1)

reader = pv.get_reader(FILENAME)

mesh = reader.read()

# REMOVE
# source = mesh.threshold(1)

# contoured = None

# if source is not None:
#     contoured = source.contour()

# conn = mesh.connectivity('specified', (0,1,2))

ref_bounds = None

# sourcez = None

if reference is not None:
    ref_bounds = reference.bounds

print(ref_bounds)

# -2.1e4
source = mesh.threshold(-2.3e4, progress_bar=True)

if source is None:
    exit("Source could not be thresholded")

source = source.clip_box(ref_bounds, invert=False, progress_bar=True)

if source is None:
    exit("Source could not be bound")

source = source.connectivity('largest', progress_bar=True)

if source is None:
    exit("Largest connectivity extraction failed")

resolution = (1000, 1000, 1000)

print(resolution)

print((resolution[0] + 1, resolution[1] + 1,
                                   resolution[2] + 1))

print((resolution[0] + 1, resolution[1] + 1,
                                   resolution[2] + 1))

bounds = source.bounds

# sk_grid = pv.ImageData()

# sk_grid.origin = (bounds[0], bounds[1], bounds[2])

# sk_grid.spacing = (
#     abs(bounds[0]-bounds[1]) / resolution[0],
#     abs(bounds[2]-bounds[3]) / resolution[1],
#     abs(bounds[4]-bounds[5]) / resolution[2]
# )

# sk_grid.dimensions = resolution[0] + 1, resolution[1] + 1, resolution[2] + 1

# sampled = sk_grid.sample(source)

print(source.points)

# volume = np.array(source.points)

# print(type(volume), volume.ndim, volume.shape)

# volume = source.cast_to_poly_points().delaunay_3d()

# if volume is not None:
#     verts, faces, _, _ = measure.marching_cubes(volume, level=0.5)

#     source = pv.PolyData(verts, faces)

resampled = None

if CONTOUR:
    s_data = pv.create_grid(source, dimensions=(
        abs(source.bounds[0]-source.bounds[1]),
        abs(source.bounds[2]-source.bounds[3]),
        abs(source.bounds[4]-source.bounds[5])
    ))

    resampled = s_data.sample(source)

    if resampled is not None:
        print(resampled.point_data.keys())

        source = resampled.contour(1, scalars=resampled['RegionId'],
                                   method='marching_cubes', progress_bar=True)
        
        print(source)

print(source)

#     sourcex = source.threshold([ref_bounds[0], ref_bounds[1]], scalars='x')
#     sourcey = sourcex.threshold([ref_bounds[2], ref_bounds[3]], scalars='y') # type: ignore
#     sourcez = sourcey.threshold([ref_bounds[4], ref_bounds[5]], scalars='z') # type: ignore

# source_mesh = sourcez

# if sourcez is not None:
#     source_mesh = sourcez.threshold(-1e4) # type: ignore

box = pv.Box(ref_bounds)

plotter = pv.Plotter()

_ = plotter.add_mesh(source, color='lightgray')

_ = plotter.add_mesh(reference, opacity=0.4)

if PLOT_BOX:
    _ = plotter.add_mesh(box, opacity=0.6)

# plotter.add_mesh_threshold(mesh, pointa=(-2e4,1), pointb=(0,1))

plotter.camera_position = 'xy'

plotter.show()