import pyvista as pv

from vtkmodules.vtkCommonDataModel import vtkIterativeClosestPointTransform

import vtk

from vmtk import pypes

FILENAME = "./data/topcow_mr_001_0000.nii.gz"

REF_FILENAME = "./data/topcow_mr_001.nii.gz"

# # Read reference mesh

reader = pv.get_reader(REF_FILENAME)
mesh = reader.read()

ref_mesh = mesh.threshold(1)

if ref_mesh is None:
    exit("Could not threshold ref mesh")

# # Read target mesh

reader = pv.get_reader(FILENAME)
tgt_mesh = reader.read()

tgt_mesh = tgt_mesh.threshold(-2.1e4)

if tgt_mesh is None:
    exit("Could not threshold target mesh")

tgt_mesh = tgt_mesh.clip_box(ref_mesh.bounds, invert=False)

if tgt_mesh is None:
    exit("Could not clip target mesh")

tgt_mesh = tgt_mesh.connectivity('largest')

if tgt_mesh is None:
    exit("Could not get connectivity of target mesh")

geo_filter = vtk.vtkGeometryFilter()
geo_filter.SetInputData(tgt_mesh)
geo_filter.Update()

poly_data = geo_filter.GetOutput()

smooth_filter = vtk.vtkSmoothPolyDataFilter()
smooth_filter.SetInputData(poly_data)
smooth_filter.SetNumberOfIterations(30)
smooth_filter.SetRelaxationFactor(0.1)
smooth_filter.FeatureEdgeSmoothingOff()
smooth_filter.BoundarySmoothingOn()  
smooth_filter.Update()

smooth_data = smooth_filter.GetOutput()

tgt_mesh = pv.wrap(smooth_data)

if tgt_mesh is None:
    exit("could not wrap target mesh")

tgt_mesh = tgt_mesh.clean()

writer = vtk.vtkXMLPolyDataWriter()
writer.SetFileName("artery.vtp")
writer.SetInputData(tgt_mesh)
writer.Write()

print("wrote")

my_pype = pypes.PypeRun("vmtknetworkextraction -ifile artery.vtp -advancementratio 1 -ofile new_artery.vtp")

# surface = my_pype.GetScriptObject('vmtknetworkextraction', '0').Surface

# print("ran vmtk network extraction")

# reader = vtk.vtkXMLPolyDataReader()
# reader.SetFileName("new_artery.vtp")
# reader.Update()
# data = reader.GetOutput()

mesh = pv.read("new_artery.vtp")

print(f"points: {mesh.n_points}")

print(f"polys: {mesh.n_cells}")

plotter = pv.Plotter()

plotter.add_mesh(tgt_mesh, color='gray')

plotter.show()

# print(data)