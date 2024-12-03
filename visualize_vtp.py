import pyvista as pv

import vtk

reader = vtk.vtkXMLPolyDataReader()
reader.SetFileName("artery.vtp")
reader.Update()
data = reader.GetOutput()

mesh = pv.wrap(data)

plotter = pv.Plotter()

plotter.add_mesh(mesh, color='gray')

plotter.show()