import pyvista as pv

import numpy as np

from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget, \
QSlider, QPushButton, QLabel

from PyQt5.QtCore import Qt

import pyvistaqt as pvqt

from vtkmodules.vtkCommonDataModel import vtkIterativeClosestPointTransform

import vtk

FILENAME = "./data/topcow_mr_002_0000.nii.gz"

REF_FILENAME = "./data/topcow_mr_001.nii.gz"

PLOT_BOX = False

CONTOUR = False

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        # Initialize Qt window

        self.setWindowTitle("Visualization")
        self.setGeometry(100, 100, 100, 100)

        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        layout = QVBoxLayout()
        central_widget.setLayout(layout)


        # Read reference mesh

        reader = pv.get_reader(REF_FILENAME)
        mesh = reader.read()

        self.ref_mesh = mesh.threshold(1)

        # Read target mesh

        reader = pv.get_reader(FILENAME)
        self.tgt_mesh = reader.read()

        self.alpha_mesh = reader.read()

        # REMOVE
        # source = mesh.threshold(1)

        # contoured = None

        # if source is not None:
        #     contoured = source.contour()

        # conn = mesh.connectivity('specified', (0,1,2))

        self.ref_bounds = None

        # sourcez = None

        if self.ref_mesh is not None:
            self.ref_bounds = self.ref_mesh.bounds

        # -2.1e4
        self.threshold(-2.3e4)

        # resolution = (1000, 1000, 1000)

        # print(resolution)

        # print((resolution[0] + 1, resolution[1] + 1,
        #                                 resolution[2] + 1))

        # print((resolution[0] + 1, resolution[1] + 1,
        #                                 resolution[2] + 1))

        # if reference is None:
        #     exit("Could not construct reference mesh")

        #     vtk.vtkUnstructuredGridTo

        # geo_filter = vtk.vtkGeometryFilter()
        # geo_filter.SetInputData(source)
        # geo_filter.Update()

        # poly_data = geo_filter.GetOutput()

        # writer = vtk.vtkXMLPolyDataWriter()
        # writer.SetFileName("artery.vtp")
        # writer.SetInputData(poly_data)
        # writer.Write()

        # myPype= pypes.PypeRun("vmtknetworkextraction -ifile ArteryObjAN129–10.vtp -advancementratio 1 -ofile centerlines/ArteryObjAN129–10.vtp")

        # icp = vtkIterativeClosestPointTransform()
        # icp.SetSource(source)
        # icp.SetTarget(reference)
        # icp.GetLandmarkTransform().SetModeToRigidBody()
        # icp.SetMaximumNumberOfLandmarks(100)
        # icp.SetMaximumMeanDistance(.00001)
        # icp.SetMaximumNumberOfIterations(500)
        # icp.CheckMeanDistanceOn()
        # icp.StartByMatchingCentroidsOn()
        # icp.Update()

        # aligned = source.transform(icp.GetMatrix())

        # bounds = source.bounds

        # sk_grid = pv.ImageData()

        # sk_grid.origin = (bounds[0], bounds[1], bounds[2])

        # sk_grid.spacing = (
        #     abs(bounds[0]-bounds[1]) / resolution[0],
        #     abs(bounds[2]-bounds[3]) / resolution[1],
        #     abs(bounds[4]-bounds[5]) / resolution[2]
        # )

        # sk_grid.dimensions = resolution[0] + 1, resolution[1] + 1, resolution[2] + 1

        # sampled = sk_grid.sample(source)

        # volume = np.array(source.points)

        # print(type(volume), volume.ndim, volume.shape)

        # volume = source.cast_to_poly_points().delaunay_3d()

        # if volume is not None:
        #     verts, faces, _, _ = measure.marching_cubes(volume, level=0.5)

        #     source = pv.PolyData(verts, faces)

        # resampled = None

        # if CONTOUR:
        #     s_data = pv.create_grid(source, dimensions=(
        #         abs(source.bounds[0]-source.bounds[1]),
        #         abs(source.bounds[2]-source.bounds[3]),
        #         abs(source.bounds[4]-source.bounds[5])
        #     ))

        #     resampled = s_data.sample(source)

        #     if resampled is not None:
        #         print(resampled.point_data.keys())

        #         source = resampled.contour(1, scalars=resampled['RegionId'],
        #                                 method='marching_cubes', progress_bar=True)
                
        #         print(source)

        #     sourcex = source.threshold([ref_bounds[0], ref_bounds[1]], scalars='x')
        #     sourcey = sourcex.threshold([ref_bounds[2], ref_bounds[3]], scalars='y') # type: ignore
        #     sourcez = sourcey.threshold([ref_bounds[4], ref_bounds[5]], scalars='z') # type: ignore

        # source_mesh = sourcez

        # if sourcez is not None:
        #     source_mesh = sourcez.threshold(-1e4) # type: ignore


        # Widgets

        t_slider = QSlider(Qt.Orientation.Horizontal, self)
        t_slider.setRange(0, 100)
        t_slider.setValue(40)
        t_slider.setSingleStep(5)
        t_slider.setPageStep(10)
        t_slider.setTickPosition(QSlider.TickPosition.TicksAbove)
        t_slider.valueChanged.connect(self.update_threshold_val)
        t_slider.sliderReleased.connect(self.threshold_change)

        self.t_slider = t_slider

        layout.addWidget(self.t_slider)

        self.submit_button = QPushButton("Submit")
        self.submit_button.clicked.connect(self.submit)

        layout.addWidget(self.submit_button)

        self.threshold_label = QLabel("* -2e4", self)

        layout.addWidget(self.threshold_label)


        # Visualize meshes

        box = pv.Box(self.ref_bounds)

        self.plotter: pv.Plotter = pvqt.QtInteractor(self) # type: ignore

        self.tgt_actor = self.plotter.add_mesh(self.tgt_mesh, color='lightgray')

        self.ref_actor = self.plotter.add_mesh(self.ref_mesh, opacity=0.4)

        if PLOT_BOX:
            _ = self.plotter.add_mesh(box, opacity=0.6)

        # plotter.add_mesh_threshold(mesh, pointa=(-2e4,1), pointb=(0,1))

        self.plotter.camera_position = 'xy'

        layout.addWidget(self.plotter) # type: ignore

        self.plotter.show()
    
    def update_threshold_val(self, value):
        self.threshold_val = value
    
    def threshold_change(self):
        self.threshold_label.setText(f"{self.threshold_val} * -2e4")
        
        self.threshold((self.threshold_val / 10) * -2e4)

        if self.tgt_mesh is None:
            exit("Source could not be thresholded")
        
        self.rerender('tgt')
    
    def threshold(self, val, mesh=None):
        mesh = self.get_mesh(mesh)

        self.tgt_mesh = mesh.threshold(val, progress_bar=True)

        if self.tgt_mesh is None:
            exit("Source could not be thresholded")
    
    def clip_box(self, mesh=None):
        mesh = self.get_mesh(mesh)

        self.tgt_mesh = mesh.clip_box(self.ref_bounds, invert=False,
                                      progress_bar=True)
    
        if self.tgt_mesh is None:
            exit("Target mesh could not be clipped")

    def find_largest(self, mesh=None):
        mesh = self.get_mesh(mesh)

        self.tgt_mesh = mesh.connectivity('largest', progress_bar=True)

        if self.tgt_mesh is None:
            exit("Largest connectivity extraction failed")
    
    def get_mesh(self, mesh):
        mesh = self.alpha_mesh if mesh is None else mesh

        if mesh is None:
            exit("Target mesh is undefined")
        
        return mesh
    
    def rerender(self, actor=None):
        if actor is None:
            self.plotter.clear()
        else:
            if actor == 'tgt':
                self.plotter.remove_actor(self.tgt_actor)

                self.tgt_actor = self.plotter.add_mesh(self.tgt_mesh)
            else:
                self.plotter.clear()
        
        self.plotter.show()

    def submit(self):
        self.threshold(self.threshold_val)
        self.clip_box()
        self.find_largest()




qt_app = QApplication([])

window = MainWindow()

window.show()

qt_app.exec_()

