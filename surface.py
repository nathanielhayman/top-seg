import pyvista as pv

import numpy as np

from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget, \
QSlider, QPushButton, QLabel, QAction, QHBoxLayout, QGridLayout, QCheckBox, \
QFileDialog, QInputDialog

from PyQt5.QtGui import QIcon

from PyQt5.QtCore import Qt

from PyQt5.QtCore import QObject

import pyvistaqt as pvqt

from vtkmodules.vtkCommonDataModel import vtkIterativeClosestPointTransform
from vtkmodules.vtkFiltersCore import vtkImplicitPolyDataDistance

import vtk

# from vmtk import pypes

FILENAME = "./data/topcow_mr_001_0000.nii.gz"

REF_FILENAME = "./data/topcow_mr_001.nii.gz"

PLOT_BOX = False

CONTOUR = False

DEFAULT_THRESHOLD = 2.1

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        # Initialize Qt window

        self.setWindowTitle("TopSeg Visualization")
        # self.setGeometry(100, 100, 100, 100)

        self.setWindowIcon(QIcon("icon.png"))

        page_layout = QVBoxLayout()

        slider_layout = QHBoxLayout()
        bsize_layout = QHBoxLayout()
        id_layout = QHBoxLayout()
        subpage_layout = QHBoxLayout()

        sidebar_layout = QVBoxLayout()
        mesh_layout = QVBoxLayout()

        sidebar_layout.setSpacing(0)
        sidebar_layout.setContentsMargins(0, 0, 0, 0)

        subpage_layout.addLayout(sidebar_layout)
        subpage_layout.addLayout(mesh_layout)

        page_layout.addLayout(slider_layout)
        page_layout.addLayout(bsize_layout)
        page_layout.addLayout(id_layout)
        page_layout.addLayout(subpage_layout)

        # central_widget = QWidget()
        # self.setCentralWidget(central_widget)

        # mesh_layout = QHBoxLayout()
        # mesh_layout.setStretch(0, 40)
        # mesh_layout.setStretch(1, 200)
        # central_widget.setLayout(mesh_layout)


        # Read reference mesh

        reader = pv.get_reader(REF_FILENAME)
        mesh = reader.read()

        self.ref_mesh = mesh.threshold(1)

        geo_filter = vtk.vtkGeometryFilter()
        geo_filter.SetInputData(self.ref_mesh)
        geo_filter.Update()

        poly_data = geo_filter.GetOutput()

        smooth_filter = vtk.vtkSmoothPolyDataFilter()
        smooth_filter.SetInputData(poly_data)
        smooth_filter.SetNumberOfIterations(30)
        smooth_filter.SetRelaxationFactor(0.1)
        smooth_filter.FeatureEdgeSmoothingOff()
        smooth_filter.BoundarySmoothingOn()  
        smooth_filter.Update()

        self.ref_mesh = pv.wrap(smooth_filter.GetOutput())

        # Read target mesh

        self.tgt_mesh = None

        self.tgt_actor = None

        self.ref_actor = None

        self.alpha_mesh = None

        self.reset_mesh = None

        self.last_mesh = None

        self.box_actor = None

        self.box_mesh = None

        self.ref_mesh_arr = []

        self.ref_mesh_arr.append(self.ref_mesh)

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
        self.threshold_val = DEFAULT_THRESHOLD * 10

        # resolution = (1000, 1000, 1000)

        # print(resolution)

        # print((resolution[0] + 1, resolution[1] + 1,
        #                                 resolution[2] + 1))

        # print((resolution[0] + 1, resolution[1] + 1,
        #                                 resolution[2] + 1))

        # if reference is None:
        #     exit("Could not construct reference mesh")

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


        # Menu

        main_menu = self.menuBar()

        file_menu = None
        edit_menu = None
        view_menu = None

        if main_menu is not None:
            file_menu = main_menu.addMenu("File")
            edit_menu = main_menu.addMenu("Edit")
            view_menu = main_menu.addMenu("View")

        clear_btn = QAction("Reset values", self)
        clear_btn.setShortcut("Ctrl+R")
        clear_btn.triggered.connect(self.action_reset)

        isolate_btn = QAction("Isolate CoW region", self)
        isolate_btn.setShortcut("Ctrl+E")
        isolate_btn.triggered.connect(self.action_isolate)

        ref_btn = QAction("Show reference CoW", self)
        ref_btn.setShortcut("Ctrl+C")
        ref_btn.triggered.connect(self.action_show_ref)

        undo_btn = QAction("Undo", self)
        undo_btn.setShortcut("Ctrl+Z")
        undo_btn.triggered.connect(self.action_undo)

        upload_btn = QAction("Upload target scan", self)
        upload_btn.setShortcut("Ctrl+Q")
        upload_btn.triggered.connect(self.action_upload_target)

        save_btn = QAction("Save as", self)
        save_btn.setShortcut("Ctrl+S")
        save_btn.triggered.connect(self.action_save)

        refup_btn = QAction("Upload reference mesh", self)
        refup_btn.setShortcut("Ctrl+W")
        refup_btn.triggered.connect(self.action_upload_ref)

        ref1_btn = QAction("Switch to ref 1", self)
        ref1_btn.setShortcut("Ctrl+1")
        ref1_btn.triggered.connect(lambda _: self.set_ref_mesh(0))

        ba_btn = QAction("Expand region on all", self)
        ba_btn.setShortcut("Ctrl+H")
        ba_btn.triggered.connect(lambda _: self.set_btype("bsize", "all"))

        bx_btn = QAction("Expand region on X", self)
        bx_btn.setShortcut("Ctrl+J")
        bx_btn.triggered.connect(lambda _: self.set_btype("bsize", "x"))

        by_btn = QAction("Expand region on Y", self)
        by_btn.setShortcut("Ctrl+K")
        by_btn.triggered.connect(lambda _: self.set_btype("bsize", "y"))

        bz_btn = QAction("Expand region on Z", self)
        bz_btn.setShortcut("Ctrl+L")
        bz_btn.triggered.connect(lambda _: self.set_btype("bsize", "z"))

        px_btn = QAction("Shift region on X", self)
        px_btn.setShortcut("Ctrl+U")
        px_btn.triggered.connect(lambda _: self.set_btype("bpos", "x"))

        py_btn = QAction("Shift region on Y", self)
        py_btn.setShortcut("Ctrl+I")
        py_btn.triggered.connect(lambda _: self.set_btype("bpos", "y"))

        pz_btn = QAction("Shift region on Z", self)
        pz_btn.setShortcut("Ctrl+O")
        pz_btn.triggered.connect(lambda _: self.set_btype("bpos", "z"))

        h_btn = QAction("Compute Hausdorff Distance", self)
        h_btn.setShortcut("Ctrl+D")
        h_btn.triggered.connect(lambda _: self.compute_hausdorff())

        if edit_menu is not None:
            edit_menu.addAction(clear_btn)
            edit_menu.addAction(isolate_btn)
            edit_menu.addAction(ref_btn)
            edit_menu.addAction(undo_btn)
            edit_menu.addSeparator()
            edit_menu.addAction(ba_btn)
            edit_menu.addAction(bx_btn)
            edit_menu.addAction(by_btn)
            edit_menu.addAction(bz_btn)
            edit_menu.addSeparator()
            edit_menu.addAction(px_btn)
            edit_menu.addAction(py_btn)
            edit_menu.addAction(pz_btn)
        
        if file_menu is not None:
            file_menu.addAction(upload_btn)
            file_menu.addAction(refup_btn)
            file_menu.addAction(save_btn)

        if view_menu is not None:
            view_menu.addAction(h_btn)
            view_menu.addAction(ref1_btn)
        
        self.view_menu = view_menu


        # Widgets

        t_label = QLabel("Threshold value: ", self)

        slider_layout.addWidget(t_label)

        t_slider = QSlider(Qt.Orientation.Horizontal, self)
        t_slider.setRange(0, 30)
        t_slider.setValue(21)
        t_slider.setSingleStep(1)
        t_slider.setPageStep(5)
        t_slider.setTickPosition(QSlider.TickPosition.TicksAbove)
        t_slider.valueChanged.connect(self.update_threshold_val)
        t_slider.sliderReleased.connect(self.threshold_change)

        self.t_slider = t_slider

        slider_layout.addWidget(self.t_slider)

        self.threshold_label = QLabel(f"{DEFAULT_THRESHOLD} * -1e4", self)

        slider_layout.addWidget(self.threshold_label)

        self.b_type = "all"

        self.b_label = QLabel(f"Isolation box expansion on all: ", self)

        bsize_layout.addWidget(self.b_label)

        self.bsize = {
            "x": 0,
            "y": 0,
            "z": 0
        }

        self.bpos = {
            "x": 0,
            "y": 0,
            "z": 0
        }
        
        self.bchange = "bsize"

        self.bshow = False

        b_slider = QSlider(Qt.Orientation.Horizontal, self)
        b_slider.setRange(-50, 50)
        b_slider.setValue(self.bsize["x"])
        b_slider.setSingleStep(1)
        b_slider.setPageStep(5)
        b_slider.setTickPosition(QSlider.TickPosition.TicksAbove)
        b_slider.valueChanged.connect(self.update_bsize_val)
        b_slider.sliderReleased.connect(self.bsize_val_change)

        self.b_slider = b_slider

        bsize_layout.addWidget(b_slider)

        self.bsize_text = QLabel("+0x", self)

        bsize_layout.addWidget(self.bsize_text)

        i_label = QLabel("Identification distance: ", self)

        id_layout.addWidget(i_label)

        self.id_distance = 4

        i_slider = QSlider(Qt.Orientation.Horizontal, self)
        i_slider.setRange(0, 500)
        i_slider.setValue(int(self.id_distance * 10))
        i_slider.setSingleStep(10)
        i_slider.setPageStep(20)
        i_slider.setTickPosition(QSlider.TickPosition.TicksAbove)
        i_slider.valueChanged.connect(self.update_id_val)

        id_layout.addWidget(i_slider)

        self.id_text = QLabel(f"{self.id_distance}px", self)

        id_layout.addWidget(self.id_text)

        self.upload_btn = QPushButton("Upload scan")
        self.upload_btn.clicked.connect(self.action_upload_target)

        sidebar_layout.addWidget(self.upload_btn)

        self.submit_btn = QPushButton("Submit")
        self.submit_btn.clicked.connect(self.submit)

        sidebar_layout.addWidget(self.submit_btn)

        self.largest_btn = QPushButton("Largest component")
        self.largest_btn.clicked.connect(self.action_get_largest)

        sidebar_layout.addWidget(self.largest_btn)

        self.align_btn = QPushButton("Align")
        self.align_btn.clicked.connect(self.align)

        sidebar_layout.addWidget(self.align_btn)

        self.remove_btn = QPushButton("Remove furthest")
        self.remove_btn.clicked.connect(self.remove_furthest)

        sidebar_layout.addWidget(self.remove_btn)

        self.contour_btn = QPushButton("Smooth")
        self.contour_btn.clicked.connect(self.contour)

        sidebar_layout.addWidget(self.contour_btn)

        self.isolate_btn = QPushButton("Isolate CoW region")
        self.isolate_btn.clicked.connect(self.action_isolate)

        sidebar_layout.addWidget(self.isolate_btn)

        self.ref_check = QCheckBox("Show reference CoW")
        self.ref_check.setChecked(False)
        self.ref_check.stateChanged.connect(self.action_show_ref)

        sidebar_layout.addWidget(self.ref_check)

        self.iso_check = QCheckBox("Show isolation region")
        self.iso_check.setChecked(False)
        self.iso_check.stateChanged.connect(self.action_show_iso)

        sidebar_layout.addWidget(self.iso_check)


        # Visualize meshes

        self.plotter: pv.Plotter = pvqt.QtInteractor(self) # type: ignore

        # self.tgt_actor = self.plotter.add_mesh(self.tgt_mesh, color='lightgray')

        # plotter.add_mesh_threshold(mesh, pointa=(-2e4,1), pointb=(0,1))

        # self.ref_actor = self.plotter.add_mesh(self.ref_mesh, opacity=0.4,
        #                                            color='purple')

        self.plotter.camera_position = 'xy'

        mesh_layout.addWidget(self.plotter) # type: ignore

        widget = QWidget()
        widget.setLayout(page_layout)
        self.setCentralWidget(widget)

        self.plotter.show()
    
    def update_threshold_val(self, value):
        self.threshold_val = value / 10
        self.threshold_label.setText(f"{self.threshold_val} * -1e4")
    
    def threshold_change(self):
        self.iso_check.setChecked(False)
        self.threshold((self.threshold_val) * -1e4)

        if self.tgt_mesh is None:
            exit("Source could not be thresholded")
        
        self.rerender(self.tgt_mesh)
    
    def update_bsize_val(self, value):
        if self.bchange == "bsize":
            if self.b_type == "all":
                self.bsize["x"] = value / 50
                self.bsize["y"] = value / 50
                self.bsize["z"] = value / 50
                self.bsize_text.setText(f"+{self.bsize['x']}x")
            else:
                self.bsize[self.b_type] = value / 50
                self.bsize_text.setText(f"+{self.bsize[self.b_type]}x")
        elif self.bchange == "bpos":
            self.bpos[self.b_type] = value / 5
            self.bsize_text.setText(f"+{self.bpos[self.b_type]}x")

    def update_id_val(self, value):
        self.id_distance = value / 10
        self.id_text.setText(f"+{self.id_distance}px")

    def set_btype(self, change="bsize", val="all"):
        self.bchange = change
        self.b_type = val

        if val == "all":
            self.bsize["y"] = self.bsize["x"]
            self.bsize["z"] = self.bsize["x"]

        if change == "bsize":
            self.b_slider.setValue(int(self.bsize[val] * 50 if val != "all" 
                                else self.bsize["x"] * 50))
        else:
            self.b_slider.setValue(int(self.bpos[val] * 10))

        self.b_label.setText(
            f"Isolation box {'expansion' if change == 'bsize' else 'shift'} on {val}: "
        )
    
    def bsize_val_change(self):
        self.box_mesh = pv.Box(self.get_bsize())
        
        if self.box_actor is not None:
            self.plotter.remove_actor(self.box_actor)
        
        if self.bshow:
            self.box_actor = self.plotter.add_mesh(self.box_mesh, opacity=0.2, 
                                                color='lightgreen')

        self.plotter.show()
    
    def get_bsize(self):
        b = self.ref_bounds

        if b is None:
            return (0, 0, 0, 0, 0, 0)

        dx = self.bsize["x"]*abs(b[1] - b[0])
        dy = self.bsize["y"]*abs(b[3] - b[2])
        dz = self.bsize["z"]*abs(b[5] - b[4])

        px = self.bpos["x"]
        py = self.bpos["y"]
        pz = self.bpos["z"]

        return (b[0]-dx+px, b[1]+dx+px, 
                b[2]-dy+py, b[3]+dy+py, 
                b[4]-dz+pz, b[5]+dz+pz)
    
    def threshold(self, val, mesh=None):
        mesh = self.get_mesh(mesh)

        self.tgt_mesh = mesh.threshold(val, progress_bar=True)

        if self.tgt_mesh is None:
            exit("Source could not be thresholded")
    
    def clip_box(self, mesh=None):
        mesh = self.get_mesh(mesh)
        
        nb = self.get_bsize()

        self.tgt_mesh = mesh.clip_box(nb, invert=False)
    
        if self.tgt_mesh is None:
            exit("Target mesh could not be clipped")

    def find_largest(self, mesh=None, n=1):
        mesh = self.get_mesh(mesh)

        # geo_filter = vtk.vtkGeometryFilter()
        # geo_filter.SetInputData(mesh)
        # geo_filter.Update()

        # data = geo_filter.GetOutput()

        # connectivity_filter = vtk.vtkConnectivityFilter()
        # connectivity_filter.SetInputData(mesh)
        # connectivity_filter.SetExtractionModeToLargestRegion()
        # connectivity_filter.Update()

        # self.tgt_mesh = pv.wrap(connectivity_filter.GetOutput())

        self.tgt_mesh = mesh.connectivity('specified', range(n))

        if self.tgt_mesh is None:
            exit("Largest connectivity extraction failed")
    
    def get_mesh(self, mesh):
        mesh = self.alpha_mesh if mesh is None else mesh

        if mesh is None:
            exit("Target mesh is undefined")
        
        return mesh
    
    def rerender(self, mesh=None, ref=False):
        mesh = self.get_mesh(mesh)
        
        if ref:
            if self.ref_actor:
                self.plotter.remove_actor(self.ref_actor)
            
            self.ref_actor = self.plotter.add_mesh(mesh, color='purple',
                                                   opacity=0.4)
        else:
            if self.tgt_actor:
                self.plotter.remove_actor(self.tgt_actor)

            self.tgt_actor = self.plotter.add_mesh(mesh, color='lightgray')

        self.plotter.set_focus(mesh.center)

        self.plotter.show()
    
    def action_reset(self):
        self.threshold_label.setText(f"{DEFAULT_THRESHOLD} * -1e4")
        self.t_slider.setValue(int(DEFAULT_THRESHOLD * 10))
        self.iso_check.setChecked(False)

        self.tgt_mesh = self.reset_mesh

        self.rerender(self.tgt_mesh)
    
    def action_isolate(self):
        self.last_mesh = self.tgt_mesh

        self.clip_box(self.tgt_mesh)

        self.rerender(self.tgt_mesh)
    
    def action_show_ref(self):
        if self.ref_actor is not None:
            self.plotter.remove_actor(self.ref_actor)

            self.plotter.show()

            self.ref_actor = None
        else:
            self.ref_actor = self.plotter.add_mesh(self.ref_mesh, opacity=0.4,
                                                   color='purple')

            self.plotter.show()

    def action_show_iso(self):
        self.bshow = not self.bshow

        self.bsize_val_change()
    
    def action_get_largest(self):
        self.last_mesh = self.tgt_mesh

        text, ok = QInputDialog().getText(self, "Find largest components", 
                                          "How many components?")

        if ok and text:
            self.find_largest(self.tgt_mesh, int(text))

            self.rerender(self.tgt_mesh)
    
    def action_upload_target(self):
        dialog = QFileDialog(self, caption="Upload target mesh",
                                filter="Scans (*.nii, *.nii.gz)")
        dialog.setFileMode(QFileDialog.FileMode.ExistingFiles)
        dialog.setViewMode(QFileDialog.ViewMode.List)
        
        if dialog.exec():
            filenames = dialog.selectedFiles()

            print(filenames[0])

            reader = pv.get_reader(filenames[0])
            self.tgt_mesh = reader.read()

            self.alpha_mesh = reader.read()

            self.threshold(DEFAULT_THRESHOLD * -1e4)

            self.reset_mesh = self.tgt_mesh

            self.last_mesh = self.tgt_mesh

            self.rerender(self.tgt_mesh)
    
    def action_undo(self):
        self.tgt_mesh = self.last_mesh

        self.rerender(self.tgt_mesh)
    
    def action_save(self):
        file_name, _ = QFileDialog.getSaveFileName(self, "Save Mesh", "", "OBJ file (*.obj)")
        
        if file_name:
            pv.save_meshio(file_name, self.tgt_mesh, file_format="obj")
    
    def set_ref_mesh(self, i):
        self.ref_mesh = self.ref_mesh_arr[i]

        self.rerender(self.ref_mesh, ref=True)
    
    def action_upload_ref(self):
        dialog = QFileDialog(self, caption="Upload target mesh",
                                filter="Scans (*.nii, *.nii.gz)")
        dialog.setFileMode(QFileDialog.FileMode.ExistingFiles)
        dialog.setViewMode(QFileDialog.ViewMode.List)

        if dialog.exec():
            filenames = dialog.selectedFiles()

            reader = pv.get_reader(filenames[0])
            mesh = reader.read()

            mesh = mesh.threshold(1)

            self.ref_mesh_arr.append(mesh)

            self.ref_mesh = mesh

            size = len(self.ref_mesh_arr)

            ref_btn = QAction(f"Switch to ref {size}", self)
            ref_btn.setShortcut(f"Ctrl+{size}")
            ref_btn.triggered.connect(lambda _: self.set_ref_mesh(size-1))

            if self.view_menu is not None:
                self.view_menu.addAction(ref_btn)

    def remove_furthest(self):
        mesh = self.get_mesh(self.tgt_mesh)

        self.last_mesh = mesh

        distances = mesh.compute_implicit_distance(self.ref_mesh)

        mesh['distance'] = distances['implicit_distance']

        self.tgt_mesh = mesh.threshold(
            [-self.id_distance, self.id_distance],
            scalars='distance'
        )

        self.rerender(self.tgt_mesh)

    def contour(self):
        mesh = self.get_mesh(self.tgt_mesh)

        geo_filter = vtk.vtkGeometryFilter()
        geo_filter.SetInputData(mesh)
        geo_filter.Update()

        poly_data = geo_filter.GetOutput()

        smooth_filter = vtk.vtkSmoothPolyDataFilter()
        smooth_filter.SetInputData(poly_data)
        smooth_filter.SetNumberOfIterations(30)
        smooth_filter.SetRelaxationFactor(0.1)
        smooth_filter.FeatureEdgeSmoothingOff()
        smooth_filter.BoundarySmoothingOn()  
        smooth_filter.Update()

        self.tgt_mesh = pv.wrap(smooth_filter.GetOutput())

        self.rerender(self.tgt_mesh)

    def align(self):
        icp = vtkIterativeClosestPointTransform()
        icp.SetSource(self.ref_mesh)
        icp.SetTarget(self.tgt_mesh)
        # icp.GetLandmarkTransform().SetModeToRigidBody()
        icp.SetMaximumNumberOfLandmarks(1000)
        # icp.SetMaximumMeanDistance(.00001)
        icp.SetMaximumNumberOfIterations(1000)
        # icp.CheckMeanDistanceOn()
        # icp.StartByMatchingCentroidsOn()
        icp.Update()

        self.ref_mesh = self.ref_mesh.transform(icp.GetMatrix())

        self.rerender(self.ref_mesh, ref=True)
    
    def compute_hausdorff(self):
        hausdorff_filter = vtk.vtkHausdorffDistancePointSetFilter()
        hausdorff_filter.SetInputData(0, self.tgt_mesh)
        hausdorff_filter.SetInputData(1, self.ref_mesh)
        hausdorff_filter.Update()

        hd = hausdorff_filter.GetOutput().GetFieldData().GetArray(0).GetValue(0)

        print(f"Hausdorff distance: {hd}")

    def submit(self):
        self.threshold(self.threshold_val * -1e4)
        self.clip_box(self.tgt_mesh)
        self.find_largest(self.tgt_mesh)
        self.align()
        self.remove_furthest()

        # geo_filter = vtk.vtkGeometryFilter()
        # geo_filter.SetInputData(self.tgt_mesh)
        # geo_filter.Update()

        # poly_data = geo_filter.GetOutput()

        # writer = vtk.vtkXMLPolyDataWriter()
        # writer.SetFileName("artery.vtp")
        # writer.SetInputData(poly_data)
        # writer.Write()

        # my_pype = pypes.PypeRun("vmtknetworkextraction -ifile artery.vtp -advancementratio 1 -ofile centerlines/artery.vtp")
        
        # reader = vtk.vtkOBJReader()
        # reader.SetFileName("artery.vtp")
        # reader.Update()
        # data = reader.GetOutput()

        # print(data)

        self.rerender(self.tgt_mesh)


qt_app = QApplication([])

window = MainWindow()

window.show()

qt_app.exec_()

