import pyvista as pv

REF_FILENAME = "./data/topcow_mr_001.nii.gz"

reader = pv.get_reader(REF_FILENAME)

mesh = reader.read()

thresholded = mesh.threshold(1)

thresholded.plot()
