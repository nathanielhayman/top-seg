# TopSeg

Main file: `surface.py`

Skeleton analysis library: `remove.c`

C-Python interface: `cutil/remove.py`

Testing for skeleton analysis: `testing.py`

## Operations & Shortcuts

### File operations
- Upload target scan (`Ctrl+Q`): Opens prompt to select a `.nii` or `.nii.gz` file to use as the target scan for performing analysis
- Upload reference mesh (`Ctrl+W`): Opens prompt to select a `.nii` or `.nii.gz` file to use as a refernce scan for performing analysis (this can be switched with other uploaded references in the "View" tab)
- Save as (`Ctrl+S`): Save the current, updated target mesh as a `.obj` file

### Edit operations
- Reset values (`Ctrl+R`): Resets the target mesh to its original state after it was uploaded
- Isolate CoW Region (`Ctrl+E`): Clips the mesh by a bounding box (select "Show isolation region" to see the bounding box)
- Show reference CoW (`Ctrl+C`): Reveals an overlay of an "ideal" CoW to reference during identification operations. This corresponds to the currently selected reference (see the "View" subsection for details)
- Undo (`Ctrl+Z`): Undoes the last operation (does not stack)
- Expand region on all/X/Y/Z (`Ctrl+H/J/K/L`): Sets the "Isolation box" slider to grow the isolation region on the specified axis
- Shift region on X/Y/Z (`Ctrl+U/I/O`): Sets the "Isolation box" slider to shift the isolation region on the specified axis

### View operations
- Compute Hausdorff Distance (`Ctrl+D`): Computes the Hausdorff Distance between the target and reference meshes. The result is logged to the console
- Switch to ref `x` (`Ctrl+<x>`): Activates the uploaded reference mesh stored in slot `x`

### Misc buttons
- Largest component: Opens a menu to take the largest `n` connected components from the target mesh
- Smooth: Performs Laplacian smoothing on the target mesh
- Align: Aligns the reference mesh with the target mesh via Laplacian deformation
- Remove furthest: Removes points on the target mesh greater than the distance specified by the "Identification distance" slider from the reference mesh
- Submit: Performs each of the described analysis operations in the following sequence: threshold, isolate region, take the largest connected component, align the meshes, remove furthest points from the reference mesh