# Examples

Every example is a plain Python script using the same `opengrid.Session` API as the desktop
application. Run one from the repository root:

```bash
conda activate opengrid
python examples/12a_airfoil_cases.py      # writes its figure, projects and meshes to examples/output/
python -m opengrid examples/output/12a_cascade_2d.ogp   # open a result in the application
```

In the application, the Script tab shows everything done in the window as the same kind of
Python, and runs new lines written at its end: a way to turn clicks into a script like these.

The projects (`.ogp`) and journals (`*_journal.py`, the replayable script of what was done)
of each run are in [`output/`](output); the meshes (CGNS, Gmsh, VTKHDF) are written there
too but not committed. Several scripts take `--gui` to open their model in the application.

## Start here: complete cases

| Script | What it builds |
|---|---|
| [`12a_airfoil_cases.py`](12a_airfoil_cases.py) | A NACA 2412 in a far field and a periodic NACA 6412 cascade from the standard shapes: clustered connectors, an unstructured domain with 40 inflation layers, boundary types in one step, a periodic pair, the cascade extruded one cell thick for a 3D solver, CGNS and Gmsh |
| [`12c_sliding_layers.py`](12c_sliding_layers.py) | Inflation layers set from the curves (`set_inflation`), tapering or running into the neighbouring curves (`slide`), and a wedge too sharp for them |
| [`12b_bodies_in_a_box.py`](12b_bodies_in_a_box.py) | A sphere and a closed swept wing in boxes: surface meshes with shared edges (`mesh_surfaces`), tetrahedra with prism layers and trailing-edge fans, wall and far-field patches, CGNS and Gmsh |
| [`04_solver_ready.py`](04_solver_ready.py) | A backward-facing step (three blocks) and a NACA 0012 C-grid with patches and zones, exported to CGNS (checked by `cgnscheck` when installed), Gmsh (read back by meshio) and VTKHDF (read back by VTK) |

## By topic

| Topic | Scripts |
|---|---|
| Geometry, projects, distributions | [`01_geometry.py`](01_geometry.py), [`02_connectors.py`](02_connectors.py) |
| Structured 2D grids | [`03_airfoil_grids.py`](03_airfoil_grids.py) (O- and C-grids with wall spacing), [`03_backward_step.py`](03_backward_step.py) (three blocks), [`11a_wall_layers.py`](11a_wall_layers.py) (orthogonal wall layers) |
| Unstructured 2D meshes | [`05a_unstructured.py`](05a_unstructured.py) (Delaunay), [`05b_advancing_front.py`](05b_advancing_front.py), [`05c_quads.py`](05c_quads.py) (quad-dominant, all-quad), [`06_inflation.py`](06_inflation.py) (boundary layers), [`12c_sliding_layers.py`](12c_sliding_layers.py) (layers from curves, sliding) |
| Structured 3D blocks | [`07a_blocks.py`](07a_blocks.py), [`07b_smoothing.py`](07b_smoothing.py) (hyperbolic extrusion, 3D elliptic smoothing) |
| Unstructured 3D volumes | [`08a_volumes.py`](08a_volumes.py) (tetrahedra), [`08b_inflation.py`](08b_inflation.py) (prism / hexahedron layers), [`08c_optimisation_mixed.py`](08c_optimisation_mixed.py) (fans, optimisation, mixed models), [`12b_bodies_in_a_box.py`](12b_bodies_in_a_box.py) |
| Multi-block topology and smoothing groups | [`09b_topology_golden.py`](09b_topology_golden.py), [`10a_multiblock_smoothing.py`](10a_multiblock_smoothing.py), [`10b_multiblock_3d.py`](10b_multiblock_3d.py), [`10e_templates_auto.py`](10e_templates_auto.py) (templates, automatic topology) |
| Interactive building blocks | [`10c_freehand_live.py`](10c_freehand_live.py) (freehand curves, live dragging), [`10d_sketch_topology.py`](10d_sketch_topology.py) (sketch topology around CAD), [`00_launch.py`](00_launch.py) (opens the application) |
| CAD import | [`09a_cad_import.py`](09a_cad_import.py) (STEP / IGES through OpenCASCADE) |
| High-order elements | [`11e_high_order.py`](11e_high_order.py), [`11f_high_order_elements.py`](11f_high_order_elements.py), [`11g_curved_elements.py`](11g_curved_elements.py) |
| Performance | [`09c_performance.py`](09c_performance.py), [`11b_performance.py`](11b_performance.py), [`11d_gpu.py`](11d_gpu.py) (experimental OpenCL) |

## Golden cases

[`golden/`](golden) holds the ten validation cases of the specification (unit square, quarter
annulus, NACA 0012 grids, cylinder O-grid, backward step, hybrid airfoil, duct, annular
sector, sphere in a box, wing stub). Each script builds its case from nothing; the test suite
runs them all and checks their metrics against reference values (`tests/python/test_golden.py`).

The numbered scripts are named after the development phase that introduced the feature; the
scripting guide [`docs/SCRIPTING.md`](../docs/SCRIPTING.md) explains the calls they use.
