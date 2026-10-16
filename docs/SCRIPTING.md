# Scripting OpenGrid

Everything the GUI does goes through one Python object, `opengrid.Session`; scripts and the
GUI's console use the same API. The full reference is `docs/API.md` (generated from the
docstrings); complete scripts for the golden validation cases are in `examples/golden/`.

```python
import opengrid

og = opengrid.Session()
og.create_line((0, 0), (1, 0), name="bottom")
og.create_connector("bottom", 21, {"type": "tanh", "begin_spacing": 1e-3}, name="c-bottom")
print(og.connector_report("c-bottom"))
```

## The model

Entities have a unique name and a kind; methods take names (or ids). They form layers, each
referring to the one below:

| Layer | Kinds | Created with |
|---|---|---|
| Geometry | `point`, `line`, `arc`, `polyline`, `spline`, `nurbs`; surfaces `plane`, `ruled_surface`, `extruded_surface`, `revolved_surface`, `nurbs_surface` | `create_line`, `create_spline`, ..., `import_airfoil_dat`, `import_csv`, `import_cad` |
| Connectors | `connector`: a curve (or chain of curves) with a number of points and a distribution | `create_connector`, `create_chain_connector`, `split_connector`, `join_connectors` |
| Domains | `domain` (structured, 4 edges), `unstructured_domain` (loops; triangles / quads, boundary layers) | `create_structured_domain`, `create_domain`, `create_unstructured_domain` |
| 3D | `block` (structured: six faces or an extrusion), `volume` (closed shell of domains: tetrahedra, boundary layers, pyramids) | `create_block`, `create_extruded_block`, `create_volume` |
| Boundary conditions | `patch` (typed group of boundary connectors / faces), `zone` | `create_patch`, `create_zone`, `pair_periodic` |

Meshes are not stored: they are computed from the parameters when asked for (`domain_grid`,
`unstructured_mesh`, `block_grid`, `volume_mesh`) and cached until something they depend on
changes. The same parameters always give the same mesh (deterministic algorithms).

**Parameters are validated, never repaired.** Invalid input raises `opengrid.ValidationError`
whose `issues` name the entity, the field and what is wrong; the model is left unchanged.

```python
try:
    og.create_domain([["c-bottom"], ["c-right"], ["c-top"], ["c-left"]])
except opengrid.ValidationError as error:
    for issue in error.issues:
        print(issue.where, issue.message)
```

## Commands, undo and journals

Every change (create, modify, rename, delete, split, import, ...) is one command:

- `og.undo()` / `og.redo()` restore exact before / after states (nothing is recomputed);
- `og.journal_text()` is a Python script that rebuilds the session from scratch, including
  undo / redo; `opengrid.replay(text)` runs it. Imports are journaled as the entities they
  created, so a replay does not need the imported file.
- `og.save("model.ogp")` / `og.load(...)` write and read the canonical project file.

The GUI writes the same journal (File → Save journal), so an interactive session can be turned
into a script and edited.

## Typical workflows

### 2D structured grid with wall spacing
```python
og.create_circle((0, 0), 0.5, name="cylinder")
og.create_circle((0, 0), 10.0, name="far")
og.create_line((0.5, 0), (10, 0), name="cut")
og.create_connector("cylinder", 129, name="c-cyl")
og.create_connector("far", 129, name="c-far")
og.create_connector(
    "cut", 81, {"type": "tanh", "begin_spacing": 1e-4, "end_spacing": 1.0}, name="c-cut"
)
og.create_domain([[("c-cyl", True)], ["c-cut"], ["c-far"], [("c-cut", True)]], name="ogrid")
og.create_patch("wall", "wall", ["c-cyl"])
og.create_patch("farfield", "farfield", ["c-far"])
og.export_cgns("cylinder.cgns")
```
Elliptic smoothing is accelerated (Anderson, `"acceleration": 5` by default; `0` for plain
relaxation) and runs in parallel.

Edges of a structured domain are listed counter-clockwise: edge 0 (j = 0, i increasing),
edge 1 (i = imax), edge 2 (j = jmax, i decreasing), edge 3 (i = 0); each edge is a chain of
connectors, `(name, True)` meaning reversed. `create_structured_domain(connectors=[...])`
assembles the edges automatically. Elliptic smoothing with wall control:
`og.set_domain_generator(name, elliptic={"walls": [{"edge": 0, "spacing": 1e-4}]})` and
`og.smooth_domain(name)`.

### Hybrid 2D mesh with boundary layers
```python
og.create_unstructured_domain(
    connectors=["c-upper", "c-lower", "c-far"],
    generator={
        "inflation": {
            "connectors": ["c-upper", "c-lower"],
            "first_height": 1e-5,
            "growth_rate": 1.2,
            "layers": 40,
        }
    },
    name="flow",
)
```

### 3D: blocks, volumes and boundary layers
- `create_block(faces=[six structured domains])` (any order and orientation) or
  `create_extruded_block(domain, "translate" | "rotate" | "normal" | "hyperbolic", ...)`.
- `create_volume(faces, generator={"max_size": ..., "inflation": {"faces": walls,
  "first_height": ..., "layers": ...}})` fills a closed shell of domains (structured faces of
  blocks included) with tetrahedra, prism / hexahedron layers and pyramids. Needs Netgen.
- Every 3D cell must meet the quality floor (scaled Jacobian >= 0.05, skewness <= 0.95 by
  default, `quality_floor` to change it); exports refuse models with cells below it.

### All-hexahedral meshes: multi-block structured
Unstructured volumes are tetrahedra-based. All-hexahedral meshes are built from structured
blocks, and the topology tools make the usual layouts quick:

```python
og.ogrid_domain("pipe", fraction=0.5, radial_points=9, radial_spacing=1e-3)  # O-H (butterfly)
for name in og.names("domain"):
    og.create_extruded_block(name, "translate", vector=(0, 0, 5), dimension=41)

og.ogrid_block("ball", fraction=0.5)  # 7 blocks: core + 6 shell blocks
og.ogrid_block("box", fraction=0.4, core=False)  # 6 blocks around a cubic hole
og.split_block("duct", "i", 20)  # split at grid index 20 along i
og.join_blocks("duct", "duct-b")
og.split_domain("d", i=10)
og.join_domains("d", "d-b")
og.split_connector("c-wall", 30)
og.join_connectors("c-wall", "c-wall-b")
```
Splits keep the existing nodes (split connectors get explicit distributions, cut lines follow
the old grid lines); split parts are regenerated by their generator from their new
boundaries.

### CAD
```python
result = og.import_cad("part.step", spacing=0.05, domains="auto")
print(result.report())  # healing log, skipped faces, free edges
og.create_volume(result.domains, name="solid")
```

### Multi-block smoothing groups
Several structured domains can be smoothed as one grid: interfaces and block corners move,
grid lines cross interfaces smoothly, walls keep their first-cell height and other boundaries
slide along their curves::

    names = og.ogrid_domain("pipe", fraction=0.55)        # core + 4 domains
    g = og.create_smoothing_group(names)                  # slide="auto", Laplace + walls
    result = og.smooth_group(g)                           # cached; members' grids follow
    print(og.group_status(g))                             # nodes, irregular points, kink
    og.create_smoothing_group(["a", "b"], walls=[{"connector": "c-wall", "spacing": 1e-3}])

Blocks work the same way (3D)::

    g = og.create_smoothing_group(blocks=og.names("block"),
                                  walls=[{"domain": "sphere-zmax", "spacing": 0.01}])
    og.smooth_group(g)          # block grids, and their faces, follow the group

`result.roles` tells what happened to each connector (interface, wall, sliding, fixed,
outside: used by a domain outside the group). A member's own elliptic settings are ignored
while it is in a group.

### Freehand curves and curve networks
    names = og.create_freehand(points)                # pointer samples -> NURBS curve(s)
    domains = og.domains_from_curves(curves)          # crossing curves -> structured domains
    og.move_handle("h1", "points", 6, (1.0, 0.9, 0))  # as a mouse drag: crossings follow

## Threads
Hot loops run in parallel (OpenMP): `opengrid.threads()` reports the count (one per core by
default, or `OMP_NUM_THREADS`), `opengrid.set_threads(n)` changes it. Results never depend on
it: the same model gives the same mesh, bit for bit, on 1 or 10 threads.

## Exports

| Method | Content |
|---|---|
| `export_cgns` | structured zones with 1-to-1 connectivity, unstructured / mixed zones with element sections, BCs from patches |
| `export_vtkhdf` | VTK's HDF5 format: cells, quality fields, patches |
| `export_gmsh` | Gmsh 4.1, physical groups = patches |
| `export_plot3d` | multi-block structured grids |

Exports check the boundary first (`og.boundary_report()`): every boundary entity in exactly
one patch, no unmeshable volume, no cell below the quality floor.

## Golden cases (`examples/golden/`)

| Script | Case |
|---|---|
| `g01_unit_square.py` | unit square, uniform and stretched |
| `g02_quarter_annulus.py` | quarter annulus, TFI exactness |
| `g03_naca0012_grids.py` | NACA 0012 O-grid and C-grid with wall spacing (elliptic) |
| `g04_cylinder_ogrid.py` | cylinder O-grid with wall spacing |
| `g05_backward_step.py` | backward-facing step, three blocks |
| `g06_airfoil_hybrid.py` | airfoil hybrid mesh (layers + triangles) |
| `g07_duct.py` | 3D duct (six-face block and extrusion) |
| `g08_annular_sector.py` | annular sector duct |
| `g09_sphere_in_box.py` | sphere in a box (prisms + tetrahedra) |
| `g10_wing_stub.py` | swept wing stub (hexahedra, prisms, tetrahedra) |

Each defines `build()` (the model, from the API only) and `check(og)` (its invariants); run one
with `python examples/golden/g04_cylinder_ogrid.py`. `tests/python/test_golden.py` runs them
all and replays their journals.
