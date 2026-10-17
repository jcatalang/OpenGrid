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
| Boundary conditions | `patch` (typed group of boundary connectors / faces), `zone` | `create_patch`, `set_boundary_type`, `create_zone`, `pair_periodic` |

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

- `og.delete(names)` refuses while other entities depend on them (the connectors on a
  curve, the domains on a connector, the blocks on a domain), naming them;
  `og.delete(names, dependents=True)` takes them along (`og.deletion_dependents(names)`
  lists them first; the GUI asks). Patches, zones and smoothing groups never block: they
  lose the members. A sketch takes its outputs along; one output on its own is refused
  (`og.remove_sketch_blocks("s", [1])` takes a block out of a sketch, which is what deleting
  a sketch's domain does in the GUI).
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
    loop = og.create_freehand(points, close_tolerance=0.01)  # ends near the start: closed
    d = og.unstructured_domain_from_curves(loop, 41)  # connectors + triangles, one step
    domains = og.domains_from_curves(curves)          # crossing curves -> structured domains
    og.domains_from_curves(curves, unstructured=True)  # regions without 4 corners: triangles
    og.move_handle("h1", "points", 6, (1.0, 0.9, 0))  # as a mouse drag: crossings follow
    og.insert_handle("h1", (1.4, 0.9, 0))            # a point (control point) added there
    og.remove_handle("h1", 6)                        # ... and taken away; one undo step each

### Extruded volumes (a 2D unstructured mesh, one or more cells thick)
    v = og.create_extruded_volume("disc", "translate", vector=(0, 0, 0.1), dimension=2)
    og.create_patch("front", "symmetry", [{"volume": v, "face": "top"}])
    og.create_patch("back", "symmetry", ["disc"])        # the domain is the bottom
    og.create_patch("wall", "wall", ["cc"])              # a connector: its side wall
    og.export_cgns("slab.cgns")                          # prisms / hexahedra, 3 boundaries

Modes: "translate" (``vector``), "rotate" (axis ``vector`` through ``axis_point``, ``angle``
in radians) and "normal" (``distance``, planar domains); ``distribution`` as for connectors.
A structured domain extrudes to a block (``create_extruded_block``).

### Changing a connector's points once domains stand on it
    og.set_dimension("c-bottom", 41)   # the connectors opposite to it in structured domains follow

A structured domain needs equal counts on opposite edges, so `set_dimension` spreads the
change through every structured domain using the connector (and so across the faces of
blocks); of several connectors on an opposite edge the longest takes it (`og.last_warnings`
says which). A sketch's connector sets its edge family in the sketch, and
`og.set_dimension("s-c0-1", 41, distribution={"type": "tanh", "begin_spacing": 1e-3})` that
edge's distribution there too (the sketch owns the rest of its parameters). `og.modify("c",
dimension=41)` changes that connector alone and is refused when a domain would become invalid.
Likewise `og.set_unstructured_generator("s-d1", smoothing="laplacian", inflation={...})` on a
sketch's block sends the settings to the sketch: the triangulation settings to its `generator`
(shared by its unstructured blocks), the inflation to walls on the block's edges.

### Boundary types in one step, and changing a domain's kind
    og.set_boundary_type(["c-bottom", "c-top"], "wall")   # they join the "wall" patch (made if absent)
    og.set_boundary_type(["c-left"], "inlet", name="in")   # a named patch of that type
    og.convert_domain("d")   # structured <-> unstructured on the same connectors (name kept)

A member belongs to one patch: `set_boundary_type` takes it out of any other, and
`clear_boundary_type(members)` takes members out of every patch. A connector shared by two
domains is an interface, found automatically: putting it in a patch is refused, saying so.
Conversion needs a single loop of four edges for the structured direction and is refused
for sketch domains (set the sketch's unstructured blocks), group members and faces of
blocks or volumes.

    og.make_periodic_pair(["c-left"], ["c-right"])        # patches periodic-1a / periodic-1b, paired (GUI: Pair as periodic)
    og.make_periodic_pair(["c-bottom"], ["c-top"], name="span")   # span-a / span-b
    og.extrude_model(0.05)                                 # the 2D model one cell thick along z:
    # blocks and volumes "<domain>-ext", the 2D patches on their sides (pairings kept), zones
    # extended, and the patches "front" / "back" as a periodic pair (ends="symmetry": planes)

`extrude_model` is what a 2D CFD case needs for a solver that reads 3D meshes: one cell
across, the ends periodic. The 2D model stays underneath (the domains are the fronts); undo
brings it back.

### Periodic pairs in 3D
    og.create_patch("back", "periodic", ["disc"])
    og.create_patch("front", "periodic", [{"volume": v, "face": "top"}])
    og.pair_periodic("back", "front")                    # the translation: centroid to centroid
    og.create_patch("side0", "periodic", [{"block": b, "face": "kmin"}])
    og.create_patch("side1", "periodic", [{"block": b, "face": "kmax"}])
    og.pair_periodic("side0", "side1", translation=(0, 0, 0),   # a sector turned about y
                     rotation_angle=math.pi / 3, rotation_axis=(0, 1, 0))

Any boundary faces pair (block faces, volume faces, the ends and side walls of extruded
volumes) when the source faces map node by node and face by face onto the target's under
x_target = R(axis, angle)(x_source - centre) + centre + translation; the boundary report says
how many nodes have no image and where the first one lands. The exports carry the pairs: CGNS
periodic connectivity (1-to-1 between block faces that map corner to corner, point lists
otherwise), Gmsh ``$Periodic``, VTKHDF field data, and the high-order meshes extend them to the
edge and face nodes. In the GUI the property panel of a periodic patch has the partner and
the transform (an empty translation is found from the centroids).

A double-click on a polyline, spline or NURBS curve (or on a connector lying on it) adds a
point there; a double-click on a handle removes it. NURBS control points are added by knot
insertion (the curve is unchanged) and removed by knot removal (exact when the curve allows
it, otherwise the log says how far the curve moved; the end control points stay).

In the GUI a freehand stroke released within the pick tolerance of its start is a closed
loop, the new curves are selected, and **Unstructured domain** meshes selected curves
directly (one connector per curve with the "points" box's points, one undo step). Select an
outer curve (or several forming a loop) together with the inner loops and the domain is the
region between them: the largest loop is the outer boundary, the others are holes (a body in
a far field). A stroke that crosses its own beginning encloses a region for "Domains from
curves"; the tails beyond the crossing are dropped. Dragging a closed curve's seam handle
moves both ends together, so the loop stays closed.

### Sketch topology
```python
og.import_cad("airfoil.step", name="cad", domains="none")
og.delete(og.names("connector"))
og.create_sketch(vertices, blocks, dimensions={(0, 1): 41, (1, 2): 31}, relax=False, name="c")
og.assign_sketch_geometry("c", og.names("nurbs"))  # curves near the sketch's edges
og.edit_sketch("c", walls={(1, 2): 1e-3, (0, 1): 2e-3}, relax=True)  # (0, 1): a wake cut
print(og.sketch_report("c"))  # where each vertex went
og.smooth_group("c-group")
og.move_handle("c", "vertices", 9, (0.8, 0.65))  # as a mouse drag
```
A sketch edge can adopt a connector made by hand (`connectors={(2, 3): "c-base"}`): the
block then uses that connector itself, so an unstructured (or any other) domain on its other
side shares its nodes, the edge family takes its point count and the group keeps it fixed.
The GUI's Sketch block tool adopts the connectors whose ends its clicks snapped to, "Assign
selected geometry to sketch" adopts selected connectors under the edges, and the log suggests
adoptions for edges drawn along a connector. Editing the adopted connector regenerates the
sketch. (The sketch's own connectors cannot be adopted by another sketch.)

Blocks can be meshed unstructured (`unstructured=[1, 3]`, or the "Unstructured sketch
blocks" toggle next to the Sketch block tool): those blocks become unstructured domains on
the same connectors, with the sketch's `generator` (triangles, quad-dominant, sizes) and the
walls on their edges as inflation layers (`walls={(2, 5): {"spacing": 1e-3, "layers": 20,
"growth": 1.2}}`). They stay out of the smoothing group; the connectors they share with
structured blocks are pinned so both sides keep the same nodes; with `extrude` they give
volumes of prisms and hexahedra beside the structured blocks.

Wall layers keep grid lines normal to a wall through the boundary layer (11a): a wall entry
`{"spacing": 1e-5, "layers": 30, "growth": 1.2}` (sketches: `walls={(1, 2): {...}}`; groups:
`{"connector": ..., ...}`; templates: `layers=30, growth=1.2` for all walls); the first 30 layers
are marched along the wall normals and kept, the smoothing does the rest. Leave about 20 grid
lines beyond the layers.
A 2D sketch with `extrude={"vector": (0, 0, 2), "points": 13}` sweeps its smoothed grid into
blocks; a 3D sketch takes `symmetry=[(a, b, c, d), ...]` for planar boundary faces that are
symmetry planes (grid lines cross them straight; their nodes move freely in the plane).
The sketch owns `c-c0-1`, `c-d0`, ..., `c-group`; edit the sketch, not them. See
`examples/10d_sketch_topology.py` (2D airfoil, sphere in a box, wing stub).

### Templates and automatic topology
```python
og.create_template("cgrid", body=["upper", "lower"], wall=1e-3)  # 2D C-grid
og.create_template(
    "bump_channel",
    bottom="bottom",
    right="right",
    top="top",
    left="left",
    blocks=(3, 1),
    wall=5e-3,
    length=1.0,
)  # 3D: the smoothed 2D grid, extruded
og.create_template(
    "wing_in_box", body=["upper", "lower"], wall=1e-5, length=2.0, taper=0.6, sweep=20.0
)  # 3D blocks (smoothed in 3D, end planes symmetry planes), wing on ruled surfaces
layout = og.cross_field(curves)  # field, singularities, separatrices
print(layout.field.index_sum(), layout.field.euler())
og.auto_topology(curves)  # a sketch of four-sided blocks
```
See `examples/10e_templates_auto.py`.

## High-order elements
The elements of every mesh cell, with the boundary nodes on the geometry (curved walls) and
the elements next to them following:
```python
mesh = og.high_order(order=3)  # triangles, quads, tets, pyramids, prisms, hexes + boundary
print(mesh.curving)  # moved nodes, scaled Jacobian of the curved elements (straight-sided),
# the repair of elements that folded or lost quality, and the rigorous validity check
flat = og.high_order(order=3, curved=False)  # straight-sided elements
print(mesh.counts(), len(mesh.points))  # shared edge / face nodes merged exactly
conn = mesh.connectivity("tetra", "cgns")  # (elements, 20) in CGNS's node order
og.export_high_order("case_p3.msh", order=3)  # Gmsh 4.1 high-order types, physical groups
og.export_high_order("case_p3.cgns", order=3)  # SIDS TETRA_20 ... HEXA_64, BCs per patch
og.export_high_order("case_p3.vtkhdf", order=3)  # VTK Lagrange cells (ParaView)
og.export_high_order("case_p3_gll.h5", order=3, nodes="gll")  # quads / hexes: GLL nodes
```
Orders 1-6 (CGNS: up to 4); exports take ``curved=False`` too. Folded elements after curving
raise a ValidationError naming them. The curved elements of 11e agglomerate p x p structured
grid cells into one element (interval counts multiples of the order):
```python
og.create_template("cgrid", body=["upper", "lower"], wall=1e-5, layers=24, interval_multiple=4)
mesh = og.curved_high_order(order=4, nodes="gll")  # (elements, 25, 3) GLL nodes, curved
print(mesh.scaled_jacobian.min())
og.export_curved_high_order("naca_p4.vtkhdf", order=4)  # VTK Lagrange cells
og.export_curved_high_order("naca_p4_gauss.h5", order=4, nodes="gauss")  # DG solution points
```
See `examples/11g_curved_elements.py`, `examples/11f_high_order_elements.py` and
`examples/11e_high_order.py`.

## GPU (experimental)
```python
opengrid.gpu_info()  # {"available": True, "device": "Apple M1 Max", ...} after a self-test
opengrid.set_gpu("on")  # 3D group smoothing: line sweeps on the GPU; "auto": >= 200k nodes
```
Off by default: on the 1M-cell ball it is 1.2x faster per iteration and stops at a relative
residual of 1e-5 (the CPU at 1e-6); see docs/PROGRESS.md (11d).

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
