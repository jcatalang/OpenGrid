# OpenGrid scripting API reference

Generated from the docstrings by `tools/gen_api_docs.py`; do not edit by hand. For concepts and
worked examples see `docs/SCRIPTING.md`; the golden cases in `examples/golden/` are complete
scripts.

```python
import opengrid

og = opengrid.Session()  # every method below is called on a Session
```

Module functions: `opengrid.replay(source)` runs a journal script and returns its session;
`opengrid.ValidationError` (with `.issues`, each an `opengrid.Issue(where, message)`) is raised
for invalid input; nothing is changed when it is raised.

## Geometry (database entities)

### `create_point(xyz, name: str | None = None) -> str`

Database point.

### `create_line(p0, p1, name: str | None = None) -> str`

Straight line from p0 to p1 (2 or 3 coordinates each).

### `create_arc(center, radius: float, start_angle: float, end_angle: float, normal=(0.0, 0.0, 1.0), x_axis=(1.0, 0.0, 0.0), name: str | None = None) -> str`

Arc of a circle; angles in radians measured from ``x_axis`` about ``normal``.

### `create_arc_3pt(start, through, end, name: str | None = None) -> str`

Arc from ``start`` through ``through`` to ``end``.

### `create_circle(center, radius: float, normal=(0.0, 0.0, 1.0), name: str | None = None) -> str`

Full circle starting on the x axis (or y axis if the normal is along x).

### `create_polyline(points, name: str | None = None) -> str`

Polyline through the points (straight segments).

### `create_spline(points, method: str = 'natural', start_tangent=None, end_tangent=None, name: str | None = None) -> str`

Interpolating cubic spline; method is 'natural', 'clamped' or 'catmull_rom'.

### `create_nurbs(degree: int, knots, control_points, weights=None, name: str | None = None) -> str`

NURBS curve: full knot vector (len(control_points) + degree + 1), optional weights.

### `create_plane(origin, u_vector, v_vector, name: str | None = None) -> str`

Plane patch x = origin + u u_vector + v v_vector, (u, v) in [0, 1]^2.

### `create_ruled_surface(curve0, curve1, name: str | None = None) -> str`

x = (1 - v) curve0(u) + v curve1(u).

### `create_extruded_surface(curve, direction, name: str | None = None) -> str`

x = curve(u) + v direction.

### `create_revolved_surface(curve, end_angle: float, axis_point=(0.0, 0.0, 0.0), axis=(0.0, 0.0, 1.0), start_angle: float = 0.0, name: str | None = None) -> str`

``curve`` rotated about the axis from start_angle to end_angle (radians).

### `create_nurbs_surface(degree_u: int, degree_v: int, knots_u, knots_v, control_points, weights=None, name: str | None = None) -> str`

Tensor-product NURBS; ``control_points`` is a list of rows (one per u index).

### `geometry(name: str)`

The C++ curve of an entity (None for points).

### `length(name: str) -> float`

Arc length of a curve (0 for points).

### `project_point(name: str, point) -> tuple[float, tuple, float]`

Closest point on a curve: (t, point, distance).

### `project_to_surface(name: str, point) -> tuple[float, float, tuple, float]`

Closest point on a surface: (u, v, point, distance).

## Connectors

### `create_connector(segments, dimension: int, distribution='uniform', mode: str = 'arc_length', breakpoints=(), propagate_spacing: bool = True, name: str | None = None) -> str`

Connector on a curve or curve chain.

``segments``: a curve name, or a list whose items are curve names or
``{"curve": name, "t0": .., "t1": ..}`` / ``(name, t0, t1)`` pieces (t1 < t0 runs a curve
backwards). ``distribution``: a type name ('uniform', 'geometric', 'tanh', 'sine') or a
dict with ``type`` and optional ``begin_spacing``, ``end_spacing`` (absolute lengths) and
``growth_ratio`` (geometric). ``breakpoints``: list of ``{"position": fraction,
"index": node, "distribution": ...}``; the breakpoint's distribution applies to the
sub-connector after it.

### `create_connectors_on_curves(curves, dimension: int, distribution='uniform') -> list[str]`

One connector per curve (a single undo step).

### `create_chain_connector(curves, dimension: int, distribution='uniform', name: str | None = None) -> str`

One connector along several curves, ordered and oriented automatically.

### `set_dimension(name: str, dimension: int) -> str`

Set a connector's number of points.

### `set_distribution(name: str, type: str, begin_spacing: float | None = None, end_spacing: float | None = None, growth_ratio: float | None = None) -> str`

Set a connector's distribution ("uniform", "geometric", "tanh", "sine").

### `connector_nodes(name: str) -> opengrid.model.connectors.ConnectorNodes`

Nodes and diagnostics of a connector (cached until the model changes).

### `connector_report(name: str) -> str`

One-paragraph summary: points, length, achieved spacings, stretching, warnings.

## Domains (2D, and faces in 3D)

### `create_structured_domain(connectors=None, edges=None, generator=None, name: str | None = None, surface=None) -> str`

Structured 2D domain.

Either ``edges``: four lists (counter-clockwise loop: j = 0 with i increasing, i = imax,
j = jmax with i decreasing, i = 0) of connector names or ``(name, reversed)`` /
``{"connector": name, "reversed": bool}`` items; or ``connectors``: a set of connector
names assembled automatically (one closed loop, or two closed loops joined by a cut).
``generator``: ``{"method": "tfi" | "hermite_tfi", "hermite_walls": (0,),
"hermite_strength": 1.0, "elliptic": None | {...}}`` (see opengrid.model.domains).

### `create_domain(edges, generator=None, name: str | None = None, surface=None) -> str`

Structured domain from explicit edges (the form recorded in journals). With
``surface`` (a surface entity), the grid is the parametric TFI on that surface.

### `create_unstructured_domain(loops=None, generator=None, name: str | None = None, connectors=None, surface=None) -> str`

Unstructured (triangle) domain.

Either ``loops``: the outer loop (counter-clockwise) followed by holes (clockwise), each a
list of connector names or ``(name, reversed)`` items; or ``connectors``: a set of
connector names assembled into closed loops automatically (largest loop = outer).
``generator``: ``{"algorithm": "delaunay", "min_angle": 25, "growth_rate": 1.2,
"max_size": None, "sources": [], "max_points": 2000000, "refine": True}``; sources are
``{"p0", "p1" (optional; line source), "size", "radius"}``. ``algorithm`` is
``"delaunay"`` (constrained Delaunay + refinement) or ``"advancing_front"`` (near
equilateral triangles that follow the size field closely; ``refine`` and ``min_angle``
apply to the Delaunay mesher only). Further keys: ``element_type`` ("triangles",
"quad_dominant", "quads"), ``smoothing`` ("none", "laplacian", "length_weighted",
"angle_based"), ``smoothing_iterations``, ``edge_swap`` and ``inflation``: boundary
layers from wall connectors, e.g. ``{"connectors": ["c-wall"], "first_height": 1e-5,
"growth_rate": 1.2, "layers": 30, "total_thickness": None, "transition": True}``.

### `set_domain_generator(name: str, **changes) -> str`

Update generator settings, e.g. ``set_domain_generator("d", method="hermite_tfi")`` or
``elliptic={"walls": [{"edge": 0, "spacing": 1e-5}]}`` (merged into the current elliptic
settings; ``elliptic=None`` switches smoothing off). Elliptic keys: ``control``
("laplace", "thomas_middlecoff", "ttm"), ``walls``, ``attractions``, ``solver``
("line_sor", "point_sor"), ``omega``, ``max_iterations``, ``tolerance`` and
``acceleration`` (Anderson depth, default 5; 0: plain relaxation).

### `set_unstructured_generator(name: str, **changes) -> str`

Change some generator settings of an unstructured domain.

### `rotate_domain(name: str, steps: int = 1) -> str`

Cycle the edges (edge ``steps`` becomes edge 0); orientation is preserved.

### `smooth_domain(name: str, progress=None) -> opengrid.model.domains.SmoothingResult`

Run the domain's elliptic smoother (cached); ``progress`` is an optional
``opengrid._core.SmoothingProgress`` for monitoring and cancelling from another thread.

### `domain_dimensions(name: str) -> tuple[int, int]`

(ni, nj) of a structured domain.

### `domain_grid(name: str, compute: bool = True) -> numpy.ndarray`

Grid points (nj, ni, 3); runs the elliptic smoother if configured and not cached
(``compute=False`` returns the TFI grid instead while smoothing is pending). A face of
a block in a smoothed 3D group is returned where the group moved it.

### `unstructured_mesh(name: str) -> dict`

``{"points" (n, 3), "triangles" (m, 3), "rejected", "warnings"}`` of an unstructured
domain (boundary points first, in loop order).

### `domain_quality(name: str) -> dict`

Per-cell quality metrics: arrays of shape (nj - 1, ni - 1) for structured domains,
(1, m) for unstructured ones (triangles first, then quads).

### `domain_status(name: str) -> str`

One-line summary of a structured or unstructured domain.

## Multi-block smoothing groups

### `create_smoothing_group(domains=(), walls=(), slide='auto', name: str | None = None, omega: float = 1.0, max_iterations: int = 2000, tolerance: float = 1e-06, acceleration: int = 20, blocks=()) -> str`

Smooth planar structured ``domains`` (or six-face ``blocks``, 3D) together as one
grid (as one grid): interfaces and block corners move, grid lines continue across
interfaces, and irregular points (3 or 5+ domains around a node) are relaxed too.
``walls``: connectors (or dicts ``{"connector", "spacing" (None: the initial
first-cell height), "decay", "relaxation", "max_source"}``) with first-cell height and
orthogonality control; they stay fixed. ``slide``: boundary connectors whose nodes
slide along their curves ("auto": every boundary connector that is not a wall; () for
none). Connectors used by domains outside the group stay fixed. Members' own elliptic
settings are ignored. The control is Laplace plus the wall sources; ``acceleration``
is the Anderson depth (default 20; 0: plain relaxation). With ``blocks``: ``walls``
and ``slide`` name block faces (domains); sliding faces move on their surface (or
plane), boundary edges between faces on different surfaces along their connectors;
faces shared with blocks or volumes outside the group stay fixed. Members' own
smoothing settings are ignored.

### `smooth_group(name: str, progress=None) -> opengrid.model.groups.GroupResult`

Run (or reuse) a smoothing group; ``progress``: optional
``opengrid._core.SmoothingProgress``. Returns grids per member domain id, residuals,
node counts (regular, irregular, sliding, fixed) and the connector roles.

### `group_status(name: str) -> str`

One-line summary of a smoothing group.

### `display_connector_nodes(name: str) -> opengrid.model.connectors.ConnectorNodes`

Connector nodes as drawn: where a smoothed group moved them (interfaces, sliding
boundaries), else the connector's own nodes.

## Blocks (structured 3D)

### `create_block(faces=None, extrusion=None, name: str | None = None, smoothing=None, quality_floor=None) -> str`

Structured block from six structured domains (``faces``, any order and orientation)
or from an ``extrusion`` of a domain: ``{"domain", "mode": "translate" | "rotate" |
"normal", "vector" (translation or rotation axis), "axis_point", "angle" (radians),
"distance" (normal), "dimension" (points along the extrusion), "distribution":
{"type", "begin_spacing", "end_spacing", "growth_ratio"}, "edges": (imin, imax, jmin,
jmax) conditions of hyperbolic extrusion ("free", "constant_x", "constant_y",
"constant_z", "periodic"), "dissipation"}``. ``smoothing``: 3D elliptic smoothing
``{"control": "laplace" | "thomas_middlecoff", "walls": [{"face": "kmin", "spacing",
"decay", "relaxation", "max_source"}], "omega", "max_iterations", "tolerance",
"acceleration" (Anderson depth, default 5; 0: plain relaxation)}``.
``quality_floor``: ``{"min_jacobian": 0.05, "max_skewness": 0.95}`` (None: the default);
cells missing it are reported by the boundary report, and exports refuse the model.

### `create_extruded_block(domain, mode: str = 'translate', name: str | None = None, smoothing=None, quality_floor=None, **settings) -> str`

Shorthand for create_block(extrusion={"domain": domain, "mode": mode, **settings}).

### `set_block_smoothing(name: str, smoothing) -> str`

Set (or remove with None) the 3D elliptic smoothing of a block.

### `block_grid(name: str) -> numpy.ndarray`

Block grid, shape (nk, nj, ni, 3).

### `block_quality(name: str) -> dict`

Per-cell quality metrics of a block, arrays of shape (nk - 1, nj - 1, ni - 1).

### `block_floor(name: str) -> dict`

Quality-floor check of a block: {"floor", "cells", "below_floor", "worst_jacobian",
"worst_skewness", "failing" (flat cell indices, i fastest)}.

### `block_status(name: str) -> str`

One-line summary of a block: points, cells, source, smoothing and quality floor.

## Volumes (unstructured / hybrid 3D)

### `create_volume(faces, generator=None, name: str | None = None) -> str`

Unstructured 3D volume (tetrahedra, pyramids on quad faces) enclosed by the given
domains (one or more closed shells). ``generator``: ``{"backend": "netgen",
"max_size": None, "grading": 0.3, "pyramids": "auto" | "always" | "never",
"inflation": None, "optimize": True}``. ``optimize``: improve cells with a scaled
Jacobian below 0.3 by moving interior nodes and pyramid apexes and swapping tetrahedra
(boundary and layer nodes never move). ``inflation``: boundary layers (prisms on
triangles, hexahedra on quads, fans at convex edges turning more than 60 degrees) from
wall faces, ``{"faces": [domains], "first_height", "growth_rate": 1.2,
"layers": 30, "total_thickness": None, "transition": True, "slide": True}``; with
``slide`` the layers slide along planar unstructured faces that meet the walls and are
used by this volume only (those faces are re-meshed around the layers).

### `set_volume_generator(name: str, **changes) -> str`

Change some generator settings of a volume.

### `volume_mesh(name: str) -> dict`

``{"points", "tets" (m, 4), "pyramids" (k, 5), "prisms" (k, 6), "hexes" (k, 8),
"faces" {domain id: boundary cells}, "inflation" (None or {"max_layers", "layers",
"height", "stops", ...}), "optimization" (None or {"moved", "swaps_23", "swaps_32",
"worst_before", "worst_after", "below_before", "below_after"}), "warnings", ...}``
(boundary nodes first, identical to the face meshes).

### `volume_quality(name: str) -> dict`

Per-cell quality arrays of a volume mesh (scaled Jacobian, angles, aspect ratio, ...).

### `volume_status(name: str) -> str`

One-line summary of a volume mesh: cell counts, points, shells, layers.

## Topology editing

### `split_connector(name: str, index: int, new_name: str | None = None) -> list[str]`

Split a connector at node ``index``: it keeps nodes 0..index, a new connector takes
the rest (explicit distributions keep every node). Domains, loops and patches using it
get both. Returns [name, new connector].

### `join_connectors(first: str, second: str) -> str`

Join two connectors sharing an end into ``first`` (all nodes kept); users must take
them consecutively. Returns the joined connector's name.

### `split_domain(name: str, i: int | None = None, j: int | None = None) -> list[str]`

Split a structured domain along grid line i (or j). The cut becomes a connector on a
polyline through the grid line; edge connectors crossing it are split. Returns [name
(part with the lower index), new domain, cut connector].

### `join_domains(first: str, second: str) -> str`

Join two structured domains sharing a whole edge into ``first``; connectors only
the shared edge used are deleted. Returns the joined domain's name.

### `split_block(name: str, axis: str, index: int) -> list[str]`

Split a six-face block at ``index`` along ``axis`` ("i", "j" or "k"): its side faces
are split and a new interior face is created. Returns [name, new block, new face].

### `join_blocks(first: str, second: str) -> str`

Join two six-face blocks sharing a face into ``first`` (side faces joined pairwise,
the shared face deleted). Returns the joined block's name.

### `ogrid_domain(name: str, fraction: float = 0.5, radial_points: int | None = None, radial_spacing: float | None = None, core: bool = True) -> list[str]`

O-grid (butterfly) a structured domain: four domains around a core whose corners are
centre + fraction * (corner - centre); radial connectors with ``radial_points`` (default
(ni + nj) / 4) and optional tanh ``radial_spacing`` at the outer boundary. ``core=False``
leaves the core as a hole. The domain is replaced; returns the new domains.

### `ogrid_block(name: str, fraction: float = 0.5, radial_points: int | None = None, radial_spacing: float | None = None, core: bool = True) -> list[str]`

O-grid (butterfly) a six-face block: six shell blocks around a core block (or only the
shell with ``core=False``), all hexahedral. Returns the new blocks.

## Patches, zones and topology

### `create_patch(name: str | None = None, patch_type: str = 'wall', members=(), periodic=None) -> str`

Named boundary-condition group of connectors; ``patch_type`` in PATCH_TYPES
(opengrid.model.patches). ``periodic`` is normally set with :meth:`pair_periodic`.

### `set_patch_type(name: str, patch_type: str) -> str`

Change a patch's boundary-condition type (wall, inlet, outlet, ...).

### `add_to_patch(name: str, connectors) -> str`

Add members (connectors in 2D, domains or block faces in 3D) to a patch.

### `remove_from_patch(name: str, connectors) -> str`

Remove members from a patch.

### `pair_periodic(source: str, target: str, translation=None, rotation_angle: float = 0.0, rotation_center=(0.0, 0.0, 0.0)) -> str`

Pair two periodic patches: target = R(angle)(source - center) + center + translation.

Without a translation and with no rotation, the translation that maps the centroid of
the source nodes onto that of the target is computed and stored explicitly.

### `unpair_periodic(source: str) -> str`

Remove the periodic pairing of a patch.

### `create_zone(name: str | None = None, members=()) -> str`

Named group of domains (volume tag).

### `add_to_zone(name: str, domains) -> str`

Add domains to a zone.

### `remove_from_zone(name: str, domains) -> str`

Remove members from a zone.

### `is_3d() -> bool`

True when the model has blocks (exports and boundary checks are then 3D).

### `topology()`

2D: connector uses of domains. 3D models: block faces (boundary / interfaces).

### `boundary_report()`

Boundary connectors (2D) or block faces (3D), interfaces (detected automatically),
boundary entities without a patch and all problems
(opengrid.model.topology.BoundaryReport).

### `unassigned_boundary_connectors() -> list[str]`

Boundary entities that belong to no patch.

## Import

### `import_airfoil_dat(path, name: str | None = None, mode: str = 'spline', method: str = 'natural') -> list[str]`

Import a Selig or Lednicer airfoil file.

mode 'spline': one spline TE -> upper -> LE -> lower -> TE (smooth leading edge);
mode 'split': '<name>-upper' and '<name>-lower' splines, both running LE -> TE, clamped
at the leading edge to the tangent of the single-loop spline (G1 at the LE);
mode 'polyline': the raw points joined by straight segments.

### `import_csv(path, name: str | None = None, mode: str = 'spline', method: str = 'natural') -> list[str]`

Import a point list as one spline (mode 'spline') or polyline (mode 'polyline').

### `import_cad(path, name: str | None = None, spacing: float | None = None, max_turn_angle: float = 15.0, domains: str = 'unstructured', domain_generator: dict | None = None, tolerance: float = 1e-06, scale: float = 1.0, heal: bool = True, sew: bool = True, split_closed: bool = True) -> opengrid.io.cad.CadImport`

Import a STEP (.step / .stp) or IGES (.iges / .igs) file (needs the OpenCASCADE
module, see opengrid.io.cad).

Creates one ``nurbs`` curve and one connector per CAD edge, and one domain per face
(``domains``: "unstructured", "auto" = structured for four-sided faces, or "none"),
on a ``nurbs_surface`` for curved faces. Connector dimensions follow ``spacing``
(default: 1/40 of the bounding-box diagonal) and ``max_turn_angle`` (degrees of
tangent turning per segment). ``tolerance`` is the healing / sewing precision (the
file's own tolerances are used when larger) and ``scale`` multiplies the coordinates
(OCCT reads lengths in millimetres). ``heal``, ``sew`` and ``split_closed`` (split
periodic faces and closed edges) switch the healing steps. ``domain_generator`` is
passed to every unstructured domain. Faces that cannot be meshed are skipped and
reported in the result (``CadImport.report()``); one undo step undoes the import.

## Export

### `export_vtkhdf(path, domains=None) -> dict`

Write domains to VTKHDF.

Without patches and zones: one UnstructuredGrid (quads and triangles, shared points
merged) with ``block_id``, ``zone_id``, cell indices and all quality metrics as cell
data. With patches or zones (the boundary must then validate): a
PartitionedDataSetCollection with one cell block per zone and one line block per patch,
grouped in the assembly as Volume/<zone> and Patches/<patch>.

### `export_cgns(path, base_name: str = 'OpenGrid') -> dict`

CGNS: structured zones (domains / blocks) and unstructured zones (unstructured domains /
volumes) in one base, families, BCs per patch and connectivity; the boundary must
validate.

### `export_gmsh(path) -> dict`

Gmsh MSH 4.1 with physical groups (zones: 2D, patches: 1D); the boundary must
validate.

### `export_plot3d(path, domains=None, binary: bool = True, three_d: bool = False) -> None`

Write domains as a multi-grid Plot3D file (2D, or 3D with nk = 1); in a model with
blocks, the blocks (3D).

## Editing, undo and the model

### `modify(name: str, **params) -> str`

Change some parameters of an entity, e.g. ``og.modify("line-1", p1=(2, 0))``.

### `rename(name: str, new_name: str) -> str`

Rename an entity (references follow; journals stay valid).

### `delete(names) -> int`

Delete entities (a name or a list); refused while other entities use them.

### `undo() -> None`

Undo the last command.

### `redo() -> None`

Redo the last undone command.

### `entity(name: str) -> opengrid.model.store.Entity`

The stored entity (id, name, kind, canonical params).

### `names(kind: str | None = None) -> list[str]`

Entity names in creation order, optionally of one kind ("line", "connector", ...).

## Projects and journals

### `new() -> None`

Start an empty model (clears the undo history).

### `save(path) -> None`

Write the project file (.ogp, canonical JSON).

### `load(path) -> None`

Replace the model by a project file. Clears undo history; journal restarts from it.

### `journal_text() -> str`

The session journal: a Python script that rebuilds the model (opengrid.replay).

### `save_journal(path) -> None`

Write the journal script to a file.
