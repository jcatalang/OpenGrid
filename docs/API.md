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

### `create_rectangle(center=(0.0, 0.0), width: float = 1.0, height: float = 1.0, name=None) -> list[str]`

Four lines counter-clockwise from the lower left: ``<name>-bottom``, ``-right``,
``-top``, ``-left`` (``name`` "rectangle-<n>" by default). One undo step; returns
the four names.

### `create_ellipse(center=(0.0, 0.0), a: float = 1.0, b: float = 0.5, name=None) -> str`

Exact ellipse (a rational NURBS) with semi-axes ``a`` along x and ``b`` along y,
starting on the +x side and running counter-clockwise.

### `create_airfoil(code: str = '0012', chord: float = 1.0, leading_edge=(0.0, 0.0), angle: float = 0.0, points: int = 101, name=None) -> list[str]`

A NACA 4-digit airfoil (``code`` "0012", "2412", ...) as two clamped splines from
the leading edge to the sharp trailing edge, ``<name>-upper`` and ``<name>-lower``
("airfoil-<n>" by default), of ``chord`` at ``leading_edge`` and ``angle`` of attack
in degrees (nose up). The two curves close on each other: select both for an
unstructured domain, or build C / O grids on them. One undo step.

### `create_box(center=(0.0, 0.0, 0.0), size=(1.0, 1.0, 1.0), name=None) -> list[str]`

The six faces of an axis-aligned box as plane patches ``<name>-xmin`` ...
``-zmax`` ("box-<n>" by default), normals outward. One undo step.

### `create_cylinder(center=(0.0, 0.0, 0.0), radius: float = 0.5, height: float = 1.0, axis=(0.0, 0.0, 1.0), caps: bool = True, name=None) -> list[str]`

A cylinder of revolution: ``center`` is the middle of the bottom face, the axis
points up through the top. Surfaces ``<name>-side`` and, with ``caps``, ``-bottom``
and ``-top`` (discs, degenerate at their centre), from the generating lines
``<name>-generatrix``, ``-bottom-radius``, ``-top-radius``. One undo step.

### `create_sphere(center=(0.0, 0.0, 0.0), radius: float = 0.5, name=None) -> list[str]`

A sphere as the surface of revolution ``<name>`` of a meridian half circle
``<name>-meridian`` (from the +z pole to the -z pole) about the z axis through the
centre ("sphere-<n>" by default). One undo step.

### `create_wing(code: str = '0012', root_chord: float = 1.0, tip_chord: float | None = None, span: float = 3.0, sweep: float = 0.0, dihedral: float = 0.0, points: int = 101, caps: bool = True, name=None) -> list[str]`

A straight tapered wing: NACA 4-digit sections at the root (y = 0) and the tip
(y = ``span``, chord ``tip_chord``, leading edge swept back by ``sweep`` degrees and
raised by ``dihedral`` degrees), chord along +x and thickness along z, as four
clamped splines ``<name>-root-upper`` ... ``-tip-lower``, the ruled surfaces
``<name>-upper`` and ``<name>-lower`` between them and, with ``caps``, the flat
``<name>-root`` and ``<name>-tip`` sections closing the body (ruled between each
section's upper and lower curves, so they meet at the leading and trailing edges).
"wing-<n>" by default; one undo step.

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

### `set_dimension(name: str, dimension: int, **changes) -> list[str]`

A connector's number of points, with the connectors opposite to it in structured
domains (and so across the faces of blocks) following, so the domains stay valid
(opengrid.model.dimension); of several connectors on an opposite edge the longest
takes the change (a warning says so). A sketch's connector sets its edge family in
the sketch, and its ``distribution`` that edge's distribution there (the sketch owns
the rest of its parameters). ``changes``: other parameters of the connector
(distribution, breakpoints, ...) set in the same step. The property panel's
connector edits go through this; ``modify`` changes one connector only and is
refused when a domain would become invalid. One undo step; returns the connectors
changed.

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

### `unstructured_domain_from_curves(curves, dimension: int, connectors=(), generator=None, distribution='uniform', name: str | None = None) -> str`

Unstructured (triangle) domain bounded by curves: one connector of ``dimension``
points on each curve, assembled together with the existing ``connectors`` into closed
loops (the largest = outer), then the domain, as one undo step (the GUI's
"Unstructured domain" on selected curves, e.g. a freehand loop). The curve ends must
meet exactly (``create_freehand`` with ``close_tolerance`` closes a stroke). Returns
the domain's name; ``generator`` as in ``create_unstructured_domain``.

### `set_domain_generator(name: str, **changes) -> str`

Update generator settings, e.g. ``set_domain_generator("d", method="hermite_tfi")`` or
``elliptic={"walls": [{"edge": 0, "spacing": 1e-5}]}`` (merged into the current elliptic
settings; ``elliptic=None`` switches smoothing off). Elliptic keys: ``control``
("laplace", "thomas_middlecoff", "ttm"), ``walls``, ``attractions``, ``solver``
("line_sor", "point_sor"), ``omega``, ``max_iterations``, ``tolerance`` and
``acceleration`` (Anderson depth, default 5; 0: plain relaxation).

### `set_unstructured_generator(name: str, generator=None, **changes) -> str`

Change generator settings of an unstructured domain (``generator``, a dict, and
``changes`` merged into the current ones). For a sketch's block the settings go to
the sketch: the triangulation settings to its ``generator`` (shared by its
unstructured blocks) and ``inflation`` to walls on the block's edges (first height,
layers, growth), which is what the property panel does there. One undo step.

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

## Freehand curves, curve networks, live editing

### `create_freehand(points, tolerance: float = 0.005, corner_angle: float | None = 60.0, name: str | None = None, close_tolerance: float = 0.0) -> list[str]`

NURBS curves fitted to a freehand stroke (pointer samples): split at corners (turns
of ``corner_angle`` degrees or more; None: never) and fitted by cubic B-splines within
``tolerance`` times the stroke length (opengrid.model.fitting). A stroke that ends
within ``close_tolerance`` (absolute; 0: never) of its start is closed into a loop
whose curves meet exactly (the GUI passes its pick tolerance), ready for
``unstructured_domain_from_curves``. One undo step; the journal records the fitted
curves. ``name`` is used for a single piece, else ``name-1``, ``name-2``...

### `domains_from_curves(curves, spacing: float | None = None, corner_angle: float = 30.0, unstructured: bool = False) -> list[str]`

Structured domains from planar curves drawn across each other (a closed loop of
curves, or a network): the curves are cut where they cross or turn sharply, every
region with four corners becomes a structured domain with connectors on the curve
pieces (point counts from ``spacing``, default 1/20 of the median piece length;
opposite sides made equal so neighbouring domains conform). Other regions are
reported in ``last_warnings`` (or become unstructured domains with
``unstructured``). One undo step; returns the domains created.

### `move_handle(name: str, field: str, index: int, xyz) -> str`

Move one defining point of a point or curve (``field``: "xyz", "p0"/"p1",
"points" or "control_points"; ``index`` -1 for single points), as a drag in the GUI
does: crossings of the curve with other curves are solved again and the connectors
meeting there follow, so domains built from a curve network stay closed. One undo
step.

### `insert_handle(name: str, xyz) -> str`

Add a defining point to a polyline or spline (on the curve where it is closest to
``xyz``, between the neighbouring points) or a control point to a NURBS curve (a knot
inserted at the closest parameter: the curve does not change), as a double-click on
a curve does in the GUI. Connectors on the curve keep their places and crossings are
solved again. One undo step.

### `remove_handle(name: str, index: int) -> str`

Remove defining point ``index`` of a polyline or spline, or control point ``index``
of a NURBS curve (knot removal: exact when the curve allows it, otherwise
``last_warnings`` says how far the curve moved; the end control points cannot be
removed), as a double-click on a handle does in the GUI. One undo step.

## Multi-block smoothing groups

### `create_smoothing_group(domains=(), walls=(), slide='auto', name: str | None = None, omega: float = 1.0, max_iterations: int = 2000, tolerance: float = 1e-06, acceleration: int = 20, blocks=(), joints: str = 'fixed', pinned=(), slide_limit: float = 0.5, symmetry=()) -> str`

Smooth planar structured ``domains`` (or six-face ``blocks``, 3D) together as one grid
: interfaces and block corners move, grid lines continue across
interfaces, and irregular points (3 or 5+ domains around a node) are relaxed too.
``walls``: connectors (or dicts ``{"connector", "spacing" (None: the initial
first-cell height), "decay", "relaxation", "max_source", "layers", "growth"}``) with
first-cell height and orthogonality control; they stay fixed. ``layers`` > 0: the
first layers off the wall are marched along the front normals (heights spacing *
growth^k, default growth 1.2) and kept, so grid lines stay normal to the wall
through them; the control then acts on the first cell beyond. ``slide``: boundary
connectors whose nodes slide along their curves ("auto": every boundary connector
that is not a wall; () for none). Connectors used by domains outside the group stay
fixed. Members' own elliptic settings are ignored. The control is Laplace plus the
wall sources; ``acceleration`` is the Anderson depth (default 20; 0: plain
relaxation). With ``blocks``: ``walls`` and ``slide`` name block faces (domains);
sliding faces move on their surface (or plane), boundary edges between faces on
different surfaces along their connectors; faces shared with blocks or volumes
outside the group stay fixed. Members' own smoothing settings are ignored.
``joints="slide"``: block corners between two sliding connectors on one curve slide
along it (3D: face corners inside one sliding surface slide on it); ``pinned``:
connectors (3D: faces) kept fixed even between members, e.g. a wake cut (a wall
between members is an internal wall: fixed, with wall control on both sides, 2D);
``slide_limit``: how far a sliding node may move, in initial intervals; ``symmetry``
(3D): planar faces that are symmetry planes (nodes move as interior nodes mirrored
across the plane, without the slide limit).

### `smooth_group(name: str, progress=None) -> opengrid.model.groups.GroupResult`

Run (or reuse) a smoothing group; ``progress``: optional
``opengrid._core.SmoothingProgress``. Returns grids per member domain id, residuals,
node counts (regular, irregular, sliding, fixed) and the connector roles.

### `group_status(name: str) -> str`

One-line summary of a smoothing group.

### `display_connector_nodes(name: str) -> opengrid.model.connectors.ConnectorNodes`

Connector nodes as drawn: where a smoothed group moved them (interfaces, sliding
boundaries), else the connector's own nodes.

## Sketch topology (floating wireframe)

### `create_sketch(vertices, blocks, dimensions=(), spacing: float | None = None, curves=(), surfaces=(), walls=(), fixed=(), pinned=(), symmetry=(), distributions=(), relax: bool = True, tolerance: float = 1e-08, max_iterations: int = 4000, extrude=None, interval_multiple: int = 1, connectors=(), unstructured=(), generator=None, name: str | None = None) -> str`

Sketch topology (wireframe of blocks, opengrid.model.sketch): ``vertices``
(points drawn loosely around the geometry) and ``blocks`` (4 vertex indices each in
2D, counter-clockwise or not; 8 in 3D, VTK hexahedron order). ``dimensions``: points
per edge, ``{(a, b): n}`` on any edge of each edge family (or ``spacing`` for the
families without one); ``curves``: ``{(a, b): curve}`` assigns edges to database
curves, ``connectors``: ``{(a, b): connector}`` makes 2D edges adopt existing
connectors (the domain shares their nodes with whatever is on their other side, the
edge family takes their point count), ``surfaces``: ``{(a, b, c, d): surface}`` 3D
boundary faces to surfaces;
``walls``: ``{(a, b): first-cell height or None}`` (edges, 2D; faces, 3D) or dicts
with ``"edge"``/``"face"``, ``"spacing"``, ``"decay"``, ``"relaxation"``,
``"max_source"``, ``"layers"``, ``"growth"`` (wall layers, as in
``create_smoothing_group``; ``{(a, b): {"spacing": h, "layers": n}}`` works too);
``fixed``: vertices that must not float; ``distributions``:
``{(a, b): distribution}`` from a to b. Generates connectors, domains (3D: also
blocks) and a smoothing group named ``<sketch>-...`` that the sketch owns; with
``relax`` (2D) vertices joining two edges on one curve slide along it until the
grid stops moving them. ``extrude`` (2D): ``{"vector": (0, 0, L), "points": n,
"distribution": {...}}`` (spacings as lengths) also makes one extruded block per
domain: the smoothed 2D grid extruded, the same in every section;
``interval_multiple`` m: interval counts rounded up to multiples of m (for
``high_order`` meshes of order m). ``unstructured``: indices of the blocks meshed
as unstructured domains on the same connectors (``generator``: their triangulation
settings; a wall on their edges becomes inflation layers with the wall's first-cell
height, layers and growth); they stay out of the smoothing group, and the connectors
they share with structured blocks are pinned. One undo step; warnings (e.g.
unassigned edges lying on geometry) are in ``last_warnings``.

### `edit_sketch(name: str, **changes) -> str`

Change sketch parameters (any of those of create_sketch) and regenerate its
outputs in place (same names: patches and other references to them stay valid).
One undo step.

### `remove_sketch_blocks(name: str, blocks) -> str`

Take blocks (indices) out of a 2D sketch and regenerate it: the remaining blocks
are renumbered (domain names follow), edges left in no block lose their walls,
pinned state, dimensions and assignments; the vertices stay. The last block cannot
go (delete the sketch). The GUI deletes a sketch's domain this way. One undo step.

### `regenerate_sketch(name: str) -> str`

Generate a sketch's outputs again (e.g. after its geometry changed).

### `sketch_report(name: str) -> str`

Vertices (where they went and whether they float), edge families with their
dimensions, errors and suggested assignments of a sketch.

### `sketch_suggestions(name: str, tolerance: float | None = None) -> list[dict]`

Unassigned boundary edges (3D: faces) of a sketch lying near a curve (surface):
``{"what", "edge" | "face", "curve" | "surface", "message"}``. Nothing is assigned;
see accept_sketch_suggestions.

### `accept_sketch_suggestions(name: str, tolerance: float | None = None) -> str`

Assign every suggested curve and surface (sketch_suggestions) in one edit.

### `assign_sketch_geometry(name: str, geometry, tolerance: float | None = None) -> str`

Assign the given curves (and 3D: surfaces) to the unassigned boundary edges (faces)
of a sketch that lie near them - within ``tolerance``, by default the ends within a
fifth of the edge length and the middle within half of it (the nearest of them; the
GUI's "assign selected geometry"); given connectors are adopted by the 2D edges
whose vertices lie at their ends. One edit.

### `create_template(kind: str, name: str | None = None, **options) -> str`

A topology template (opengrid.model.templates): the far-field geometry it needs and
a sketch with no user-placed vertices. 2D: ``"ogrid"`` (body, radius, blocks,
radial_points, wall), ``"butterfly"`` (boundary, core, radial_points), ``"cgrid"``
(body with a sharp trailing edge, radius, downstream, ring, radial_points=(inner,
outer), wall, wake), ``"hgrid"`` (bottom, right, top, left, blocks=(nx, ny),
normal_points, wall, wall_sides), ``"body_in_channel"`` (body, bottom, right, top,
left, margin, radial_points, wall). 3D: ``"pipe"``, ``"cylinder_in_channel"``,
``"bump_channel"``, ``"wing_in_box"`` (the 2D template's options plus ``length`` and
``span_points``: extruded along z; ``wing_in_box`` also ``taper`` and ``sweep``) and
``"sphere_in_box"`` (center, radius, box, cube, points, radial_points, outer_points,
wall). All take ``spacing`` (default: a twentieth of the median boundary edge),
``relax``, ``layers`` / ``growth`` (wall layers on every wall; a C-grid's wake cut
then takes the wall spacing) and ``interval_multiple`` (interval counts rounded up to
multiples, for ``high_order``). One undo step; returns the sketch.

### `auto_topology(curves, size: float | None = None, spacing: float | None = None, relax: bool = False, name: str | None = None) -> str`

Automatic 2D block topology (opengrid.model.crossfield): a cross field of the region
bounded by the closed loops of ``curves`` (triangulated with mesh ``size``, default
1/200 of their total length), its separatrices traced from the singularities and
corners to the boundary (created as polylines ``separatrix-<k>``) and a sketch whose
blocks are the four-sided faces they cut (point counts from ``spacing``). One undo
step; returns the sketch (edit it like any other).

### `cross_field(curves, size: float | None = None)`

The cross field of a planar region and its separatrices, without creating anything:
a Layout with ``field`` (points, triangles, complex values z = exp(4 i theta),
corners, singularities, ``index_sum()`` and ``euler()``) and ``separatrices``.

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

### `create_volume(faces=None, generator=None, name: str | None = None, extrusion=None) -> str`

Unstructured 3D volume (tetrahedra, pyramids on quad faces) enclosed by the given
domains (one or more closed shells), or - with ``extrusion`` instead of ``faces`` -
prisms and hexahedra swept from an unstructured domain (see
``create_extruded_volume``). ``generator``: ``{"backend": "netgen",
"max_size": None, "grading": 0.3, "pyramids": "auto" | "always" | "never",
"inflation": None, "optimize": True}``. ``optimize``: improve cells with a scaled
Jacobian below 0.3 by moving interior nodes and pyramid apexes and swapping tetrahedra
(boundary and layer nodes never move). ``inflation``: boundary layers (prisms on
triangles, hexahedra on quads, fans at convex edges turning more than 60 degrees) from
wall faces, ``{"faces": [domains], "first_height", "growth_rate": 1.2,
"layers": 30, "total_thickness": None, "transition": True, "slide": True}``; with
``slide`` the layers slide along planar unstructured faces that meet the walls and are
used by this volume only (those faces are re-meshed around the layers).

### `create_extruded_volume(domain, mode: str = 'translate', name: str | None = None, generator=None, **settings) -> str`

Volume of prisms (from the triangles) and hexahedra (from the quads) swept from an
unstructured ``domain``: ``mode`` "translate" (``vector``), "rotate" (axis ``vector``
through ``axis_point``, ``angle`` in radians) or "normal" (``distance`` along the plane
normal of a planar domain), with ``dimension`` points across the layers and a
connector-style ``distribution`` (as ``create_extruded_block``; the base is mirrored
when the sweep runs against its normal). Its boundary for patches: the domain (the
bottom), ``{"volume": name, "face": "top"}`` and the domain's connectors (the side
walls), so a 2D model's patches carry over to the one-cell-thick 3D mesh some solvers
want. ``generator`` only sets the quality floor. Structured domains extrude to blocks.

### `set_volume_generator(name: str, **changes) -> str`

Change some generator settings of a volume.

### `volume_mesh(name: str) -> dict`

``{"points", "tets" (m, 4), "pyramids" (k, 5), "prisms" (k, 6), "hexes" (k, 8),
"faces" {domain id, ("top", domain id) or ("side", connector id): outward boundary
cells}, "extrusion" (None or {"mirrored", "layers"}), "inflation" (None or
{"max_layers", "layers", "height", "stops", ...}), "optimization" (None or {"moved",
"swaps_23", "swaps_32",
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

### `convert_domain(name: str) -> str`

Turn a structured domain into an unstructured one on the same connectors, or an
unstructured domain (one loop that assembles into four edges) into a structured one;
the name stays and zones and 3D patches follow. Refused for a sketch's domains (set
the sketch's unstructured blocks), members of a smoothing group and faces of blocks or
volumes. One undo step.

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

Add members (connectors in 2D, domains, block faces or extruded-volume ends
``{"volume", "face"}`` in 3D; connectors are also the side walls of extruded volumes)
to a patch.

### `remove_from_patch(name: str, connectors) -> str`

Remove members from a patch.

### `set_boundary_type(members, patch_type: str, name: str | None = None) -> str`

Give boundary members (connectors in 2D; domains, ``{"block", "face"}``,
``{"volume", "face"}`` and side-wall connectors in 3D) a boundary type: they join
the patch called ``name`` (created when absent), else the first patch of that type,
else a new patch named after the type - leaving the patches they were in (a member
belongs to one patch). The GUI's "Boundary type" menu. One undo step; returns the
patch name.

### `clear_boundary_type(members) -> list[str]`

Take boundary members out of every patch: they are unassigned again (the GUI's
"Clear type"). Returns the patches that changed; one undo step.

### `pair_periodic(source: str, target: str, translation=None, rotation_angle: float = 0.0, rotation_center=(0.0, 0.0, 0.0), rotation_axis=(0.0, 0.0, 1.0)) -> str`

Pair two periodic patches: target = R(axis, angle)(source - center) + center +
translation (the rotation by ``rotation_angle`` radians about ``rotation_axis``
through ``rotation_center``, right-hand rule).

Without a translation and with no rotation, the translation that maps the centroid of
the source nodes onto that of the target is computed and stored explicitly. 2D: the
source connectors must map node by node onto the target's; 3D (block faces, volume
faces, the ends and side walls of extruded volumes): the faces must map node by node
and face by face. The boundary report names what does not match.

### `make_periodic_pair(source, target, translation=None, rotation_angle: float = 0.0, rotation_center=(0.0, 0.0, 0.0), rotation_axis=(0.0, 0.0, 1.0), name: str | None = None) -> list[str]`

Two periodic patches from the two sides (connectors in 2D; domains, block faces or
extruded-volume ends in 3D) and their pairing, in one step (the GUI's "Pair as
periodic"): patches ``<name>-a`` and ``<name>-b`` (``periodic-1``, ``periodic-2``, ...
when no name is given). Without a
translation and no rotation, the translation maps the centroid of the source side
onto the target's. Returns [source patch, target patch]; one undo step.

### `extrude_model(thickness: float, layers: int = 1, ends: str = 'periodic') -> list[str]`

The one-cell-thick (``layers`` cells) 3D model of a 2D one: every domain extruded
by ``thickness`` along z (blocks from structured domains, volumes of prisms and
hexahedra from unstructured ones, named ``<domain>-ext``), the 2D patches carried
over to the side faces with their periodic pairings, zones extended, and the fronts
and backs in the patches "front" / "back": a periodic pair (``ends="periodic"``) or
two symmetry planes (``ends="symmetry"``). One undo step; returns the new names.

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

## High-order elements

### `high_order(order: int = 2, curved: bool = True, repair: bool = True)`

Order-``order`` elements (1 .. 6) of every mesh cell: triangles and quadrilaterals
(2D), tetrahedra, pyramids, prisms and hexahedra (3D), the same cells as the linear
mesh with equispaced Lagrange nodes, shared nodes merged exactly, and boundary
elements (lines / triangles / quadrilaterals) per patch (with patches the boundary
must validate). ``curved``: boundary nodes moved onto the geometry (curves,
surfaces) and the displacement extended into the neighbouring elements, which are
checked for validity (ValidationError for folded elements; ``repair``: the elements
that fold or lose most of their quality are first optimised, see ``curving``
["repair"]); False: straight-sided elements (straight edges, planar triangles,
bilinear quadrilateral faces). Returns an
ElementMesh (``points``, ``cells`` / ``faces`` per kind in Gmsh node order,
``connectivity(kind, fmt)``, ``element_nodes(kind, nodes)``, ``curving`` report).

### `export_high_order(path, order: int = 2, nodes: str = 'equispaced', base_name: str = 'OpenGrid', curved: bool = True, repair: bool = True) -> dict`

Write ``high_order(order, curved, repair)`` by extension: ``.msh`` (Gmsh 4.1
high-order element types, physical groups, periodic links), ``.cgns`` (SIDS types up
to order 4, one unstructured zone per zone, BCs per patch), ``.vtkhdf`` (VTK Lagrange
cells) or ``.h5``
(element node coordinates in the node set ``nodes``: "equispaced", or "gll" / "gauss"
for quadrilaterals and hexahedra).

### `curved_high_order(names=None, order: int = 2, nodes: str = 'equispaced', snap: bool = True)`

Curved order-``order`` elements of structured domains (2D) or blocks (3D): p x p
(x p) grid cells per element (interval counts must be multiples of p), the grid's
nodes are the equispaced node set, so element edges follow the grid lines and curved
walls. ``nodes``: "equispaced" (file formats), "gll" (Gauss-Lobatto-Legendre:
continuous Galerkin / spectral elements, DGSEM with boundary points) or "gauss"
(Gauss-Legendre: DG solution points); ``snap``: project GLL nodes on boundary sides onto
their curves (surfaces in 3D). Returns a HighOrderMesh (``xyz`` (elements, (p + 1)^d,
3) in tensor order, ``scaled_jacobian``, ``merged()``, ``owners``).

### `export_curved_high_order(path, names=None, order: int = 2, nodes: str = 'equispaced', snap: bool = True) -> dict`

Write ``curved_high_order(...)``: ``.vtkhdf`` (VTK Lagrange cells, equispaced nodes)
or ``.h5`` (any node set: element node coordinates, Jacobians, merged numbering).

## Editing, undo and the model

### `modify(name: str, *, new_name: str | None = None, **params) -> str`

Change some parameters of an entity, e.g. ``og.modify("line-1", p1=(2, 0))``, and /
or rename it (``new_name``) as one undo step; returns the (new) name. A sketch:
edit_sketch, which regenerates it.

### `rename(name: str, new_name: str) -> str`

Rename an entity (references follow; journals stay valid).

### `delete(names, dependents: bool = False) -> int`

Delete entities (a name or a list). What depends on them (the connectors on a
curve, the domains on a connector, the blocks and volumes on a domain, ...) goes too
with ``dependents=True``; without it such a deletion is refused naming them (see
:meth:`deletion_dependents`). Patches, zones and smoothing groups never block: they
lose the members (a group left empty goes, a periodic patch losing its partner is
unpaired). A sketch takes its outputs along; an output on its own is refused (remove
its block with :meth:`remove_sketch_blocks`). Returns the number deleted.

### `deletion_dependents(names) -> list[str]`

The entities that deleting ``names`` would take along (``delete(...,
dependents=True)``), in deletion order; raises for what cannot be deleted at all.

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
