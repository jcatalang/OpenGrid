# Open issues and things worth knowing

A running list of limitations, caveats and design changes collected phase by phase, to review
when the project is complete. Items fixed in a later phase are removed (see `PROGRESS.md` for the
history). Phase tags show where each item came from.

## Not verified
- **CI and other platforms** [0–9c]:
  - GitHub Actions has never run (no remote).
  - The Linux build and the uv-based setup are untested; the conda recipe was built and
    tested on osx-arm64 only (it should build on Linux unchanged, but has not been run there).
- **Real input devices** [1–3]: real mouse and trackpad interaction is untested (handlers are
  tested in world coordinates), and trackpad pinch-to-zoom is not handled.
- **No solver has been run** on any exported file [4]. CGNS is checked with cgnscheck and h5py,
  Gmsh files with Gmsh 4.12 and meshio, VTKHDF with VTK's reader.
- **meshio 5.3.5 cannot read VTKHDF, Plot3D or CGNS** [3, 4], so those spec read-back checks use
  VTK, cgnscheck or h5py instead.
- **Cancel time** [3]: "smoothing cancel returns within 0.5 s" was not timed. Cancellation is
  checked once per iteration (about 1 ms on the demo grids).
## Geometry, connectors, GUI basics
- **GUI geometry** [1]:
  - Geometry is created in the z = 0 plane only (the API accepts 3D).
  - No periodic/closed splines.
  - A modify followed by a rename in the property panel takes two undo steps.
- **nanobind leak report** [1]: disabled at exit (Qt/VTK shutdown order). The API-only suite runs
  leak-free with `OG_LEAK_WARNINGS=1`.
- **Very small far-end spacings** [2]: below about 1e-7·L they can't reach 1e-8 relative
  accuracy (a floating-point limit).
- **Spacing propagation** [2]: one level deep, and only fills tanh/sine ends without a spacing;
  geometric and uniform connectors never adopt spacings.
- **Missing edit tools** [2, 4]:
  - Connector split/join and dimensioning by average spacing are not implemented (the spec
    lists them).
  - Because of that, a connector shared by two boundary types must be split by hand.
- **GUI property panel** [5b]: the status label is clipped.
- **GUI** [2026-10-17 round]: the 2D tools draw in the z = 0 plane whatever the view; 3D
  picking hits drawn lines and faces (a wireframe block is picked on its grid lines, not in
  the open middle of a face); a drag preview rebuilds the merged curve / connector layers on
  its first frame (a noticeable pause on models with thousands of connectors); the busy
  dialog of long operations has no Cancel button except where the operation can stop (the
  smoother); the model tree is rebuilt on every command (fine to a few thousand entities).
- **3D cut plane** [7a]: fixed normal (x) at the model centre.

## Structured 2D
- **Hermite TFI** [3]: can fold C-type grids; the C-grid demo uses plain TFI.
- **Elliptic smoothing (2D)** [3]:
  - Planar z = const domains only.
  - A grid cut is a fixed boundary (no periodic smoothing across it).
  - Wall control only on j = 0 / j = jmax (rotate the edges for other walls).
- **Automatic edge assembly** [3]: doesn't handle C-type loops (a connector on both sides of a
  cut); use explicit edges.
- **Overlapping domains** [3]: the GUI draws them on top of each other.

## Exports
- **2D CGNS** [4]: written with PhysicalDimension 2 (x, y). Some solvers expect 3D data with
  one cell in z; extrude for that (a block from a structured domain, Phase 7; a volume of
  prisms and hexahedra from an unstructured domain, 2026-10-17, with the connectors as side
  walls and ``{"volume", "face": "top" | "bottom"}`` as the ends).
- **Extruded volumes** [2026-10-17]: translate, rotate and normal sweeps only (no hyperbolic
  marching, no smoothing, no inflation layers: the layers are the distribution's); "normal"
  needs a planar domain; a rotation axis through the domain folds it (reported, cells flat or
  inverted); the curved high-order export leaves the side walls straight-sided (a swept
  connector curve is not a surface of the model). The ends pair as periodic faces (back = the
  domain, front = the volume's top, see Periodicity).
- **Plot3D** [5c, 8a]: structured only; it refuses unstructured domains and volumes.
- **Periodicity** [4, 7a, 2026-10-17]:
  - 2D periodic pairs support translation only in CGNS (a 2D base has no component for a
    rotation about z). 3D pairs write rotations as RotationAngle = axis x angle, which readers
    take as a rotation vector (the SIDS wording is "angles about the axes": the same thing for
    an axis-aligned rotation axis; for a general axis it is the rotation-vector reading).
  - 3D periodic faces must match node by node and face by face: the two patches need the same
    mesh (the ends of an extruded volume do; two blocks need matching face dimensions and
    distributions). Mismatches are reported (how many nodes have no image, how far the first
    one lands), never repaired. Block faces pair as CGNS 1-to-1 periodic connections only when
    they map corner to corner; otherwise, and for volume faces, the nodes go into point-list
    connections with the periodic property.
  - Periodic boundaries (2D connectors, 3D faces) are not curved by the curved high-order
    export, so the paired nodes stay images of each other; a periodic boundary on a curved
    surface stays straight-sided.

## Unstructured 2D
- **Delaunay quality near abrupt boundary-spacing jumps** [5a]: not guaranteed (boundary nodes
  are fixed and never split); quality checks run away from such jumps.
- **Advancing front** [5b]:
  - No angle guarantee (angles are measured, not proven).
  - Its fallback to Delaunay for edges it cannot close exists but was never triggered in tests.
  - About 2× slower than Delaunay at 1 M triangles.
- **All-quad meshes** [5c]:
  - Need an even number of segments on each boundary loop (two odd loops with an even total
    are not handled).
  - Loops touching at a point are rejected.
- **Quad-dominant meshes** [5c]: not all-quad (2–8 % triangles remain).
- **Quad quality** [5c]: measured, not guaranteed (worst scaled Jacobian 0.45–0.68 in tests and
  demos).
- **Phase 6 changed Phase 5c output** [6]: all-quad meshes keep sharp boundary corners in their
  coarse boundary now, so 5c meshes differ slightly (NACA demo: 10 756 vs 10 760 quads, min SJ
  0.665 vs 0.683).

## 2D inflation
- **Open wall ends** [6]: layers taper to zero height where a wall connector meets a non-wall
  one, because that connector's nodes are fixed.
- **Layer cells are not smoothed** [6].
- **Worst cells** [6]: the first-layer fan triangles at sharp trailing edges (SJ ≈ 0.36).
- **Internal settings not exposed in Python** [6]: the taper slope (0.5), fan angles and
  proximity factor.
- **Transition on by default** [6]: "stop at core size" means `layers` is a maximum. Columns can
  stop earlier at the core size, near other bodies, or where a cell would fold; turn
  `transition` off for exactly `layers` layers (collision stops still apply).

## Structured 3D
- **Block faces** [7a]: collapsed (degenerate) block faces are not supported.
- **Volume accuracy** [7a]: chord error of the face discretisation. The duct demo is within
  2.9e-5 of analytic at 41 points per 30°; 1e-6 needs ≈ 600 points per 60°.
- **Hyperbolic extrusion** [7b]:
  - An explicit marching variant, not the implicit Chan–Steger scheme.
  - Concave regions rely on dissipation; inverted layers are reported, not repaired.
- **Surface smoothing** [7b]: no wall control.

## Unstructured 3D
- **Pyramids** [8a]:
  - Pyramids on quad faces meeting at sharp edges leave thin tetrahedra between them (minimum
    dihedral about 2.6° on a box).
  - So by default ("auto") pyramids are only used where a structured block shares the face;
    other quads are split into triangles.
- **Faces on surfaces** [8a]: meshed in a parameter plane with a constant metric correction;
  strongly varying surface metrics are only approximated.
- **Optimisation** [8c]:
  - Only interior nodes and pyramid apexes move, and only core tetrahedra are swapped; layer
    cells and boundary nodes are never touched. So the worst cell of a mesh with layers (a thin
    layer transition cell) is usually unchanged.
  - The objective trades cells: on the wing, core tetrahedra below 0.3 drop from 4 745 to 3 691,
    but the worst core pyramid drops from 0.34 to 0.23 (apexes move to help their tetrahedra).
  - Only 2-3 / 3-2 swaps (no general edge removal or 4-4 swaps); compass search, not a
    gradient method. Settings (threshold 0.3, 4 passes) are internal.
  - Netgen's own optimisation still runs first.

## 3D inflation
- **Ridge fans** [8c]:
  - Fans start at a 60° turn (as in 2D), so 90° body edges and box corners now get fans and
    corner caps too. This changed Phase 8b results: the 8b wing went from 43 % to 82 % of wall
    nodes with all 20 layers.
  - A ridge that ends inside a single sector (possible only with a higher threshold than the
    turn of its neighbouring edges) does not grow at its end node ("ridge end"); the layers
    taper there.
  - Corner-cap cells are only checked for positive volume, not for the 0.2 Jacobian floor.
  - Fan columns keep a fixed (radial) direction for all layers.
  - Fans need the turn to be below about 180°; cusps (zero-angle edges) are untested.
- **Thin transition cells** [8b, 8c]:
  - Where columns stop inside anisotropic wall cells, the transition tetrahedra and pyramids
    are thin; a quad split from its centroid gives cells whose quality is bounded by layer
    height / cell width wherever the centre is. On the wing (with fans) 830 layer tetrahedra
    and 39 layer pyramids have a scaled Jacobian below 0.1 (Phase 8b: 3 215 and 715), mostly
    where the layers taper near the tip.
  - Near the wing tip the layers still taper: the tip cap's small triangles reach the core size
    after about 15 layers (transition stop) and the slope rule spreads that onto the wing.
  - Netgen's tetrahedra next to such steps can be poor too.
- **Collision tests** [8b]:
  - Floating point with a 1e-9 tolerance, not exact predicates.
  - Pyramid caps on the layer front are added after the collision tests and never checked
    themselves.
- **Sliding** [8b]:
  - Only along planar unstructured faces used by this volume alone.
  - The volume re-meshes such a face, so its view of the face differs from the face's own mesh,
    which is what the GUI shows.
- **Internal settings** [8b, 8c]: not exposed in projects (as in 2D): minimum Jacobian 0.2,
  normal scale cap 2, slope 0.5, proximity 1, visibility 0.2, 5 smoothing sweeps, fan threshold
  60°, fan step 30°.

## 3D quality floor
- **Layers lost to the floor** [8c addendum]:
  - Thin layers over stretched wall cells can't step down within the floor, so stops spread:
    a wall ending at fixed nodes (a non-sliding neighbouring face) keeps no layers, and ridge
    ends (with a threshold above 60°) stop their neighbours too.
  - The demo wing needed 61 instead of 25 spanwise points (with 25: mean 11 of 20 layers and one
    core cell below the floor).
  - Lowering `quality_floor` trades this back.
- **Core cells** [8c addendum]: cells on a thin triangle of the layer front can only be improved
  through their free node; when none moves enough the volume is reported (it can't be
  exported) rather than repaired.
- **2D has no floor** [8c addendum]: only 3D cells are checked.
- **GUI** [8c addendum]: the floor is editable for volumes; for blocks only through the API.
- **Boundary report cost** [8c addendum]: it now meshes every volume (cached) to check the floor.

## CAD import
- **Faces with poles or apexes get no domain** [9a]: spheres, cones and similar faces have a
  degenerate edge, and meshing in their (u, v) parameters breaks at the pole. They are reported;
  the user builds those domains by hand (e.g. a cubed-sphere layout as in the 08b demo).
- **Surface frames use a constant metric** [8a, 9a]: unstructured domains on curved CAD faces
  are meshed in (u, v) scaled by mean tangent lengths; strongly non-uniform parametrisations
  (some B-spline faces from real CAD) give stretched triangles. Not tested on industrial files.
- **Only test parts were imported** [9a]: OCCT-written STEP / IGES and a netgen.occ part. Real
  CAD exports (assemblies, units other than mm, many tiny faces, sliver faces) are untested; no
  face merging (virtual topology) exists, so every CAD face becomes a domain.
- **Units** [9a]: coordinates arrive in OCCT's session unit (millimetres); use `scale` to
  convert. Names, colours and assembly structure (XDE) are not read.
- **"auto" structured faces** [9a]: only single-wire faces with exactly four edges; a side made
  of several CAD edges is not recognised, and the propagated dimensions can raise counts on
  distant faces.
- **Build-time option** [9a]: `OG_WITH_OCCT` is ON and silently becomes a no-op without
  OpenCASCADE (a configure message says so); `occt` comes into the conda env as a Netgen
  dependency and is now listed explicitly.

## Topology editing
- **Split parts are regenerated** [9b]: the cut follows the old grid line exactly, but each
  part is regenerated by its generator from its new boundaries. Linear TFI on straight-sided
  domains reproduces the old interior; elliptic smoothing, Hermite TFI and curved domains give
  a different (still valid) interior. Attraction lines are dropped on split / join.
- **Six-face blocks only** [9b]: extruded blocks cannot be split, joined or O-gridded (their
  side faces are not domains); faces shared with another block cannot be split (the
  neighbour's face would have to be split too), so neighbouring blocks are split one at a time
  only where the cut does not cross a shared face. Smoothing walls are dropped (face names
  change). Patches naming block faces directly block the operation.
- **O-grid geometry** [9b]: the core is a scaled copy of the corner polygon / polyhedron with
  straight edges; there is no smoothing of the O-grid layout, and domains on surfaces are
  refused (the radial lines would leave the surface).
- **Join leaves curves** [9b]: joined connectors and deleted shared connectors leave their
  database curves (and the polylines of split cuts) in the model.
- **Explicit distributions in the GUI** [9b]: shown as type "explicit"; the positions cannot be
  edited (changing the type or the point count replaces them).

## Scripting and validation
- **Golden cases are coarse** [9b]: sized to run in seconds (whole suite about 13 s); they check
  invariants, not solution accuracy. The NACA O-grid's smallest scaled Jacobian is 0.14 at the
  trailing edge.
- **Journal replay needs the same version** [1–9b]: journals call the API by name; there is no
  versioning of the API itself.
- **Earlier commits lacked the distribution sources** [9b]: a `.gitignore` rule (`dist/`)
  excluded `core/src/dist` and its header directory from Phase 2 on; fixed and committed in 9b,
  so commits before it do not build from a fresh clone.
- **09a_bracket.step changes on every demo run** [9a]: OCCT writes the date into the STEP
  header.

## Performance and packaging
- **3D target missed** [9c]: 5M tetrahedra take 176–186 s (target 120 s; two runs) on an M1 Max:
  Netgen about 157 s even with its thread pool (serial Delaunay insertion about 50 s), then our
  optimiser 13 s, quality and checks about 5 s. Fewer Netgen optimisation passes would meet the
  time but drop the worst scaled Jacobian from 0.17 to 0.015, so they were not used. Partitioned
  parallel meshing or another backend would be needed.
- **Serial parts** [9c]: Delaunay / advancing-front insertion, 2D and 3D inflation, the
  optimiser's node moves and swaps (order-dependent; 13 s of the 5M case) and Netgen's insertion.
  The 2D targets are met anyway (1M cells in 1.8–4.8 s including export).
- **"Converges interactively"** [9c] was read as: converges to the default tolerance in under
  2 s with live progress. Measured on a 500 x 200 O-grid with wall control: 1.06 s.
- **Smoothed grids changed slightly** [9c]: zebra line ordering (needed for parallelism) and
  Anderson acceleration (default depth 5) reach the same fixed point by a different path, so
  smoothed grids differ from earlier versions at the tolerance level (about 1e-3 of the local
  spacing on the test O-grid). Projects saved earlier get `acceleration: 5` when loaded;
  `acceleration: 0` gives plain relaxation (but still zebra-ordered).
- **Determinism across threads** [9c]: checked (identical digests) for 1–10 threads on the
  elliptic smoothers, volume meshes and quality metrics; Netgen's thread pool gave identical
  meshes in our checks, which Netgen itself does not promise.
- **VTKHDF export dominates 2D timings** [9c]: about 1.7–2.4 s of each 1M-cell case
  (HDF5 writes from Python; not parallelised).
- **App bundle** [9c]: arm64 only; signed ad hoc, not notarised (Gatekeeper asks on first launch
  of a downloaded copy); large (a whole conda environment with Qt, VTK, OCCT, Netgen). Building
  it needs conda, rattler-build and conda-pack, and network access to conda-forge.
- **No AppImage** [9c]: the conda package is the Linux deliverable.
- **Package metadata** [9c]: version 0.0.1, license `LicenseRef-OpenGrid` and a placeholder home
  page; set them before publishing.

## Multi-block smoothing (Phase 10)
- **Laplace control only in groups** [10a]: no Thomas-Middlecoff or attraction functions.
  Clustering survives only at walls (first-cell height control); clustering on other
  boundaries (e.g. a refined far field) is smoothed out. Make such connectors walls to keep
  their first spacing.
- **Sliding is local** [10a]: a sliding node stays within half its initial intervals of its
  starting position (unbounded sliding piled nodes into 90° boundary corners: min scaled
  Jacobian 0.08 on the cylinder demo). So orthogonality on sliding boundaries is improved,
  not reached: 2.6° error on the bump channel's top wall (5.5° fixed). End nodes of boundary
  connectors slide only with `joints="slide"` (10d: between two connectors on one curve in
  2D, inside one surface in 3D); otherwise an interface ending on a boundary keeps its end.
- **Topology decides the worst cells** [10a]: on the cylinder demo the group smooths the
  diagonals (kink 90° → 1°) but the worst cell drops from 0.45 to 0.36, in the channel
  corners, where a 90° boundary corner is split between two domains. Multi-block layouts
  give such corners a single block; templates (10e) will.
- **Kink near irregular points** [10a]: grid lines must turn there; the kink metric ignores
  three nodes around them (up to 13° there on the circle O-grid).
- **Group default wall spacing** [10a]: None keeps the initial first-cell height at each node
  (the single-domain smoother interpolates it from the side boundaries).
- **Planar 2D groups** [10a]: 2D groups need domains in a z = const plane (blocks: 10b).
- **Coarse-to-fine smoothing not built** [11c]: grid sequencing was prototyped on the 1M-cell
  ball (every other node, displacement interpolated to the fine grid): 1.5x fewer fine
  iterations at best, a coarse level that stalls itself, and a different result (slide limits
  are relative to the start). Group runs often end at a residual floor (relative ~1e-6, e.g.
  the 33-point ball runs to the 2000-iteration limit) that no coarse start removes; likely
  sliding nodes clamped at their slide limits. Neither a stagnation criterion nor multigrid
  was built.
- **Curved high-order elements: structured grids only** [11e]: `curved_high_order` makes
  quadrilaterals and hexahedra from structured domains and blocks; no curved triangles,
  tetrahedra, prisms or pyramids (curving boundary-layer prisms needs the inner layers curved
  too, and untangling), and the curved elements go to VTKHDF and HDF5 only (not Gmsh / CGNS).
  Interval counts must be multiples of the order (`interval_multiple` on sketches and
  templates); elements inherit the grid's quality, and the curved elements' scaled Jacobian
  can be low at a sharp trailing edge (0.064 on the C-grid at order 4). GLL snapping in 3D
  moves face interiors onto their surfaces but leaves nodes on block edges as interpolated;
  the equispaced geometry between grid nodes is the degree-p interpolant, not the CAD itself.
- **Combinations not covered** [12b]: the CAD bracket, wing_in_box / wing stub (10 s+ each)
  and rotational periodicity through the exports are not among the combination cases (3D
  periodic patches have their own tests, `tests/python/test_periodic_3d.py`); nor is the GUI
  (12a covers its interactions on a small model).
- **High-order overlay cost** [12a]: the overlay recomputes the whole high-order mesh on every
  model change (0.5 s for the C-grid template at order 4, 10 s for a 32k-hexahedron ball);
  switch it off while editing large models. It draws the patches' faces in 3D, or every cell
  edge without patches. Pyramids' curved edges are drawn (VTK files write them linear).
- **Viewport picking** [12a]: a click picks connectors before curves and domains; curves under
  a connector are edited through the connector's handles or the model tree.
- **Trackpad gestures** [GUI]: in 2D navigation the zoom window is a left-drag; a three-finger
  drag is one only with macOS's "Three Finger Drag" (System Settings -> Accessibility ->
  Pointer Control -> Trackpad Options); in 3D navigation it is Alt + drag and a plain drag
  rotates. Two-finger scroll zooms (as a wheel); pinch gestures are not mapped by OpenGrid. The launcher app depends on the `opengrid` conda environment (moving or removing
  it breaks the app; rerun the installer). A running OpenGrid keeps the code it started with
  (the status bar says when it changed: restart); the self-contained bundle in packaging/dist
  is not rebuilt with the code (build_app.sh, on request).
- **Curved boundary elements** [11g, 12c]: `high_order` curves the boundary by default, repairs
  the elements that fold or lose most of their quality by moving their free high-order nodes,
  and checks the low ones with Bezier bounds. Repair moves no linear vertex and no node on the
  geometry, so elements whose straight cells are nearly flat under a curved face (one block
  wrapped onto a sphere) stay folded and are reported; the O-grid version does not fold. A
  coarse tetrahedron on a curved layer stack keeps a loss (0.97 -> 0.70 after repair). The
  extension is a stiffened graph Laplacian (weights 1 / length^3, chosen on four cases), not
  elasticity, over a region of 3 boundary-element lengths; CG (tolerance 1e-10) above 20k
  unknowns, 92 - 900 iterations seen. Faces without a surface are curved from their own grid
  (degree-p interpolation: the pipe wall 2.6e-7 off its cylinder at order 4, against 1e-11 for
  projected surfaces); give curved faces a surface for exact nodes. Periodic patches and
  internal walls stay straight; VTKHDF draws pyramids linear (curved ones too). Projection is
  the closest point on the connector's curves or the face's surface (no parameter-interval
  check). Validity is rigorous (Bezier bounds, degree elevation) for elements whose sampled
  scaled Jacobian is below 0.3 and through order 4; hexahedra of orders 5-6 and pyramids are
  sampled only, and elements sampled above 0.3 are not bound-checked (a fold between samples
  there is not excluded).
- **Straight-sided element details** [11f]: quadrilateral faces are bilinear, planar only
  when their four vertices are coplanar (hexahedra of smooth grids have slightly warped
  faces). Orders 1-6 (CGNS defines up to 4). GLL / Gauss node sets for quadrilaterals and
  hexahedra only (HDF5); no simplex node sets (warp-and-blend, Fekete, Witherden-Vincent):
  solvers make those from the vertices. VTKHDF writes pyramids as linear pyramids (VTK 9.7 has
  no Lagrange pyramid; exact geometry, their interior nodes left out). CGNS: a 2D base takes
  translational periodicity only (as the linear 2D export), a 3D base the rotations too;
  PENTA_40 / PENTA_75 node orders follow the SIDS
  descriptions (Gmsh does not write them; checked against all its statements, not against a
  second implementation). Gmsh 4.15 cannot evaluate prisms of order >= 3 itself, and its
  multithreaded Jacobian evaluation crashed on the 35k-cell order-4 file (single-threaded it
  reads and evaluates it); a solver's reader is the real check. The element mesh covers the
  whole model (no subset); with patches the boundary must validate.
- **GPU smoothing does not pay yet** [11d]: the experimental OpenCL path (`opengrid.set_gpu`,
  off by default) runs the 3D group smoother's line sweeps 2.2x faster than the CPU on the 1M-
  cell ball, but the whole iteration only 1.2x (Anderson mixing and sliding projections stay
  on the CPU), its colouring converges more slowly (Jacobi within a colour), and its float32
  residual floor (~1.3e-6) is above the default 1e-6: GPU runs stop at 1e-5 and differ from
  the CPU grid by up to 7 % of the local cell size. Short lines (the 10d wing stub: ~9 nodes)
  are latency-bound on the device (no gain). OpenCL is deprecated on macOS (Metal is the
  native API); only the 3D group smoother has a GPU path; tested on one device (M1 Max).
- **3D speed target missed** [10b, 11b, 12d]: a 1M-cell group of blocks takes 18.5 s after
  12d (27 s before 11b; target 10 s): 328 iterations at 56 ms, half of it in the line sweeps
  and a third in the Anderson mixing. The sweeps are bound by scattered memory access along
  the j and k lines (about 1.5 KB of cache-line traffic per node per iteration), the mixing
  by the traffic of its history (depth 20 over 3.2 M unknowns; single precision now).
  Depth, over-relaxation, coarse-to-fine (11c) and the GPU (11d) were all measured and do
  not pay; a per-direction data layout for the sweeps was not built. Single-threaded speed
  is the limit (the 10d wing stub: 43.7 s on one thread, 9.2 s on ten). `OG_PROFILE=1`
  prints the phase times of a 3D group run.
- **3D wall grids are not smoothed** [10b]: wall faces stay fixed, and wall control copies
  their grid into the first layers. On the sphere in a box the cubed sphere's own grid turns by
  36° at its patch seams, so grid lines near the wall kink by up to 31° there (1.5° without
  wall control). The 3D kink metric also leaves out three lines along each interface face's
  edges.
- **Box corners** [10b]: the sphere-in-a-box layout gives each box corner to three blocks;
  their corner cells drop to a scaled Jacobian of 0.33 (0.58 per block). As in 2D, the layout
  decides the worst cells (a box corner should belong to one block).
- **Sliding in 3D is local too** [10b]: a node moves at most half its smallest initial
  neighbour distance (in parameter space); sliding faces need a surface entity or a plane
  (curved faces without a surface stay fixed). Seams between faces on different surface
  entities of one smooth surface (the cubed sphere) slide along their connector, not across.
- **Extruded blocks** [10b]: not allowed in groups (six-face blocks only).
- **Real mouse dragging untested** [10c]: drags and strokes were driven through the
  world-space handlers offscreen; the VTK observers that take the left button from the camera
  ran only in the screenshot session (no manual check with a mouse or trackpad).
- **Handles** [10c]: only for selected points, lines, polylines, splines, NURBS control
  points and connector break points; arcs, surfaces and NURBS weights have none. Grid nodes
  are never dragged (they are derived).
- **Curve networks** [10c]: planar curves only; a structured face needs four corners
  (turns > 30 degrees); opposite sides made of several network edges are not balanced (the
  face is reported if the counts differ); a smooth closed curve (no corners) is not cut into
  four sides; meeting points of three curves follow only the first two when a curve is
  dragged. "Block from curves" (3D) was not built.
- **Freehand** [10c]: fitted once (later drags move control points); pointer samples in the
  z = 0 plane; corners below 60 degrees are smoothed over. A stroke is closed into a loop only
  when released within the pick tolerance (8 pixels) of its start [2026-10-17]: a loop whose
  end misses its start by more is open, and "Unstructured domain" says which end meets
  nothing; the closure is C0 (the loop may have a small kink where the stroke began).
- **Adding and removing curve points** [2026-10-17]: a point added to an interpolating
  spline lies on the old curve, but the new spline through one more point differs slightly
  elsewhere (a few per cent of the chord); NURBS control points are added by knot insertion
  (exact) and removed by Tiller's knot removal, which is exact only when the curve is smooth
  enough there (otherwise the curve moves and the log says by how much); the control points
  without an interior knot (the first and last two of a cubic) cannot be removed, they fix
  the end points and tangents. Lines and arcs have no insertable points.
- **Curve networks: rings** [2026-10-17]: a curve joining an outer loop to an inner one makes
  the region between them a face whose boundary runs along that curve twice; it is reported,
  not meshed (delete the joining curve and mesh the loops as a domain with a hole, or select
  the connectors and use "Unstructured domain").
- **Sketch vertices on curves in 3D stay where drawn** [10d]: 3D vertices slide on planes and
  surfaces only. Sliding along an edge between two surfaces (e.g. where the wing stub's
  outflow meets a side wall) was built (one chain per curve, as in 2D) but made the 3D solver
  stall (1500 iterations without converging in some passes) and was removed. The wake blocks
  that suffered from it (0.12) reach 0.62 since the side walls are symmetry planes (10e
  follow-up); the outflow vertices still keep their drawn positions. Not retried.
- **Sketch boundary vertices without geometry stay where drawn** [10d] (2D boundary edges
  with no curve, box corners, curve ends); the result depends on the drawing only through
  them, the assigned geometry and the topology. Vertices on pinned edges and internal walls
  (a wake cut) stay where drawn too: the passes have nothing to move them by.
- **Vertex passes cost** [10d]: each pass smooths the group on a coarser copy (half the
  intervals); 6-15 passes per edit: 2.7 s for the airfoil (3 edits), 44 s for the wing stub
  (159k nodes). Previews skip them; set `relax=False` while setting a sketch up.
- **Sketch assignment by proximity** [10d]: `assign_sketch_geometry` assigns the nearest
  given curve within 20 % (ends) / 50 % (middle) of the edge length; a sketch drawn very
  loosely can get an edge assigned to the wrong curve (seen with a 50 % end tolerance on
  the wing: a mid-chord spanwise edge took the trailing-edge line). Check the report.
- **Wall sources restart** [10d]: warm-started group runs restart wall control from zero
  (only positions are warm-started); the sketch passes run without walls.
- **Sketch connector adoption** [2026-10-17]: 2D edges only; a sketch cannot adopt another
  sketch's connectors (two sketches do not share a boundary); the vertex passes use a
  sketch-owned copy of an adopted connector at half resolution; adopted connectors are kept
  fixed by the group (their nodes belong to both sides), so the sketch's vertices at their
  ends do not relax; a wall on an adopted edge gets wall control on the sketch's side only.
- **Unstructured sketch blocks** [2026-10-17]: 2D only; all unstructured blocks of a sketch
  share one generator; a wall on an unstructured block needs a first-cell height (its layers
  take the sketch wall's layers and growth; the structured wall's decay / relaxation have no
  unstructured meaning); the vertex passes relax the topology as if every block were
  structured; an unstructured block does not take part in the smoothing group, so a
  connector between the two kinds stays where the connector puts it.
- **Re-dimensioning** [2026-10-17]: when the opposite edge of a structured domain is a
  chain of connectors, the longest one takes the change (a warning names it); a connector
  with explicit breakpoints past the new count is refused by its own validation; a change
  that reaches a sketch's connector is refused (set it in the sketch); an over-constrained
  loop (every opposite connector already fixed and still unequal) is refused naming the
  domain.
- **Deleting** [2026-10-17]: removing a block from a sketch renumbers the remaining blocks,
  so a domain's name can shift (``s-d1`` becomes ``s-d0``); 3D patches naming such a domain
  by id keep the entity, which now stands for another block. A sketch's connector or curve
  cannot be deleted on its own (delete the domain of its block, or the sketch). Blocks can
  only be removed from 2D sketches. Deleting a curve a sketch edge is assigned to is refused
  (unassign first).
- **Boundary panel** [2026-10-17]: a 2D model's report is recomputed at every change (the
  connector nodes are cached, so it is cheap); a 3D model's report meshes every block and
  volume, so the panel only marks itself stale and waits for Refresh. Empty patches are kept
  (listed as such) rather than deleted when their last member leaves. Grid snapping of
  clicks is 2D only (a 3D click has no depth on a grid plane).
- **Extruded 2D model** [2026-10-17]: along z only, one thickness, the end patches are named
  front / back (refused when those names are taken); a structured domain's edge must lie in
  one patch (a block face carries one patch: split the domain where the type changes);
  structured domains become blocks without smoothing, unstructured ones volumes of prisms
  and hexahedra; the 2D domains stay as the fronts (deleting them deletes the extrusion).
- **Suggested values** [2026-10-17]: the clustering spacings (a quarter of the uniform
  spacing), the growth ratio (1.2) and the first layer height (a twentieth of the boundary
  spacing) are starting points, not y+ estimates: no flow data is known; the inflation hint
  takes the mean boundary spacing of the whole domain (volume: of its faces), not of the
  wall connectors alone.
- **Domain conversion** [2026-10-17]: structured -> unstructured keeps no generator settings;
  unstructured -> structured needs one loop whose connectors assemble into four edges with
  matching counts (the usual edge assembly), holes are refused; the entity id changes
  (zones and 3D patches follow, blocks and volumes using the domain refuse the conversion).
- **Sketch tool attachment** [2026-10-17]: a block drawn onto an existing sketch's vertices
  joins that sketch; two separate sketches never share connectors (a sketch adopts only
  connectors made by hand).
- **Reference grid** [2026-10-17]: an axis view shows the plane facing it, a tilted view the
  three axis planes through the round point nearest the view centre (the lines on the world
  axes are coloured only when that point lies on them); the step follows the view's extent
  (1, 2, 5 x 10^k with 10, 4 or 5 minor divisions), not the model's units; the edge labels
  of a perspective view are read at mid depth.
- **Sketch GUI** [10d]: blocks are added four clicks at a time (2D, in the z = 0 plane);
  vertices, blocks and assignments are otherwise edited by script or the property panel
  (dimensions, walls, relax); no 3D axis handles; deleting single blocks has no tool.
- **Sketch kink metric** [10d]: the airfoil group reports a 22° interface kink three cells
  from the trailing edge, where the C-lines turn by the 14° wedge (8° six cells away).
- **Wing in box: the 3D-smoothed variants differ from the 2D grid** [10e, follow-up]: the
  template now extrudes the smoothed 2D C-grid (min scaled Jacobian 0.62, target met). The
  six-face variant (`three_d_blocks=True`, and taper / sweep, which need it) smooths in 3D with
  the end planes as symmetry planes: 0.48 (taper 0.6: 0.55, sweep 20°: 0.48, both: 0.60), and
  its sections differ from the 2D grid by up to 0.24 in the far-field blocks (0.07-0.10 near
  the outflow): the block corners on the end planes keep the 2D-relaxed positions and the
  far-field edges slide on the extruded arc and lines within their slide limit. Not found:
  why the end-plane sections are not exactly the 2D solution.
- **Symmetry planes** [10e follow-up]: planar faces only (a flagged face on a curved surface
  is an ordinary sliding face), listed explicitly (`symmetry`; templates set their end planes);
  2D groups have none (their boundaries slide as before).
- **Wall orthogonality without wall layers is controlled at the first cell only** [10e
  follow-up, measured]: Steger-Sorenson sources with Hilgenstock feedback make the first-cell
  height and the angle of the first grid line exact at convergence; further out the lines bend
  (1e-5 wall: mean 2.4° at layer 8, 10° at layer 16, 19° at layer 29 on the C-grid). Wall
  layers (11a) keep them normal: 0.05-0.09° mean through 30 layers.
- **Wall layers: corners between two walls are refused** [11a]: where two walls with layers
  meet at a convex corner of the domain (a duct or cavity corner) their lines cross in the
  corner region; this is reported as overlapping layers, not built (no corner blending).
- **Wall layers: front ends lean** [11a]: where a wall with layers meets a non-wall boundary
  (an inflow, a far field) the last line sees one facet and leans by up to about half a cell's
  turning angle (2.25° on the test arc; 0 on straight walls).
- **Wall layers: spacing steps** [11a]: neighbouring walls with different spacings (or layer
  counts) step the layer front where they meet; with the C-grid's default 2x wake spacing the
  iteration diverged, so templates with layers give the wake cut the wall spacing. No height
  blending across a junction.
- **Wall layers: what they replace** [11a]: inside the layers the heights are the geometric
  series spacing * growth^k, not the connectors' own distribution; the layers are fixed during
  the smoothing (they do not adapt to the outer grid); lines need about layers + 20 nodes for
  a good outer grid (C-grid at 1e-5: 30 layers with 31 inner radial points gave 0.21, with 61
  0.47). The deviation reported in the layers is the angle to the central-difference tangent;
  at the leading edge, where the front turns 14° per node, it reaches 1.1°. Single-domain
  elliptic smoothing (outside groups) has no layers.
- **Sharp trailing edge** [10e follow-up]: C-grid lines turn from the wedge faces (±8°) into
  the wake within ~0.005 chords of the trailing edge: interface kink 20° at a 1e-3 wall
  spacing, 53-60° at 1e-5 - 1e-6 (the group status reports it).
- **Accelerated smoothing can diverge or stall** [10e follow-up, 11a]: undamped Anderson mixing
  blew up on the C-grid at wall spacings <= 1e-5 (and returned the TFI grid), and stalled near
  1e-2 before jumping to NaN on layered grids. Now guarded: restart from the best iterate on a
  blow-up, a non-finite residual or 200 iterations without a 10 % drop; plain iteration after
  two restarts (warning). The plain iteration needs more iterations (C-grid at 1e-6: ~2100-
  3100). 3D wall control alone at 1e-5 (the tapered wing) does not converge in 4000
  iterations (min scaled Jacobian 0.43); with wall layers it converges in 1770.
- **Automatic topology: two-element airfoil fails** [10e]: the cross field is right (index sum
  -1) but in the slot between main element and flap the separatrices from the flap's leading
  edge and the main trailing edge leave a five-cornered face and a triangle. No repair of
  faces with other corner counts exists; the case is a strict xfail in the tests.
- **Automatic topology is heuristic** [10e]: singularities next to the boundary are snapped
  onto it (within 4 local mesh sizes; -1/4 joining a junction, +1/4 a convex-corner mark),
  boundary separatrices start straight and evenly spaced; tuned on the four test cases. The
  layouts are valid but not optimised: the bracket's fillet corners make thin blocks (min
  scaled Jacobian 0.13); no layout simplification (merging thin blocks, moving singularities)
  is done. The field is harmonic plus a unit-length smoothing, not an optimised frame field.
- **Templates fix their boundary vertices** [10e]: sliding them drifted (C-grid sides,
  towards the trailing edge) or had a neutral mode (rotation along a circle); the C-grid's and
  O-grid's vertices are where the template puts them. Channel templates relax (2D).
- **Template GUI** [10e]: O-grid, C-grid and butterfly from the selection; the channel and 3D
  templates are scripting only (they need the four sides named).
- **Groups and topology edits** [10a, 10b]: members follow splits, joins and O-grids; a 2D
  group's member order changes when an O-grid replaces a member.
- **Picture bug in 9b** [10a]: `09b_topology.png` coloured cells by the wrong quality values
  (structured cell order); fixed and regenerated in 10a.

## Design changes and deviations from the spec or the plan
- **Phase splits:**
  - Phase 5 → 5a (Delaunay), 5b (advancing front), 5c (quads, smoothing, CGNS unstructured).
  - Phase 7 → 7a (foundation), 7b (smoothing, hyperbolic extrusion).
  - Phase 8 → 8a (volumes), 8b (3D inflation), 8c (optimisation, ridge fans, mixed models).
  - Phase 9 → 9a (CAD import), 9b (scripting, topology helpers incl. O-grid / butterfly for
    multi-block hexahedral meshes), 9c (performance, packaging).
- **Inflation scope** [6]: inflation is part of the unstructured domain generator, not a separate
  entity.
- **Hyperbolic extrusion** [7b]: explicit marching scheme (see above).
- **Tetrahedral backend** [8a]:
  - `TetBackend` is a Python class driving Netgen's Python module, instead of the C++ interface
    planned in ARCHITECTURE §7. That is Netgen's supported interface, and it keeps Netgen an
    optional run-time dependency.
  - `netgen` was added to `environment.yml`.
- **3D inflation is a volume generator setting** [8b], like 2D inflation, not a separate
  entity. Walls are given as faces (domains); the spec says "wall patches".
- **Optimisation is our own, not mmg** [8c]: mmg3d handles tetrahedra only (our meshes mix
  prisms, hexahedra and pyramids) and isn't packaged for osx-arm64. Spec and ARCHITECTURE had
  planned mmg for optimisation.
- **Environment** [8a]: Netgen installed into the `opengrid` conda env (VTK, PySide6, pyvista
  and NumPy unchanged).
