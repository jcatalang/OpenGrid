# Open issues and things worth knowing

A running list of limitations, caveats and design changes collected phase by phase, to review
when the project is complete. Items fixed in a later phase are removed (see `PROGRESS.md` for the
history). Phase tags show where each item came from.

## Not verified
- **CI and other platforms** [0–9b]:
  - GitHub Actions has never run (no remote).
  - The Linux build and the uv-based setup are untested.
- **Real input devices** [1–3]: real mouse and trackpad interaction is untested (handlers are
  tested in world coordinates), and trackpad pinch-to-zoom is not handled.
- **No solver has been run** on any exported file [4]. CGNS is checked with cgnscheck and h5py,
  Gmsh files with Gmsh 4.12 and meshio, VTKHDF with VTK's reader.
- **meshio 5.3.5 cannot read VTKHDF, Plot3D or CGNS** [3, 4], so those spec read-back checks use
  VTK, cgnscheck or h5py instead.
- **Cancel time** [3]: "smoothing cancel returns within 0.5 s" was not timed. Cancellation is
  checked once per iteration (about 1 ms on the demo grids).
- **3D performance target** [8a]: "5 M tetrahedra < 2 min" not measured. The 118 k-tetrahedron
  volume takes about 1 s in Netgen plus quality, but most of the 8a demo's 11 s is surface
  meshing and shell assembly in Python.

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
  one cell in z; extrude with a block (Phase 7) for that.
- **Plot3D** [5c, 8a]: structured only; it refuses unstructured domains and volumes.
- **Periodicity** [4, 7a]:
  - 2D periodic pairs support translation only in CGNS.
  - Periodic pairing of 3D faces is not supported.

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
- **3D elliptic smoothing** [7b]: point SOR only, so it is slow on large blocks (about 850
  iterations for 18 000 cells).
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
