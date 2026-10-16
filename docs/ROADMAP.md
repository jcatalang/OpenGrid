# OpenGrid — Roadmap

Each phase ends with passing tests, a demo under `examples/`, an output under
`examples/output/`, an entry in `docs/PROGRESS.md`, and a stop for review.
Acceptance criteria below are automated (pytest / GoogleTest) unless marked *(manual)*.

---

## Phase 0 — Kickoff and scaffold
- `cmake --preset release` + build + `ctest` pass; the GoogleTest suite has ≥1 test.
- `pip install --no-build-isolation -e .` succeeds; `import opengrid._core` works and `add(2, 3) == 5`.
- GUI smoke test (pytest-qt, offscreen): the main window builds; the model tree has the top-level nodes
  Database, Connectors, Domains, Blocks, Patches, Zones; the console evaluates `1+1`.
- `ruff check`, `clang-format --dry-run --Werror`, `tools/check_no_ternary.py` are clean.
- CI workflow exists for macos-14 and ubuntu-latest. *(Running it needs a remote.)*

## Phase 1 — 2D geometry kernel + project I/O
- Curves: point, line, arc, polyline, cubic spline (natural/clamped/Catmull-Rom), NURBS.
- The arc length of a circle arc and a spline is within 1e-10 relative of the reference (analytic / adaptive Gauss–Kronrod at tight tolerance).
- After arc-length reparametrization, `|s(t_k) - k·L/N|/L < 1e-10`.
- Closest-point projection: the distance error is < 1e-12·L for points on and off curves, including near endpoints and high curvature.
- Every command has an undo/redo test: apply → undo gives a store equal to the initial one; redo gives a store equal to the applied one.
- Project round trip: save → load → save is byte-identical. Journal replay reproduces the same `.ogp`.
- Airfoil `.dat` (Selig and Lednicer) and CSV import. Malformed files produce reports with line numbers.
- Demo `examples/01_geometry.py`: NACA 0012 + farfield circle → `examples/output/01_geometry.ogp`.

## Phase 2 — Connectors and distributions
- Distributions: uniform, geometric, tanh and sine with one- and two-sided spacing (Vinokur), breakpoints/sub-connectors, arc-length vs parametric mode.
- Achieved end spacings vs requested: relative error < 1e-8 over a parameter sweep.
- Node positions are strictly monotone. The max stretching ratio is reported and checked against a threshold.
- Extreme ratios (e.g. spacing 1e-12·L, or impossible requests) fall back with a warning and never produce NaN/Inf. This is fuzz-tested.
- Spacing propagation across shared nodes: changing one connector updates its flagged neighbours, and undo restores them.
- Demo: airfoil connectors with LE/TE clustering.

## Phase 3 — Structured 2D domains + first export
- Domains from 4 connector-chain edges, with auto-assembly from a selection and orientation handling. A count mismatch produces a report naming the edges and counts.
- Linear TFI is exact on the quarter annulus (interior nodes vs analytic polar grid when edges are arcs/lines with matching distributions, error < 1e-12). Hermite TFI is also provided.
- Elliptic: Winslow and Poisson (TTM) with Steger–Sorenson wall control. The residual drops ≥ 6 orders on the golden cases. Wall spacing is achieved within 5% and wall orthogonality to within 2°.
- Worker-thread run with progress and cancel. Cancel returns within 0.5 s.
- Quality metrics (2D set). All Jacobians are > 0 on golden cases.
- VTKHDF export read back by `vtk` and `meshio` with matching cell counts and coordinates. Plot3D 2D multiblock round trip.
- Demos: NACA 0012 C-grid and O-grid with prescribed wall spacing; backward-facing step multiblock.

## Phase 4 — Patches, zones, solver-ready export
- Validator: every boundary connector is in exactly one patch, interfaces are auto-detected, and violations are reported by name.
- Structured 1-to-1 connectivity detection between domains, with point-range transforms verified.
- VTKHDF with patches (collection or fallback), CGNS (zones, 1-to-1, BC_t) validated by `cgnscheck` if available, and Gmsh 4.1 read by meshio with the correct physical groups.
- The patch face counts in exported files equal the model counts.

## Phase 5a — Unstructured 2D: constrained Delaunay + refinement (split approved after Phase 4)
- Unstructured domains from an outer loop and inner loops of connectors (automatic loop assembly
  and orientation); boundary nodes are exactly the connector nodes (no boundary Steiner points).
- Robust orient2d/incircle (filtered + exact expansion arithmetic), Bowyer–Watson insertion,
  segment recovery by flips, constrained Lawson flips, domain classification by flood fill.
- Delaunay refinement (circumcenters, rejected when they would encroach on a boundary segment)
  driven by a size field: boundary spacing with growth rate, maximum size, point/line sources.
- Tests: every boundary segment present exactly once, total area equals the geometric area
  (1e-12 relative), no inverted triangles, min angle >= 20.7 degrees away from small input
  angles, determinism (identical output on repeated runs), conformity with structured neighbours.
- Exports: triangles in VTKHDF and Gmsh (CGNS unstructured zones later).
- Demo: NACA 0012 in a circular farfield (triangles), and a hybrid structured O-ring + triangles.

## Phase 5b — Advancing-front triangulation (alternative generator)
- Same input/output contract as the Delaunay mesher: boundary nodes preserved, every boundary
  segment present exactly once, exact area, no inverted triangles, deterministic output.
- Min angle ≥ 25° on uniform and graded test domains; median edge/target size within 5% of 1
  on a uniform disc.
- Never fails silently: an edge the front cannot close is handed to the Delaunay mesher, with
  a warning.
- Selectable per domain (`generator.algorithm`) in the API, projects, journal and GUI.

## Phase 5c — Quad-dominant / all-quad, smoothing and edge swapping; CGNS unstructured zones
- Element types `triangles`, `quad_dominant`, `quads` per domain (API, projects, journal, GUI).
- Every mesh: boundary nodes are the connector nodes, each boundary segment present exactly once,
  exact area, scaled Jacobian > 0 for every cell, deterministic.
- All-quad: no triangles at all. A loop with an odd number of segments is reported (an all-quad
  mesh of it does not exist).
- Smoothing (Laplacian, length-weighted, angle-based) improves mean quality and never lowers the
  worst cell; edge swapping; untangling repairs inverted cells.
- CGNS: unstructured zones (element sections, BCs as BAR_2 sections at EdgeCenter) and vertex
  connectivity with structured zones and for periodic pairs; cgnscheck clean.
- 1 M cells < 10 s for every element type.

## Phase 6 — 2D inflation (hybrid)
Implemented inside unstructured domains (`generator.inflation`); the core is filled by the
Phase 5 generators (triangles, quad-dominant or all-quad).
- First-cell height and growth rate achieved (relative error < 1e-6 on straight walls, < 2% on curved walls).
- No overlaps (exact segment-intersection checks), positive Jacobians, and sharp and blunt TE handled.
- Proximity/collision termination works (two-element airfoil gap). The last-layer size matches the core size field within 20%.
- Demos: NACA 0012 hybrid; two-element airfoil.

## Phase 7 — 3D structured *(split after Phase 6: 7a foundation, 7b smoothing + hyperbolic extrusion)*
- **7a:**
  - surfaces with projection;
  - domains on surfaces / in 3D;
  - blocks from six faces (automatic orientation) and by translate/rotate/normal extrusion;
  - 3D TFI and quality;
  - 3D topology and patches;
  - exports (VTKHDF hexahedra, CGNS 3D, Plot3D 3D, Gmsh);
  - GUI cut plane / shrink;
  - demos: annular sector duct, extruded wing section.
- **7b:**
  - surface elliptic smoothing;
  - 3D elliptic smoothing with wall control;
  - hyperbolic extrusion;
  - demo: 3D cylinder O-grid.
- Surfaces: planar, ruled, revolved, extruded, NURBS, with projection accuracy < 1e-10·L.
- Blocks from 6 faces, 3D TFI, 3D elliptic smoothing with wall control, hyperbolic extrusion, and translate/rotate/normal extrusion.
- Block volume vs analytic (annular sector) within 1e-6 relative for fine grids. Positive Jacobians. Consistent 1-to-1 interfaces.
- Export VTKHDF (hexes), CGNS 3D structured, Plot3D 3D.
- Demos: annular sector duct; extruded wing section; 3D cylinder O-grid.

Phase 7 acceptance status: volume vs analytic, positive Jacobians, 1-to-1 consistency and all exports verified (7a); smoothing and hyperbolic extrusion (7b).

## Phase 8 — 3D unstructured and hybrid *(split: 8a surface mesh + shell + backend, 8b inflation, 8c optimisation + mixed blocks)*
- **8a:**
  - unstructured domains in 3D planes and on surfaces;
  - volumes from closed shells (watertightness, orientation, inner bodies);
  - tetrahedral backend interface with Netgen;
  - pyramids at quad faces;
  - tet/pyramid quality;
  - 3D unstructured exports (VTKHDF, Gmsh, CGNS);
  - GUI.
- **8b:** prism/hex inflation from wall faces of volumes:
  - marching normals with smoothing and visibility;
  - termination (layers, thickness, transition to the core size, slope, proximity, collision,
    quality);
  - transition cells (pyramids, tetrahedra) where columns stop;
  - layers sliding along planar faces (e.g. symmetry planes), re-meshed around the layers;
  - pyramids on the quad front;
  - export of mixed cells.

  Criteria: first height and growth achieved (exactly on flat walls), positive cells, volume
  equal to the enclosed volume, element-type counts. Demos: sphere in a box, hybrid wing stub.
- **8c:** optimisation (mmg or own), mixed structured/unstructured models in CGNS, pyramid
  quality at sharp edges. As built:
  - own optimiser for mixed meshes (node moves of interior nodes and pyramid apexes, 2-3 / 3-2
    swaps of core tetrahedra), generator key `optimize` (on by default);
  - ridge fans in 3D inflation (convex edges turning more than 60°: fan columns, zero-area fan
    faces, corner caps), which replace the stops and tapers at sharp edges;
  - one 3D CGNS writer for blocks and volumes (structured + unstructured zones, vertex
    connectivity between them).

  - (addendum) quality floor for every 3D cell (scaled Jacobian >= 0.05, skewness <= 0.95 by
    default): enforced in inflation, rescued by the optimiser, reported otherwise (exports
    refuse).

  Criteria: no generated layer cell misses the floor; optimisation never lowers the worst cell
  and reduces cells below the threshold;
  boundary and layer nodes unchanged; fans give every node on interior ridges all layers, a
  closed front and layer volume equal to the volume between wall and front; mixed CGNS passes
  cgnscheck with matching connection coordinates. Demo: wing with fans and optimisation,
  structured block + unstructured duct volume.
- Surface meshing on parametric surfaces, a watertight shell check (each edge shared by exactly 2 faces, consistent orientation).
- `TetBackend` with Netgen first and mmg for optimisation.
- Prism/hex inflation with ridge/corner handling, pyramids at quad faces.
- Element-type counts as expected, no negative volumes, volume vs analytic within 1e-6, inflation height/growth achieved.
- Demos: sphere in a box with prism layers; wing stub hybrid. Export VTKHDF, CGNS (mixed sections + BCs), Gmsh.

## Phase 9 — CAD import, scripting polish, performance *(split: 9a CAD import, 9b scripting + topology helpers, 9c performance + packaging)*

### 9a — STEP / IGES import (OCCT)
- Optional backend `backends/occt` (`OG_WITH_OCCT`, skipped when OpenCASCADE is not found) and
  module `opengrid._occt`; `og_core` does not depend on it.
- STEP and IGES read; healing (ShapeFix), sewing at the larger of the requested and the file's
  tolerances, small-edge removal, splitting of periodic faces and closed edges; each step logged.
- Edges → B-spline curves whose ends are bitwise the vertices; converted curves match OCCT's
  edge lengths to 1e-9 and lie on the part's surfaces to 1e-9. Faces → B-spline surfaces trimmed
  to the face's (u, v) bounds.
- `og.import_cad`: one `nurbs` curve + connector per edge (dimension from spacing and turning
  angle), one domain per face (unstructured, or structured for four-sided faces with shared
  opposite dimensions), one undo step, journal replays without the file. Faces that cannot be
  meshed (poles), free and non-manifold edges are reported.
- Tested: STEP solid, unsewn IGES sewn into a closed shell, sphere faces reported, imported
  solid filled with tetrahedra; GUI action. Demo: filleted bracket with a hole — solid and
  external flow with boundary layers, exported.

### 9b — Scripting polish, topology helpers
- Scripting guide (`docs/SCRIPTING.md`) and an API reference generated from the docstrings
  (`docs/API.md`, `tools/gen_api_docs.py`); a test fails when a public method is undocumented,
  unlisted or the reference is stale.
- The ten golden cases of SPEC.md as standalone scripts (`examples/golden/`, public API only),
  each with invariant checks (exact TFI on the annulus, wall spacings, areas / volumes against
  exact polygons, interfaces, quality floor); a test runs them all and replays their journals.
- Split / join of connectors (nodes kept by explicit distributions), structured domains (cut
  along a grid line, the parts reproduce linear TFI exactly) and six-face blocks (split and
  rejoined blocks have the original points); users (loops, inflation walls, patches, zones,
  volumes) follow. Undo is one step; journals replay.
- O-grid (butterfly) of a structured domain (core + four) and of a six-face block (core + six
  shells, or the shell only): a circle section goes from min scaled Jacobian 0.05 to 0.71, a
  solid sphere from inverted cells to 0.58 — all-hexahedral multi-block meshes.
- GUI: Topology menu. Demo: all-hex pipe and ball (CGNS multi-block, cgnscheck clean), block
  split / join, the golden-case report.

### 9c — Performance and packaging
- Performance targets from SPEC.md measured by a benchmark script in Release on M1
  (`examples/09c_performance.py`): 2D 1M cells < 10 s (structured, Delaunay, advancing front,
  hybrid; mesh + quality + export), 500 x 200 elliptic smoothing with wall control converged in
  < 2 s, 3D 5M tetrahedra via Netgen < 120 s (reported; bounded by Netgen's own time).
- OpenMP for hot loops (`OG_WITH_OPENMP`, `util/parallel.hpp`) with determinism preserved:
  identical results for any thread count (tested). Anderson acceleration of the 2D and 3D
  elliptic smoothers; volume quality metrics in C++; optimiser and Netgen transfer overheads
  removed (identical meshes).
- Packaging: conda package (rattler-build recipe, tested on osx-arm64; the Linux deliverable)
  and a macOS arm64 app bundle + DMG built from it (conda-pack).

## Phase 10 — Multi-block meshing *(proposed 2026-10-15, beyond SPEC.md; approved 2026-10-15, done: 10a-10e)*

Goal: structured multi-block grids smoothed as one grid, edited live with the mouse. The
user draws curves freehand or sketches a rough block topology; block interfaces, corners and boundary nodes move during a single smoothing of the
whole grid; templates and automatic topology cover common cases. Each sub-phase builds on the
one before. Everything stays a Command (undo, journal), deterministic for any thread count, and
the existing patches, exports and quality checks keep working on the result.

### 10a — Global multi-block smoothing, 2D
- A `smoothing group` entity: a set of structured domains smoothed together. Domains in a group
  are evaluated by the group: TFI from their connectors first (the initial guess), then one
  elliptic smoothing of all domains together. The connectors' nodes on interfaces and at
  corners are only that initial guess.
- **Interfaces float:** the domains' nodes are merged into one global node array (connectivity
  from `model/topology.py`). A node on an interface is relaxed with a stencil that continues
  across into the neighbour (ghost rows with the index transform), so grid lines stay smooth
  across the interface instead of turning there.
- **Irregular points:** corners where 3 or 5+ domains meet have no (i, j) stencil. They move by
  a length- and angle-weighted average of their neighbours (Spekreijse-type treatment).
- **Boundary nodes slide** along their database curve (curve projection, warm-started), under
  an orthogonality condition (Neumann slip), or stay fixed. This is chosen per connector, and
  wall control (first-cell height and orthogonality) keeps working on walls.
- **Solver:** the zebra line relaxation and Anderson acceleration from 9c, applied to the
  global state. Results are bit-identical on 1–10 threads.
- **Acceptance:**
  - A domain split in two (9b) and smoothed as a group matches the smoothed single domain to
    the solver tolerance, so the interface leaves no trace.
  - On the O-grid of a circle (core + 4, four valence-3 corners), the largest *kink* of a grid
    line at an interface (its turn there minus the mean turn at its neighbours, ignoring three
    nodes next to irregular points, where lines must turn) drops from the per-domain value
    (34°) to < 0.5°, and the minimum scaled Jacobian does not decrease. (Planned as "turn <
    2°"; the turn includes the curvature of the lines, 2.8° per node on a 128-node circle.)
  - Sliding nodes stay on their curve within 1e-10·L.
  - Exported 1-to-1 interfaces stay valid (cgnscheck clean).
  - Undo is one step; journals replay.
- **GUI:** create a group from the selected domains, set slide/fixed per boundary connector,
  live residual plot.
- **Demo:** a cylinder in a channel (four domains, wall control, sliding channel walls) and
  the circle O-grid, comparing per-domain smoothing with group smoothing (built:
  `examples/10a_multiblock_smoothing.py`; the multi-element airfoil was replaced by these two,
  which show interfaces, irregular points, walls and sliding in a small setup).

### 10b — Global multi-block smoothing, 3D
- Groups of six-face blocks. Interface faces, interface edges and corners float. Irregular
  edges (3 or 5+ blocks around an edge) and irregular vertices get the averaged update from
  10a, extended to 3D.
- **Boundary nodes slide:** boundary faces slide on their surface (surface projection
  warm-started from the previous (u, v)), and boundary edges slide on their curves. Wall
  control works on any boundary face, as in 7b.
- **Acceptance:**
  - A split block smoothed as a group matches the smoothed single block.
  - On the ball O-grid (core + 6 shells), the largest turn of a grid line at an interface is
    < 3° and the min scaled Jacobian is ≥ the 9b value (0.58).
  - Sliding nodes stay on the sphere within 1e-9·R.
  - A 1M-cell group converges in < 10 s on the M1 (Release); results are identical for any
    thread count. *(Built: identical results met; 1M cells take 27 s - missed, see
    OPEN_ISSUES.)*
- **Demo:** sphere in a box and cylinder in a channel, all-hexahedral, CGNS multi-block.
  *(Built: `examples/10b_multiblock_3d.py` - sphere in a box (exported to CGNS) and the solid
  ball; the cylinder in a channel is the 2D case of 10a.)*

### 10c — Live mouse editing and freehand drawing
- **Drag handles** (Select tool; a drag on empty space still pans) for:
  - point entities;
  - the defining points of lines, arcs, polylines and splines;
  - the control points of NURBS curves (control polygon shown; weights in the property
    panel);
  - connector break points.
- **Live updates while dragging:** the dependent connectors, domains, blocks and groups
  regenerate as the handle moves.
  - During the drag: a fast preview (TFI only, or the last converged smoothing as the initial
    guess with a few sweeps), throttled to the frame rate.
  - On release: full regeneration (smoothing to tolerance). One drag is one undo step, and the
    journal records the final position only.
  - Grid nodes themselves are not dragged. They are derived from the geometry (and from 10d
    sketch vertices), so a dragged node would be overwritten by the next regeneration.
- **Freehand drawing:** a Freehand tool records the pointer stroke and fits a cubic NURBS curve
  to it.
  - Fit: least squares with chord-length parameters; knots are added where the error exceeds
    the tolerance (default 0.5 % of the stroke length, set in the tool bar).
  - Corners are detected from the turning angle and kept as C0 breaks (optional).
  - Strokes start and end with snapping to existing points and curve ends.
  - A NURBS tool also creates a curve by clicking its control points.
- **Blocks from curves:**
  - "Domain from curves": select 4 curves (or one closed curve, split at its four sharpest
    corners), and connectors (dimension from the average spacing) plus one structured domain
    are created in one undo step. Three-sided selections are reported, with the option of an
    unstructured domain.
  - "Domains from curve network": strokes that cross are split at their intersections (curve
    intersection by subdivision, then Newton). The faces of the resulting planar network
    become domains: four-sided faces structured, the others reported or unstructured on
    request. Dimensions propagate across shared edges and conflicts are reported. The domains
    can go straight into a 10a smoothing group.
  - In 3D: "Block from curves" builds a six-face block from 12 edge curves (or 4 base curves
    plus an extrusion).
- **Acceptance:**
  - Dragging a spline point updates the dependent structured domain at ≥ 10 frames/s for a
    200×100 grid (preview mode, measured in an offscreen test with simulated drag events).
  - A drag undoes in one step and the journal replays to the same grid.
  - Freehand: synthetic strokes (sampled circle, NACA section, S-curve with noise) are fitted
    within tolerance with ≤ 12 control points per quarter circle, and corners at ≥ 60° turns
    are found.
  - Curve network: a hand-drawn 3×2 lattice of crossing strokes gives 6 structured domains with
    conforming 1-to-1 interfaces (cgnscheck clean after export), and a network containing a
    five-sided face reports it by name.
  - Trackpad and mouse dragging checked by hand *(manual)*.
- **Demo:** an airfoil-like body and its block layout drawn freehand, turned into domains,
  group-smoothed (10a), with the drag recorded as before/after screenshots. *(Built:
  `examples/10c_freehand_live.py` - a channel with a bump drawn freehand; "Block from curves"
  in 3D was not built, see OPEN_ISSUES.)*

### 10d — Sketch topology (floating wireframe)
- A `sketch` entity: a coarse block lattice defined by vertices (points) and blocks (4 vertex
  ids in 2D, 8 in 3D). It is drawn loosely around the geometry, as a wireframe.
- **Dimensions:** set on one edge of each edge family and propagated through the lattice.
  Conflicting dimensions are reported by name.
- **Association with geometry:** sketch edges are assigned to database curves and sketch faces
  to surfaces. An unassigned vertex or edge within a tolerance of geometry is offered for
  assignment (reported, never snapped silently).
- **Generation:** straight-edged TFI from the lattice, assigned nodes projected onto their
  geometry, then the 10a/10b group smoothing with sliding.
- **Output:** ordinary connectors, domains and blocks (explicit distributions) that the sketch
  owns and regenerates. Patches, exports, quality and the topology helpers work on them
  unchanged. Editing them directly is refused, with a message pointing to the sketch.
- **Acceptance:**
  - Moving sketch vertices by up to 20 % of the local block size changes the final grid by
    < 1 % of the local spacing. The result depends on the topology, not on where the sketch
    was drawn.
  - 2D: an airfoil in a C-topology sketch gives a min scaled Jacobian ≥ 0.5. 3D: a sphere in a
    box gives ≥ 0.5.
  - Projected nodes lie on their geometry within 1e-9·L.
  - Undo is one step per edit; journals replay; the `.ogp` round-trips.
- **GUI:** place and drag sketch vertices in the viewport (2D; 3D with axis handles), assign
  geometry from the selection, regenerate with live progress.
- **Demo:** a 2D airfoil and a 3D wing stub, each sketched by hand around CAD geometry imported
  in 9a. *(Built: `examples/10d_sketch_topology.py`, plus the sphere in a box; airfoil min
  scaled Jacobian 0.68, sphere 0.69, 20 % vertex moves change the grid by 2e-4 spacings. The
  wing stub reaches 0.12 in its wake blocks (vertices on 3D curves do not slide, see
  OPEN_ISSUES); 3D vertices have no axis handles.)*

### 10e — Templates and automatic topology
- **Templates:** parametrised generators that create a 10d sketch from the selected geometry:
  - 2D: O-grid and C-grid around a closed body, H-grid channel, body in channel.
  - 3D: pipe / duct (butterfly), sphere in box, cylinder in channel, bump channel, wing in box
    (C-H).
- **Automatic 2D topology:** a cross-field (frame-field) decomposition of a multiply connected
  planar region into four-sided blocks. Separatrices are traced from the field's singularities
  and converted into a sketch, which the user can edit afterwards.
- **No automatic 3D topology:** that is open research. 3D relies on the templates and on
  sweeping or extruding 2D sketches.
- **Acceptance:**
  - Every template meshes its golden geometry with min scaled Jacobian ≥ 0.4 and no
    user-placed vertices.
  - The 2D cross-field topology produces a valid all-quad block decomposition (every block
    four-sided, conforming) for: an airfoil in a box, a two-element airfoil, a plate with three
    holes, and the 9a bracket section.
  - Singularity counts satisfy the index theorem (tested).
- **GUI:** Topology → Templates menu, and an Automatic topology action (2D).
- **Demo:** each template on its geometry, and automatic topology for the two-element airfoil
  and the plate with holes. *(Built: `examples/10e_templates_auto.py`. Nine templates meet
  0.4; wing in box (C-H) reached 0.25 - missed, then 0.62 in the follow-up (extruded 2D
  C-grid; tapered + swept 3D variant 0.60). Automatic topology: valid layouts for the
  airfoil in a box, the plate with three holes and the bracket section; the two-element
  airfoil fails (a five-cornered face in the slot) - missed; the index theorem holds in all
  four. See OPEN_ISSUES.)*

References: J. F. Thompson, Z. U. A. Warsi, C. W. Mastin, Numerical Grid Generation (1985),
ch. VI–VII (multiblock); W. A. Spekreijse, J. Comput. Phys. 118 (1995) 38–61 (irregular points);
P. R. Eiseman, "Multiblock grid generation with automatic zoning", NASA CP-3291 (1995);
N. Kowalski, F. Ledoux, P. Frey, IMR 2012 (cross fields for block decomposition); D. Bommes
et al., Computer Graphics Forum 32 (2013) 51–76 (quad layouts).

## Phase 11 — CFD boundary layers, performance, GPU, high-order meshes *(requested 2026-10-16, beyond SPEC.md)*

Goal: grids for wall-resolved CFD (grid lines normal to walls through the boundary layer at
y+ ~ 1), faster smoothing (CPU first, then the GPU where it pays), and high-order elements for
continuous and discontinuous Galerkin solvers. Same rules as before: Commands, determinism
(per backend for the GPU), validators that report.

### 11a — Wall layers
- Per wall (groups, sketches, templates): `layers` N and `growth`; the first N layers off the
  wall are marched along the front normals (heights spacing * growth^k) and kept; the elliptic
  smoothing controls the first cell beyond them. 2D and 3D, internal walls (wake cuts) and
  symmetry planes included; overlapping or folding layers are reported.
- **Acceptance:** C-grid template at a 1e-5 first cell with 30 layers: mean deviation of the
  grid lines from the wall normal < 0.2° in every layer (> 10° at layer 16 with wall control
  only), first cell exact; the 3D tapered wing likewise; min scaled Jacobian >= 0.4.
  *(Built: `examples/11a_wall_layers.py`: 0.05-0.09° mean in all 30 layers in 2D and 3D, max
  1.1° at the leading edge; min scaled Jacobian 0.47 (2D) and 0.41 (3D).)*

### 11b — CPU performance
- Load-balanced parallel loops (the static schedule left 69 % of the threads idle on the wing
  stub), allocation-free curve and surface evaluation (half of the busy time was malloc/free
  in NURBS evaluation), the remaining serial sections of the group smoothers.
- **Acceptance:** results bit-identical for 1-10 threads; the 10d wing stub group smoothing at
  least 2x faster (39 s); the 9c benchmark not slower. *(Built: `examples/11b_performance.py`:
  wing stub 9.2 s (4.3x), identical grids on 1, 4 and 10 threads; 9c elliptic cases 1.6x
  faster (2D 1.10 -> 0.68 s, 3D block 0.96 -> 0.58 s), the rest unchanged; 1M-cell ball
  26.2 -> 22.3 s.)*

### 11c — Coarse-to-fine group smoothing
- Grid sequencing (smooth a coarsened copy, interpolate, finish on the fine grid) and/or
  multigrid for the group smoothers.
- **Acceptance:** at least 2x fewer fine-grid iterations on the 1M-cell ball (10b) and the
  wing stub, same converged grid (to the tolerance). *(Assessed with a prototype, not built:
  on the 1M-cell ball the fine level needed 218 instead of 328 iterations with the same
  absolute stopping level (1.5x), but the coarse level itself stalled (2000 iterations, 17 s;
  total 30.5 s against 21.0 s), and the grid moved by 1 % of the domain because slide limits
  are measured from the starting positions. A smaller ball (n = 33) stalls at a residual
  floor (relative 1.5e-6 after 2000 iterations) with or without the coarse start. The
  iteration counts are set by that floor, not by smooth error modes.)*

### 11d — GPU assessment (and backend if it pays)
- Precision study: the elliptic update in float32, in double-float (two float32) positions with
  float32 stencil arithmetic on differences, and in float64, on y+ ~ 1 grids. An OpenCL
  prototype of the interior sweep against the CPU kernels after 11b/11c on this M1 Max (no
  float64 on the GPU). Decision from the measurements; if the GPU pays: an optional backend
  with runtime detection (device, float64 support, a self-test against the CPU) and CPU
  fallback, results equal to the CPU within tolerance. *(Assessment done, `tools/gpu/`:
  float64 hardware is not needed - positions as two float32 with float32 arithmetic on
  differences are accurate to 2e-7 of the local cell height at a 1e-6 wall, plain float32 is
  off by 6 %. Kernels on a 1M-node block: GPU 7.4x (line relaxation) and 9x (point Jacobi)
  faster than the same kernels in float64 on the 10 CPU threads; launch overhead cancels the
  gain below ~200k nodes. Decision: build the backend for large 3D groups.)* *(Built
  (`core/src/gpu`, `examples/11d_gpu.py`): detection, self-test, modes, the line sweeps of
  the 3D group smoother on the GPU. Measured on the 1M-cell ball: device sweeps 2.2x faster
  than the CPU's, the whole iteration 1.2x (Anderson and sliding projections stay on the
  CPU), and the GPU colouring needs more iterations - to the CPU's accuracy the GPU run is
  slower. Experimental, off by default.)*

### 11e — High-order elements (CG / DG)
- Curved high-order elements from the meshes: quadrilaterals and hexahedra of order p from the
  structured grids (agglomerating p x p sub-cells of a p-times finer grid, nodes moved to the
  chosen node set, boundary nodes on the CAD geometry), later triangles, tetrahedra and
  prisms (straight-sided elements curved at the boundary and untangled). Node sets:
  equispaced (file formats), Gauss-Lobatto-Legendre (CG, spectral elements; DGSEM), Gauss-
  Legendre (DG solution points); warp-and-blend points for simplices. Validity: the Jacobian
  determinant at the nodes / quadrature points (scaled Jacobian of curved elements). Exports:
  VTK Lagrange cells (VTKHDF), Gmsh high-order .msh, CGNS where its element types allow.
  *(Built for quadrilaterals and hexahedra (`core/src/highorder`, `model/highorder.py`,
  `examples/11e_high_order.py`): p x p (x p) grid cells per element, equispaced / GLL /
  Gauss node sets, GLL nodes snapped onto curves and surfaces, Jacobian validity, VTKHDF
  Lagrange cells (point order checked against VTK) and HDF5 node sets; `interval_multiple`
  on sketches and templates. Boundary error on a circle at the rate h^(p+1) (h^(p+2) for even
  p). Not built: triangles, tetrahedra and prisms, Gmsh and CGNS high-order types.)*

### 11f — Straight-sided high-order elements of every cell type
- Requested by the user (2026-10-16): high-order elements as the cells of the mesh with nodes
  inside (straight edges, planar faces), for triangles, quadrilaterals, tetrahedra, pyramids,
  prisms and hexahedra, with Gmsh and CGNS export. Acceptance: node orders equal to each
  format's own reference (Gmsh, VTK, CGNS) for every kind and order; shared nodes merged
  exactly; exports read back with every element's nodes on its linear map; Gmsh and
  `cgnscheck` accept the files. *(Built (`io/lagrange.py`, `model/elements.py`,
  `io/cgns_elements.py`, `examples/11f_high_order_elements.py`): orders 1-6 (CGNS 1-4),
  Gmsh / CGNS / VTKHDF / HDF5, boundary elements per patch, zones, periodic links. The curved
  structured elements of 11e became `curved_high_order`.)*
