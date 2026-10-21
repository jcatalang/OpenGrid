# OpenGrid — Progress log

## Phase 0 — Kickoff and scaffold (2026-10-11) — approved
- Docs: `docs/ARCHITECTURE.md` (layout, licenses, build, C++/Python boundary, command/undo/journal,
  `.ogp` v1 schema, tet backend decision, spec risks), `docs/ROADMAP.md`.
- Scaffold: CMake + presets + scikit-build-core; `og_core` library; nanobind module `opengrid._core`
  (`version()`, `add()`); PySide6 + pyvistaqt main window (viewport, model tree, property panel,
  log, Python console); GoogleTest + pytest (incl. pytest-qt GUI smoke test); clang-format,
  clang-tidy, ruff configs; `tools/check_no_ternary.py`; GitHub Actions for macos-14 + ubuntu-latest.
- Demo: `examples/00_launch.py` → `examples/output/phase0_main_window.png`.

### Verified on this machine (macOS 15.2, M1, Apple Clang 16, conda-forge env)
- Env: Python 3.12.15, PySide6 6.11.2, VTK 9.7.0, pyvista 0.49.0, nanobind 3.1.0, CMake 4.4.4.
- `ctest`: 3/3 pass (release and `ci` preset with `-Werror`). `pytest`: 12/12 pass, including
  offscreen GUI smoke tests and self-tests of the ternary checker.
- ruff check/format, clang-format, clang-tidy (0 findings in project code), check_no_ternary: clean.
- On-screen launch: the window, viewport, tree, console and log work; HiDPI screenshot at 2800×1800.

### Not verified
- GitHub Actions (no remote configured), the Linux build, the uv-based setup, and trackpad/mouse
  interaction beyond a single launch.

### Notes
- conda 23.5 on this machine only has the classic solver (slow env creation). Consider `conda update -n base conda`.
- conda-forge `clang-tools` lacks clang builtin headers, so clang-tidy needs
  `--extra-arg=-resource-dir=$(clang++ -print-resource-dir)` (see CLAUDE.md).

## Phase 1 — 2D geometry kernel + project I/O (2026-10-11) — approved
- **C++ core (`core/…/geom`, `util`)**: `Curve` base class (eval, d1, d2, breakpoints, arc
  length, arc-length inversion, closest-point projection, adaptive tessellation). Line, Arc,
  Polyline, CubicSpline (natural / clamped / Catmull-Rom, Hermite form on chord-length knots)
  and Nurbs (Piegl & Tiller A2.3/A4.2). Adaptive Gauss–Kronrod 7/15 quadrature.
- **Bindings**: all curve classes, `eval_many`/`tessellate_points` as NumPy arrays, and `project`
  returning `(t, point, distance)`. C++ `invalid_argument` becomes Python `ValueError`.
- **Model and commands (Python)**:
  - `Store` with transactions.
  - `Create`/`Modify`/`Rename`/`Delete`/`CreateEntities` commands, plus `CommandStack` with
    exact snapshot-based undo/redo.
  - A journal of replayable `og.create_*` lines.
  - `opengrid.Session` as the scripting API, and `opengrid.replay`.
  - Validators report entity, field and index (or file and line).
- **Project files**: `.ogp` v1 save/load. Save, load and save again is byte-identical.
- **Import**: Selig and Lednicer airfoil `.dat` (as a single spline, split upper/lower G1
  splines, or a polyline), and CSV/whitespace point lists (with or without a header). Errors
  carry line numbers.
- **GUI**:
  - Tools for point, line, arc (3 points), polyline and spline, with rubber-band and spline
    previews and snapping to points and curve ends (10 px).
  - Picking via core projection, with shift/ctrl for multi-select. Tree ↔ viewport selection sync.
  - Property panel editing every entity type (as undoable commands).
  - Undo/redo with labels, delete, New/Open/Save/Save As, imports, journal export, and a
    modified flag in the title.
  - 2D navigation style (left-drag pan, wheel/scroll zoom).
- **Demo**: `examples/01_geometry.py` builds a NACA 0012 (two clamped splines, G1 at the LE)
  and a farfield circle (R = 20c). It writes `01_geometry.ogp`, `01_geometry_journal.py` and
  `01_geometry.png`, and `--gui` opens the result (`phase1_main_window.png`).

### Verified
- `ctest`: 32/32, including the `ci` preset with `-Werror`. Covered:
  - Arc length within 1e-10 of an independent Richardson chord-sum reference for 11 curves
    (lines, arcs, polylines, three spline types, NACA spline, NURBS circle, polynomial and
    rational NURBS), and exact circle/NURBS-circle lengths (1e-12).
  - Reparametrisation within 1e-10·L and monotone.
  - Projection on curves (< 1e-12·L), off curves along normals, and beyond ends.
  - Finite-difference derivative checks, and tessellation tolerance.
- `pytest`: 99/99. Covered:
  - Undo/redo (two cycles) and journal replay for 15 command variants.
  - Rollback on failure, project round trip and invalid-project reports, importers, and
    property panel parsing/errors.
  - GUI tools, snapping, picking and file workflows (offscreen).
- ruff, ruff format, clang-format, clang-tidy (one NOLINT inside the nanobind macro), and
  check_no_ternary: clean.
- Demo: perimeter 2.0395489 chords (NACA 0012 reference ≈ 2.0396), farfield length exact, and
  the journal replays to an identical project. GUI launched on screen.

### Not verified / known limitations
- Real mouse and trackpad interaction were not exercised by tests (handlers are tested in world
  coordinates). Pinch-to-zoom on the trackpad is not handled. Only one on-screen launch was done.
- The GUI creates geometry in the z = 0 plane only (the API accepts 3D). No periodic/closed
  splines yet. A modify followed by a rename in the property panel takes two undo steps.
- nanobind's exit-time leak report is disabled (Qt/VTK shutdown order). The API-only suite
  runs leak-free with it enabled (`OG_LEAK_WARNINGS=1`).
- CI not run (no remote).

## Phase 2 — Connectors and distributions (2026-10-11) — approved
- **C++**:
  - Distributions: uniform; geometric (begin, end or growth ratio); tanh and sine, one- and
    two-sided. Discrete end-spacing matching, and fallbacks with warnings.
  - Brent root finder, `CurveChain` (arc-length and parametric access), and `connector_nodes`
    (sub-connectors/breakpoints, arc positions, max stretching ratio).
- **Python**:
  - `connector` entity with full validation: curves exist and are curves, the chain is
    continuous, the distribution is consistent, breakpoints increase.
  - Dependency protection on delete/modify, and connectors in `.ogp` and the journal (curve
    names resolved at the time of the command).
  - Spacing propagation at shared nodes, with a per-connector switch.
  - Session API: `create_connector`, `create_connectors_on_curves`, `create_chain_connector`
    (automatic ordering/orientation), `set_dimension`, `set_distribution`, `connector_nodes`
    and `connector_report`.
- **GUI**:
  - Connectors drawn with their nodes and coloured per interval by stretching ratio
    (coolwarm, red at the threshold). Toggles for colouring and nodes; the threshold is set
    from a dialog.
  - Connectors listed in the tree with "⚠" when they have warnings or exceed the threshold.
    Warnings are logged on create/modify.
  - Property panel for dimension, mode, distribution, breakpoints (one line each) and
    propagation, plus live length, achieved spacings, max ratio, propagation source and
    warnings.
  - "Connectors on curves" and "Chain connector" actions with a default-dimension box. Picking
    prefers a connector over its curve.
- **Demo**: `examples/02_connectors.py` sets up a NACA 0012 with LE 1e-3 / TE 5e-4 tanh
  clustering (121 points per surface), a wake whose TE spacing is propagated, and a farfield
  with a breakpoint. Outputs: `02_connectors.ogp`, journal, `02_connectors.png` and
  `phase2_main_window.png`.

### Verified
- `ctest`: 50/50, including the `ci` preset with `-Werror`.
  - Achieved end spacings within 1e-8 relative over sweeps of type × points (3 to 1001) ×
    spacing (begin down to 1e-12, far end down to 1e-7).
  - Monotonicity, constant geometric ratio, deterministic output, and documented fallbacks for
    8 infeasible cases with warnings.
  - Tiny spacings down to 1e-300 stay finite, and 3000 seeded random requests produce no
    NaN or non-monotone result.
  - Chain continuity and gap errors, breakpoint nodes and spacings, and parametric vs
    arc-length mode.
- `pytest`: 136/136. Covered:
  - Connector API, validation (14 error cases), propagation and its undo.
  - Dependency protection, renames, project round trip with identical nodes, and journal
    replay.
  - GUI creation, picking, editing, threshold warnings and toggles.
- Demo: every connector end spacing within 1e-13 relative (requirement 1e-8). The breakpoint
  node sits at 0.5 of the circumference to 12 digits, and the journal replay gives identical
  nodes.
- ruff, ruff format, clang-format, clang-tidy (0 findings) and check_no_ternary: clean.

### Not verified / known limitations
- The far-end spacing below about 1e-7·L can't reach 1e-8 relative accuracy (a
  floating-point limit; see ARCHITECTURE).
- Spacing propagation is one level deep and only fills tanh/sine ends that have no spacing.
  Geometric and uniform connectors never adopt spacings.
- Connector split/join and dimensioning by average spacing are not in Phase 2 (the spec lists
  them among the GUI edit tools); they're candidates for Phase 3 prep.
- Interaction was not tested with a real mouse or trackpad. CI not run (no remote).

## Phase 3 — Structured 2D domains + first export (2026-10-11) — approved
- **C++**:
  - `tfi_linear` (Soni arc-length blending) and `tfi_hermite` (wall-orthogonal).
  - `smooth_elliptic`: Laplace / Thomas–Middlecoff / TTM controls, Steger–Sorenson wall
    control with Hilgenstock feedback, point SOR or alternating line SOR, progress/cancel from
    another thread.
  - `quality_2d`: scaled Jacobian, skewness, aspect ratio, min/max angle, FV
    non-orthogonality, neighbour area ratio, area.
- **Python**:
  - `domain` entity (explicit edges or automatic assembly: simple loop or O-type) with full
    validation and dependency protection (curve → connector → domain, transitive).
  - Canonical shared-node coordinates; generator settings; cached smoothing.
  - Session API: `create_structured_domain`, `create_domain`, `set_domain_generator`,
    `rotate_domain`, `domain_grid`, `smooth_domain`, `domain_quality`, `domain_status`,
    `export_vtkhdf`, `export_plot3d`.
  - Writers: VTKHDF (merged conforming quads, block id, cell indices, all quality metrics,
    block names) and Plot3D (multi/single grid, 2D or 3D with nk = 1, ASCII or Fortran binary)
    with a reader.
- **GUI**:
  - Domains drawn as grid lines or coloured by a chosen quality metric.
  - Quality dock (histogram, statistics, show worst cell) and Smoother dock (live residual
    plot, run/cancel on a worker thread).
  - Domain editor (TFI type, Hermite walls, elliptic settings, wall spacing, TTM attractions,
    rotate edges, run smoother).
  - "Structured domain" from selected connectors, and File → Export VTKHDF / Plot3D. Domains
    appear in the tree with "⚠" while smoothing is pending or not converged.
- **Demos**:
  - `examples/03_airfoil_grids.py`: NACA 0012 O-grid (201×81) and C-grid (321×81), wall
    spacing 1e-4 c.
  - `examples/03_backward_step.py`: three-block backward-facing step with shared connectors.
  - Outputs: `.ogp`, journals, `.png`, `phase3_main_window.png`, and `.vtkhdf`/`.xyz`
    (git-ignored, regenerated by the demos).

### Verified
- `ctest`: 63/63 (including `-Werror`). Covered:
  - TFI exact on a quarter annulus with non-uniform radial and angular distributions
    (< 3e-12), boundaries reproduced exactly, and Hermite with no walls equal to linear.
  - Hermite orthogonal at a curved wall (< 1°, against > 10° for linear TFI).
  - Laplace leaves a uniform rectangle unchanged; Thomas–Middlecoff leaves a separable
    clustered grid unchanged (< 1e-9).
  - Laplace on the annulus converges by 8 orders with both solvers.
  - Wall control on a bump channel: spacing within 5% and angle within 2° (achieved: exact to
    4 digits and 0.000°).
  - Cancel works, also from another thread.
  - Quality metrics on known cells (squares, a 30° shear, an inverted cell, growth ratio).
- `pytest`: 168/168. Covered:
  - Domain assembly, validation, dependencies, rotation, bitwise shared-node matching.
  - Hermite and wall-controlled smoothing through the API, cache/undo/determinism, cancel.
  - Round trip and journal replay with domains.
  - VTKHDF read back with VTK's `vtkHDFReader` (point merge, conforming edge count, cell data,
    field data, layout, determinism), and Plot3D round trips (ASCII/binary × 2D/3D, single
    grid).
  - GUI: creation, quality dock, worst cell, generator editing, threaded smoothing and cancel,
    rotate button, exports.
- Demos:
  - O-grid: converged in 891 iterations / 0.8 s, wall spacing 1.0000e-4 everywhere,
    orthogonal to 0.00°, min scaled Jacobian 0.14 (in the sharp-TE corner cell).
  - C-grid: converged in 2110 iterations / 3.2 s, spacing exact, 0.03° from orthogonal, min
    scaled Jacobian 0.53.
  - Step: area exactly 45 (1e-12), all Jacobians 1, the 212 shared points merged.
  - All journals replay to identical projects.
- ruff, ruff format, clang-format, clang-tidy (0 findings) and check_no_ternary: clean.

### Not verified / known limitations
- **meshio cannot read VTKHDF or Plot3D (5.3.5)**, so the spec's meshio read-back check isn't
  possible. VTK's reader is used instead.
- **"Cancel returns within 0.5 s"** was not timed. Cancellation is checked once per iteration,
  which takes about 1 ms on the demo grids.
- **Hermite TFI can fold C-type grids** (documented); the C-grid demo uses plain TFI.
- **Elliptic smoothing** is for planar domains in z = const. The grid cut is treated as a fixed
  boundary (no periodic smoothing across it). Wall control is available on j = 0 / j = jmax
  only (rotate edges for other walls).
- **Automatic assembly** doesn't handle C-type loops (a connector used on both sides of a cut);
  use explicit edges.
- **The GUI** shows overlapping domains on top of each other (the demo project holds both
  airfoil grids at the same place). Real mouse/trackpad interaction is untested. CI not run.
- Phase 1/2 bugs found and fixed while doing this phase:
  - With connectors present, snapping in the creation tools crashed (regression test added).
  - Line SOR was unstable with explicit sources.

## Phase 4 — Patches, zones, solver-ready export (2026-10-11) — approved
- **Model**:
  - `patch` and `zone` entities: validation, dependency protection, journal and `.ogp`
    sections.
  - Topology analysis and boundary validator: unassigned or multiply-assigned boundary
    connectors, interfaces placed in patches, non-manifold connectors, domains in several
    zones.
  - Interfaces and cuts are detected automatically; 1-to-1 ranges and transforms are computed.
- **API**: `create_patch`, `set_patch_type`, `add_to_patch`, `remove_from_patch`,
  `create_zone`, `add_to_zone`, `remove_from_zone`, `boundary_report`,
  `unassigned_boundary_connectors`, `export_cgns`, `export_gmsh`, and `export_vtkhdf` (now a
  PartitionedDataSetCollection with patch blocks).
- **C++**: CGNS structured writer (optional build), checked with `cgnscheck`.
- **GUI**:
  - Patches and zones in the tree.
  - Property-panel editing of patch type and members, and of zone members.
  - Boundary dock: create patches/zones from the selection, add/remove, validate, with a
    summary plus an issue and interface list.
  - "Colour connectors by patch" (patch palette, unassigned in red, interfaces in grey; domain
    wireframes faded).
  - File → Export CGNS / Gmsh.
- **Demo**: `examples/04_solver_ready.py`, run on the backward-facing step (3 blocks, 2
  interfaces, 2 zones) and the NACA 0012 C-grid (wake cut). Outputs: CGNS, Gmsh (git-ignored),
  VTKHDF (git-ignored), `.ogp`, journals, and `phase4_main_window.png`.

### Verified
- `ctest`: 63/63 (`-Werror`); clang-tidy: 0 findings.
- `pytest`: 182/182. Covered:
  - Boundary report and interfaces (step, O-grid cut, C-grid cut), validator messages, the
    non-manifold case.
  - Entity dependency/rename/undo, round trip and replay.
  - Gmsh read by meshio: correct physical group names, tags and dimensions; patch face counts
    and zone cell counts exact; conforming node count.
  - VTKHDF read by `vtkHDFReader`: block names, assembly, per-block cell counts, zone ids.
  - CGNS: `cgnscheck` with no warnings or errors; zones, families and BC counts read back with
    h5py; transforms [1, 2] (step, O-cut) and [−1, −2] (C-wake); PointRange/PointRangeDonor
    end points coincide in the grid; 32-character name limit.
  - GUI patch workflow, invalid export reported, patch colouring, property-panel editing.
- Demo: for both cases the patch face counts agree between the model, Gmsh/meshio and
  VTKHDF/VTK, cgnscheck is clean, and the journals replay.

### Not verified / known limitations
- **meshio cannot read CGNS here** (no reader in 5.3.5). CGNS is checked with cgnscheck and
  h5py instead.
- **No solver was run** on the exported files.
- **Patches hold whole connectors.** A connector shared by a wall and a symmetry boundary has
  to be split (connector splitting is not implemented yet).
- **CGNS 2D** is written with PhysicalDimension 2 (x, y only). Some solvers expect 3D data with
  one cell in z (the Phase 7 extrusion will provide that).
- **Gmsh** classifies all nodes in one surface entity block; meshio reads it, but Gmsh itself
  was not run (not installed).
- CI not run (environment.yml now also lists `zlib`).

### Phase 4 addendum (after review): periodic pairing and Gmsh check
- Periodic patches are paired from one side (`og.pair_periodic(source, target, translation=None |
  vector, rotation_angle, rotation_center)`; automatic centroid translation when none given).
  Validator: each periodic patch paired exactly once, node-by-node connector matching (either
  direction), mismatch distances reported.
- Exports: CGNS periodic GridConnectivity1to1 (both directions, inverse translation, exponents;
  cgnscheck clean; rotational periodicity rejected in 2D); Gmsh `$Periodic` (affine + node
  correspondences, verified via meshio); VTKHDF field data `periodic_partner` / `periodic_affine` /
  `periodic_source`. Project load handles pairings in a second pass.
- Gmsh 4.12 (/Applications/Gmsh.app, or $OPENGRID_GMSH) reads all exported meshes without errors or
  warnings, with matching node/element counts (tests skip if Gmsh is absent).

## Phase 5a — Unstructured 2D: constrained Delaunay + refinement (2026-10-11) — approved, committed 5f2838c
- **C++**: exact predicates, Bowyer–Watson + Sloan recovery + Lawson flips + classification,
  Delaunay refinement (circumcenters / off-centers, boundary never split, diametral-lens
  rejection), size field (boundary spacing + growth, max size, point/line sources), triangle
  quality metrics.
- **Python**:
  - `unstructured_domain` entity (validation of loops/orientation/holes, automatic loop
    assembly, generator settings, cached meshes).
  - API: `create_unstructured_domain`, `set_unstructured_generator`, `unstructured_mesh`;
    `domain_quality` and `domain_status` cover both kinds.
  - Topology, zones and patches cover both kinds. Exports handle mixed cells (Gmsh, VTKHDF).
- **GUI**: unstructured domains drawn (wireframe or quality colouring), picking, tree count,
  generator editor (min angle, growth, max size, max points, sources), "Unstructured domain"
  from the selected connectors, Quality dock and worst cell for triangles.
- **Demo** `examples/05a_unstructured.py`:
  - NACA 0012 in a farfield: 9866 triangles, min angle 28.15°.
  - Hybrid structured O-ring (wall 1e-4 c, elliptic) + triangles (7806, min angle 28.61°),
    conforming along the shared ellipse.
  - Gmsh files read clean by Gmsh 4.12. Outputs: `05a_unstructured.png`,
    `phase5a_main_window.png`, projects, journals.

### Verified
- `ctest`: 73/73 (`-Werror`); clang-tidy: 0 findings. Covered:
  - orient2d exact on 20 000 near-collinear lattice cases (checked against 128-bit integers);
    incircle exact on cocircular points.
  - Size field equal to brute force; conformity (every segment exactly once, on its left).
  - Exact area, no inverted triangles, ≥ 20.7° on uniform and graded domains, holes,
    determinism, invalid input, the abrupt-spacing case terminating with < 5 % bad triangles.
- `pytest`: 202/202. Covered:
  - Annulus area to 1e-12 vs the polygon, ≥ 20.7°, conformity, boundary points identical to
    the connector nodes.
  - Determinism and caching, size controls, loop validation and assembly, dependencies,
    round trip, replay.
  - Hybrid conformity: merged node count, every interior edge shared twice, Gmsh read clean,
    VTKHDF mixed cells.
  - CGNS and Plot3D refusals; GUI create, pick, edit, quality.
- Performance: 1.39 M triangles in 1.0 s (target 1 M in < 10 s).

### Not verified / known limitations
- **Quality next to abrupt boundary-spacing jumps** is not guaranteed (boundary nodes are
  fixed); quality checks run away from such jumps.
- **CGNS unstructured export**: Phase 5c.
- **Not in 5a:** advancing front (5b); quads and smoothing (5c).

## Phase 5b — Advancing-front triangulation (2026-10-11) — approved
- **C++:** `advancing_front()` and `flip_to_delaunay()` (neighbour-array Lawson flips),
  `smooth_interior()` (guarded Laplacian).
- **Size field:** a k-d tree replaces the bucket grid.
- **Predicates:** allocation-free exact `orient2d` path.
- **Python / GUI:**
  - `generator.algorithm` ("delaunay" | "advancing_front"), validated, saved, journaled, undoable.
  - Older projects default to Delaunay.
  - GUI algorithm combo; min angle and refine are disabled for the advancing front.
- **Demo:** `examples/05b_advancing_front.py`, NACA 0012 in R = 20 c, growth 1.15.

  | | Triangles | Min angle | Mean min angle | Under 30° |
  |---|---|---|---|---|
  | Delaunay | 9866 | 28.15° | 46.8° | 0.21% |
  | Advancing front | 7490 | 29.15° | 52.5° | 0.11% |

  - Hybrid O-ring + advancing front outside: 5958 triangles, min angle 31.31°, conforming.
  - Gmsh 4.12 reads both files clean.
  - Outputs: `05b_advancing_front.png`, `phase5b_main_window.png`, projects, journals.

### Verified
- `ctest`: 80/80 (`-Werror`); clang-tidy: 0 findings. Advancing-front tests:
  - square with hole: exact area, ≥ 25°, mean ≥ 45°, no warnings;
  - graded L-shape: ≥ 20.7°, deterministic;
  - abrupt spacing jump: valid mesh;
  - invalid input rejected;
  - uniform disc: median size ratio 1 ± 0.05, ≥ 30°;
  - flips: fan → locally Delaunay, area and orientation preserved, idempotent.
- `pytest`: 206/206 (+4: advancing-front annulus, algorithm switch/undo/save/old projects,
  sources, GUI combo).
- Benchmarks (scratch probes, not in CI), advancing front vs Delaunay:

  | Case | Triangles | Median size ratio | Min angle | Mean min angle |
  |---|---|---|---|---|
  | Graded cylinder | 6770 / 9058 | 1.01 / 0.91 | 32.8° / 28.1° | 51.4° / 47.0° |
  | 1.1 M-triangle disc (2.4 s vs Delaunay 1.1 s) | — | 1.00 | 32.8° | 59.4° |

### Not verified / known limitations
- **No angle guarantee:** the advancing front's angles are measured, not proven. The
  Delaunay-fallback path for leftover edges is covered by code but never triggered in the tests.
- **Speed:** the advancing front is about 2× slower than Delaunay at 1 M triangles.
- **Status label:** clipped in the GUI property panel (pre-existing).

## Phase 5c — Quad-dominant / all-quad, smoothing, edge swapping, CGNS unstructured (2026-10-11) — approved
- **C++:**
  - `quality_cells` (mixed metrics); `recombine` (greedy + augmenting paths);
    `smoothing` (Laplacian / length-weighted / angle-based, guarded, edge swapping, untangling);
  - `domain_mesher` (pipeline; all-quad via coarse boundary + subdivision);
  - `write_cgns_2d` (unstructured zones, BAR_2 BCs, vertex connectivity, periodic).
- **Python / GUI:**
  - New generator keys `element_type`, `smoothing`, `smoothing_iterations`, `edge_swap`,
    validated, saved and journaled; older projects keep their meshes.
  - Mesh dict gains `quads`.
  - Quality, picking, worst cell, model tree and status handle quads.
  - Gmsh/VTKHDF write mixed cells.
  - CGNS writes hybrid and unstructured meshes; Plot3D still refuses unstructured domains.
  - GUI editor: element type, smoothing method, iterations, edge swapping.
- **Demo** `examples/05c_quads.py`: NACA 0012 in R = 20 c (advancing front + angle-based
  smoothing).

  | Element type | Cells | Min SJ | Mean SJ |
  |---|---|---|---|
  | Triangles | 7490 | 0.624 | 0.913 |
  | Quad-dominant | 3641 quads + 208 triangles | 0.668 | 0.888 |
  | All-quad | 10760 quads | 0.683 | 0.920 |

  - Hybrid O-ring + all-quad outside: 7208 quads, min SJ 0.665.
  - cgnscheck, Gmsh 4.12 and meshio read everything clean.
  - Outputs: `05c_quads.png`, `phase5c_main_window.png`, CGNS / msh / VTKHDF, projects,
    journals.

### Verified
- `ctest`: 88/88 (`-Werror`); clang-tidy: 0 findings. Covered:
  - right-triangle pairs → exact squares;
  - quad-dominant annulus: conforming, exact area, ≥ 90% of the area in convex quads;
  - all-quad annulus with both triangulators: no triangles, conforming, exact area, min SJ > 0.3;
  - odd loop rejected; determinism;
  - triangles without smoothing = Phase 5a output;
  - all three smoothing methods: worst kept, mean improved, boundary fixed;
  - untangling repairs a node thrown across its ring.
- `pytest`: 215/215 (+9). Covered:
  - quad-dominant / all-quad annulus: area, SJ > 0, conformity, boundary nodes;
  - even-loop check; generator validation, save/load, replay; smoothing;
  - CGNS hybrid for all three element types (structure, coordinates, BAR_2 ranges,
    EdgeCenter, interface point lists matching coordinates, cgnscheck);
  - periodic unstructured channel (triangles and quads) in CGNS (translated point lists) and
    Gmsh.
- 1 M-cell benchmark: every element type in 1.2–3.4 s.

### Not verified / known limitations
- **All-quad parity:** needs an even number of segments on each loop. Two odd loops (even
  total) are not handled; loops touching at a point are rejected.
- **Quad-dominant:** is not guaranteed to be all-quad (2–8% triangles).
- **No angle guarantee for quads:** quality is measured (min SJ 0.45–0.68 in tests and demos),
  not proven.
- **Plot3D:** still structured only.
- **GUI:** the property panel's status label clipping (pre-existing) remains.

## Phase 6 — 2D inflation layers (2026-10-11) — approved
- **C++:**
  - `inflation.cpp`: advancing layers with smoothed normals, fans, merges, transition, slope,
    crossing and proximity termination, cell checks, and corner-aware parity for all-quad
    cores.
  - `domain_mesher` runs the inflation, then meshes the core.
  - `loop_corners.hpp`; the untangler gains a quality threshold (0.3).
- **Python:**
  - `generator.inflation = {"connectors", "first_height", "growth_rate", "layers",
    "total_thickness", "transition"}`, validated (connectors must be in the loops).
  - Stored as ids; the journal writes names; save/load and undo work.
  - The mesh reports `layers`, and the status shows the layer count.
- **GUI:** inflation walls, first height, layer growth, max layers, total thickness and a
  "stop at core size" checkbox in the unstructured editor.
- **Demo** `examples/06_inflation.py` (first height 1e-5 c, growth 1.2, at most 60 layers;
  the transition stops at 43):

  | Case | Cells | Min SJ | Mean SJ |
  |---|---|---|---|
  | NACA 0012, triangle core | 6748 triangles + 7735 quads | 0.360 (trailing-edge fan triangles) | 0.957 |
  | NACA 0012, all-quad core | 16971 quads + 78 layer triangles | 0.357 | 0.959 |
  | Two-element (30% flap, 20°, small gap) | 8492 triangles + 11905 quads | 0.360 | 0.959 |

  - First-layer height error is at most 0.25%.
  - Layers stop short of each other in the gap.
  - cgnscheck and Gmsh 4.12 read everything clean.
  - Outputs: `06_inflation.png`, `phase6_main_window.png`, CGNS / msh / VTKHDF, projects,
    journals.

### Verified
- `ctest`: 95/95 (`-Werror`); clang-tidy: 0 findings. Inflation tests:
  - straight wall: every layer exactly h₀gᵏ (1e-6);
  - cylinder with triangle and all-quad cores: first height within 2%, transition before the
    layer limit, valid and conforming, exact area;
  - sharp and blunt trailing edges: fans with layers behind the edge, valid;
  - two-element gap: valid, core cells in the gap;
  - last layer within 20% of the core size;
  - invalid options rejected; determinism.
- `pytest`: 220/220 (+5). Covered:
  - first-layer height; validation (missing first height, connector not in the loops,
    layers < 1);
  - save/load/replay, rename, undo;
  - Gmsh conformity (only wall and farfield edges unmatched), cgnscheck for triangle and
    all-quad cores;
  - GUI fields.
- All Phase 1–5 demos still run. Phase 5c all-quad meshes changed slightly (see ARCHITECTURE).

### Not verified / known limitations
- **Open wall ends:** layers taper to zero height where a wall meets a non-wall connector
  (fixed connector nodes).
- **Layer smoothing:** layer cells are not smoothed.
- **First-layer fan triangles** at sharp trailing edges are the worst cells (SJ ≈ 0.36).
- **Python options:** `max_slope`, the fan angles and the proximity factor are not exposed.
- **3D:** no inflation yet (Phase 8).

## Phase 7a — 3D structured foundation (2026-10-11) — approved
- **C++:**
  - `geom/surface.cpp` (5 surface types, projection);
  - `structured/tfi3d.cpp`, `extrude.cpp`, `surface_tfi.cpp`;
  - `quality/quality3d.cpp`;
  - `write_cgns_structured_3d`;
  - bindings in `bindings/og_3d.cpp`.
- **Python:**
  - Surface entities and API (`create_plane`, `create_ruled_surface`,
    `create_extruded_surface`, `create_revolved_surface`, `create_nurbs_surface`,
    `project_to_surface`).
  - Domains with `surface` or in 3D.
  - `block` entity (`create_block(faces=...)`, `create_extruded_block`, `block_grid`,
    `block_quality`, `block_status`) and the 3D topology.
  - Face/block-face patch members; zones with blocks.
  - 3D exports (CGNS, VTKHDF, Gmsh, Plot3D); `blocks` project section.
- **GUI:** blocks and surfaces in the viewport, block actions, extrusion editor, cut plane,
  shrink, 3D quality and worst cell.
- **Demo** `examples/07a_blocks.py`:
  - Annular sector duct of two six-face blocks (21 × 41 × 21 each):
    - min scaled Jacobian 1.000;
    - interface `duct-a.jmax | duct-b.jmin`;
    - volume −2.9e-5 relative to analytic (chord error at 41 points per 30°).
  - Extruded wing section (O-ring around a NACA 0012, 201 × 41 × 21, 160 000 hexahedra):
    - wake cut detected as a block-internal interface;
    - min scaled Jacobian 0.144 (2D ring's trailing edge).
  - cgnscheck clean, Gmsh 4.12 reads clean, journals replay.
  - Outputs: `07a_blocks.png`, `phase7a_main_window.png`, CGNS / VTKHDF / msh / xyz, projects,
    journals.

### Verified
- `ctest`: 101/101 (`-Werror`); clang-tidy: 0 findings. Covered:
  - projection of on-surface points within 1e-10·size for all five surfaces; exact
    off-surface distances; analytic vs finite-difference derivatives; NURBS exact cylinder;
  - block TFI reproduces an annular sector to 1e-13, volume within 1e-6 (1201 θ points);
    mismatched faces rejected;
  - unit-cube quality; translate/rotate/normal volumes, mirroring;
  - surface TFI stays on the cylinder.
- `pytest`: 228/228 (+8). Covered:
  - surfaces (projection, dependency, rebuild after a curve change, replay);
  - six-face block (volume within 1e-6 at 601 θ points, any face order, radii exact);
  - validation (not six faces, both or neither of faces/extrusion, duplicates, faces not
    forming a box, Hermite on a 3D face);
  - two-block interface: every CGNS 1-to-1 point mapped through the transform onto the donor
    point, cgnscheck;
  - extrusions (exact volumes, rotate volume, base interface between up and down extrusions,
    replay);
  - 3D patches and exports (Gmsh element counts, Gmsh reads clean, VTKHDF hexahedra,
    Plot3D 3D round trip, unassigned and connector members reported);
  - project round trip; GUI (block from selection, extrusion editing, cut, shrink, quality,
    worst cell).
- All earlier demos still run.

### Not verified / known limitations
- **Not in 7a (moving to 7b):** surface and 3D elliptic smoothing, hyperbolic extrusion, the
  3D cylinder O-grid demo.
- **Block faces:** periodic 3D faces and collapsed faces are not supported.
- **Domain quality in 3D:** non-planar face domains use planar 2D metrics.
- **Cut plane:** fixed normal (x).

## Phase 7b — 3D elliptic smoothing, surface smoothing, hyperbolic extrusion (2026-10-11) — approved
- **C++:** `elliptic3d.cpp`, `surface_elliptic.cpp`, `hyperbolic.cpp`, with bindings
  `smooth_elliptic_3d`, `smooth_on_surface`, `extrude_hyperbolic`.
- **Python / GUI:**
  - Block `smoothing` setting (`set_block_smoothing`) and extrusion mode `hyperbolic` (edge
    conditions, dissipation).
  - Elliptic generator allowed on surface domains.
  - Block status reports smoothing convergence and warnings.
  - GUI editors for block smoothing and the new extrusion fields.
- **Demo** `examples/07b_smoothing.py`:
  - **3D cylinder O-grid** (hyperbolic, 129 × 21 × 61, 153 600 cells in 0.08 s):
    - first cell 1.000e-3, outer radius exactly 10.5;
    - min scaled Jacobian 0.9997, max non-orthogonality 1.41°;
    - seam detected as a 1-to-1 interface.
  - **Bump channel**: 3D elliptic smoothing with wall control brings the first-spacing error
    from 0.150 to 0.000 and the wall angle from 19.5° to 0.00°.
  - cgnscheck clean. Outputs: `07b_smoothing.png`, `phase7b_main_window.png`, CGNS / VTKHDF /
    msh, projects, journals.

### Verified
- `ctest`: 107/107 (`-Werror`); clang-tidy: 0 findings. Covered:
  - Laplace keeps a uniform cube (1e-12);
  - Thomas–Middlecoff converges to the same grid from an exact and a perturbed start (1e-8);
  - wall control: spacing within 2%, wall angle < 1°, faces untouched;
  - surface smoothing stays on the cylinder (1e-11) and recovers a perturbed patch;
    walls rejected;
  - hyperbolic: flat-plate heights exact (1e-12); cylinder O-grid radial (1e-9); mirroring;
    open periodic edges rejected.
- `pytest`: 232/232 (+4). Covered:
  - hyperbolic O-grid through the API (radii, first cell, interface, cgnscheck, replay);
  - block smoothing with wall control (spacing < 2%, angle < 1.5°, faces unchanged, undo,
    validation);
  - surface domain smoothing on a cylinder; 3D validation; save/load.
- All Phase 1–7a demos still run.

### Not verified / known limitations
- **Hyperbolic extrusion** is an explicit variant, not the implicit Chan–Steger scheme.
  Concave regions rely on dissipation; inverted layers are reported, not repaired.
- **3D elliptic smoothing** uses point SOR only (no line relaxation): about 850 iterations for
  18 000 cells.
- **Surface smoothing** supports no wall control.

## Phase 8a — Unstructured 3D volumes (2026-10-11) — approved
- **Environment:** `netgen` (conda-forge, LGPL-2.1) added to `environment.yml` and installed;
  VTK 9.7, PySide6 6.11, pyvista 0.49 and NumPy 2.5.4 are unchanged.
- **Python:**
  - Frames for unstructured domains in 3D planes and on surfaces (`surface` key, API
    `create_unstructured_domain(..., surface=)`); 3D face quality.
  - `volume` entity (`create_volume`, `set_volume_generator`, `volume_mesh`,
    `volume_quality`, `volume_status`): shell assembly, watertightness and orientation, inner
    bodies, pyramid modes, Netgen backend, vectorised tet/pyramid quality.
  - Combined 3D topology for blocks and volumes. Exports: VTKHDF/Gmsh mixed cells; CGNS
    unstructured 3D zones (new C++ writer and binding).
- **GUI:** volumes in the tree and viewport (cut, shrink, quality, worst cell), "Volume from
  faces", volume generator editor.
- **Demo** `examples/08a_volumes.py`: finite cylinder in a channel; the sides lie on surfaces.
  - "near" refinement volume: 117 689 tetrahedra, min dihedral 18.1°.
  - "far" volume: 41 574 tetrahedra, min dihedral 12.8°.
  - Conforming across the 6 refinement-box faces; volume within 1.3e-5 of analytic (cylinder
    as a polygon).
  - cgnscheck clean. Outputs: `08a_volumes.png`, `phase8a_main_window.png`, CGNS / VTKHDF /
    msh, project, journal.

### Verified
- `ctest`: 107/107; clang-tidy: 0 findings.
- `pytest`: 239/239 (+7). Covered:
  - 3D domains: plane area exact, cylinder patch on the surface, non-planar rejected, replay;
  - volume with an inner body: volume exact (35), positive cells, boundary nodes bitwise,
    tetrahedron faces paired twice except exactly the boundary faces, report OK;
  - shell errors; pyramids always/never with exact volume;
  - two volumes sharing a face: identical interface nodes, CGNS connectivity, cgnscheck;
  - Gmsh/VTKHDF counts; GUI.
- All Phase 1–7 demos still run.

### Not verified / known limitations
- **Pyramids at sharp edges:** pyramids on boundary quads meeting at sharp edges give thin
  tetrahedra ("auto" avoids pyramids except where blocks share a face; 8c).
- **Surface frames:** a constant metric correction; strongly varying metrics are approximated.
- **Not in 8a:** CGNS for mixed blocks + volumes and mmg optimisation (8c); 3D inflation (8b).
- **Performance:** "5 M tets < 2 min" not measured (the demo's 118 k tetrahedra including
  quality take about 1 s; most of the demo's 11 s is surface meshing and assembly in Python).

## Phase 8b — 3D inflation layers (2026-10-12) — approved
- **C++** `inflate_3d` (`unstructured/inflation3d.cpp`, binding `_core.inflate_3d`):
  - prism/hexahedron layers from wall faces of a closed shell;
  - normal smoothing and visibility, scaled normals;
  - stops by layers, thickness, transition, slope, proximity, collision (BVH), visibility and
    quality;
  - transition cells (pyramids, tetrahedra, centroid splits);
  - sliding nodes in planes; side faces and per-node results.
- **Python:**
  - Volume generator key `inflation` = `{faces, first_height, growth_rate, layers,
    total_thickness, transition, slide}`, validated, with names in journals.
  - Sliding planar faces are re-meshed around the layers.
  - Front quads get pyramids.
  - Cell types tets/pyramids/prisms/hexes everywhere (`CELL_TYPES`, `volume_cells`): quality
    (dihedral angles for all types, scaled Jacobian capped at 1), divergence-theorem volumes,
    VTKHDF/Gmsh (wedge 6)/CGNS (PENTA_6, HEXA_8), scene, picking, tree.
- **GUI:** inflation rows in the volume editor (walls, first height, growth, layers,
  thickness, stop at core size, slide).
- **Demo** `examples/08b_inflation.py`:
  - **Sphere in a box** (cubed-sphere unstructured patches on three revolved surfaces):
    - 37 200 prisms in 15 layers on every wall node, 21 714 tetrahedra;
    - min scaled Jacobian 0.51 (prisms) / 0.14 (tetrahedra);
    - volume equal to the enclosed volume (faceted sphere +2.0e-5 vs analytic).
  - **NACA 0012 wing stub** (structured wing surfaces, unstructured tip, sliding symmetry
    plane):
    - 28 183 hexahedra, 2 347 prisms, 4 009 pyramids, 46 242 tetrahedra, up to 20 layers;
    - the 24 trailing-edge nodes stop and the layers taper towards them;
    - volume equal to the enclosed volume.
  - cgnscheck is clean for both. Outputs: `08b_inflation.png`, `phase8b_main_window.png`,
    CGNS / VTKHDF / msh, projects, journals.

### Verified
- `ctest`: 113/113 (+6 `Inflation3D`). Covered:
  - exact first height normal to flat walls;
  - exact total height at an untouched node;
  - prism/hex counts;
  - fixed-node taper;
  - sliding nodes stay in their plane;
  - proximity stops;
  - invalid input.
- clang-tidy clean.
- `pytest`: 246/246 (+7). Covered:
  - prism counts and achieved heights;
  - hexahedra with pyramid caps;
  - sliding floor: area of the re-meshed floor, VTK cell volumes positive and summing to the
    enclosed volume, Gmsh counts, cgnscheck;
  - journal and project round trips;
  - the slide option off;
  - validation, total thickness, GUI editor.
- All demos 01–08b run.

### Not verified / known limitations
- **Sharp ridges:** no 3D fans. Columns on sharp convex ridges (e.g. trailing edges) stop and the
  layers taper over a large part of the wall (slope 0.5).
- **Thin transition cells:** in anisotropic surface cells (e.g. 0.006 × 0.06 near a trailing
  edge) the transition tetrahedra and pyramids are very thin (dihedral angles near 0° / 180°;
  volumes positive).
- **Collision tests:** floating point with a tolerance, not exact predicates. Pyramid caps on the
  front are placed after the collision tests.
- **Sliding:** only on planar unstructured faces used by one volume. The volume's view of such a
  face differs from the face's own mesh.
- **Not in 8b (8c):** mixed blocks + volumes in CGNS; optimisation.

## Phase 8c — Optimisation, ridge fans, mixed models (2026-10-12) — approved

### Built
- **Mixed models in CGNS:**
  - One C++ writer `write_cgns_3d` (the structured-only and unstructured-only writers are now
    wrappers): structured zones first, then unstructured ones.
  - Faces shared by a block and a volume (or two volumes) get vertex connectivity both ways,
    with (i, j, k) triples on the block side.
  - Python `write_cgns_3d` replaces `write_cgns_blocks` / `write_cgns_volumes`; BCs come from
    `analyze_3d`.
- **Optimisation** (`unstructured/optimize3d.cpp`, own implementation, not mmg — see
  ARCHITECTURE §7):
  - Signed scaled Jacobian for all cell types.
  - Node moves of interior nodes and pyramid apexes (deficit-sum objective, worst-cell guard).
  - 2-3 / 3-2 swaps of core tetrahedra.
  - Generator key `optimize` (default on), stats in `mesh["optimization"]`, GUI checkbox.
- **Ridge fans in 3D inflation:**
  - Sectors per node, fan columns along convex edges turning more than 60°, zero-area fan quads,
    corner caps around centre columns, ridge ends.
  - Collapsed-bottom cells built from their faces (sideways prisms in a fan's first layer).
  - Sliding planes re-meshed from cancelled edge chains.
  - `fan_columns` reported; fan threshold and step are internal settings.
- **Demo** `examples/08c_optimisation_mixed.py`:
  - **Wing** (the 08b case): 82 % of wall nodes with all 20 layers (8b: 43 %); 18 of 25
    trailing-edge nodes (8b: none). 42 894 hexahedra, 3 668 prisms, 3 364 pyramids, 55 236
    tetrahedra. Optimisation: core tetrahedra below 0.3 from 4 745 to 3 691, below 0.1 from 640
    to 336.
  - **Duct:** structured block (512 hexahedra) + unstructured volume around a structured cube
    with 10 hexahedral layers, edge fans and corner caps; no column stops; worst hexahedron
    0.966; volume equal to analytic to 1e-12. CGNS has 2 zones joined by vertex connectivity;
    cgnscheck is clean.
  - Outputs: `08c_optimisation_mixed.png`, `phase8c_main_window.png`, CGNS / VTKHDF / msh,
    projects, journals.
- **08a/08b outputs regenerated:** new generator key; the 08b wing now has fans (85 % of its
  wing-surface nodes with all layers).

### Verified
- `ctest`: 119/119 (+5 `Optimize3D`, +1 `Inflation3D.FansAtSharpRidges`; the box-body test now
  also checks fans). Covered:
  - ideal cells score 1;
  - untangling an interior node;
  - 2-3 and 3-2 swaps, and their blocking by fixed or non-tetrahedral cells;
  - fans with all layers on interior ridges, ridge ends, a closed front and exact layer volume.
- Strict CI build and clang-tidy clean.
- `pytest`: 251/251 (+5):
  - mixed CGNS (cgnscheck, connection coordinates matched through h5py), VTKHDF and Gmsh;
  - optimisation (boundary and layer nodes unchanged, volume kept, stats consistent);
  - validation and journals.
- All demos 01–08c run.

### Not verified / known limitations
- The worst cell of a mesh with layers is usually a layer transition cell, which optimisation
  doesn't touch.
- Optimisation trades pyramid quality for tetrahedra (wing: worst core pyramid 0.34 → 0.23).
- Near the wing tip the layers still taper: the tip cap reaches the core size first.
- See `OPEN_ISSUES.md` (3D inflation, optimisation).

## Phase 8c addendum — quality floor for 3D cells (2026-10-12) — approved
Requested at the 8c review: very distorted or skewed 3D cells must not exist, and quality must be
tracked.

### Built
- `quality/quality_cells3d.cpp`: scaled Jacobian, dihedral angles and equiangle skewness of all
  3D cell types (the optimiser's metric moved here); `_core.quality_cells_3d`, equal to
  `volume_quality` to 1e-15.
- `quality_floor` (default scaled Jacobian >= 0.05, skewness <= 0.95) in the volume generator and
  block parameters; GUI rows for volumes.
- Inflation:
  - every layer cell must meet the floor;
  - soft stops (size transition, slope) are cancelled when their cells would miss it.
- Optimiser: a rescue pass (max-min of the floor margin) for core cells below the floor.
- Reporting: `mesh["quality"]`, `og.block_floor()`, status texts, warnings, boundary-report
  issues (exports refuse).
- Demo wing refined to 61 spanwise points.

### Results
- **Wing:** 0 of 210 980 cells below the floor (worst scaled Jacobian 0.050, skewness 0.936).
  Mean 17.9 of 20 layers, mostly stopping on reaching the core size; trailing-edge nodes mean
  17.2.
- **Without the floor** (the 8c wing, 25 spanwise points): 447 cells below a scaled Jacobian
  of 0.05, down to 0.003.
- **Sphere, duct, 08a, blocks:** unchanged and above the floor.

### Verified
- `ctest` +2 (rescue, metrics); every default-floor inflation test now checks all cells against
  the floor. The taper and ridge-end mechanics are tested without it, and with it (all cells
  pass, columns stop instead).
- `pytest` +3:
  - floor tracked and matching `volume_quality`;
  - unmet floor → status, report issue, export refused;
  - validation and journal;
  - block floor (sheared extrusion).
- All demos run with their exports, so every demo mesh meets the floor.

## Phase 9a — STEP / IGES import through OpenCASCADE (2026-10-12) — approved
Phase 9 split into 9a (CAD import), 9b (scripting polish, golden cases from scripts, split / join
and O-grid / butterfly helpers for multi-block hexahedral meshes, as agreed at the 8c review) and
9c (performance, packaging).

### Built
- `backends/occt` (`og_occt`, option `OG_WITH_OCCT`) and module `opengrid._occt.read_cad`: STEP /
  IGES read, ShapeFix, sewing (when free edges exist; at the larger of the requested and the
  file's tolerances), small-edge removal, splitting of periodic faces and closed edges; B-spline
  curves (ends set to the vertices) and trimmed B-spline surfaces; topology with wires, edge use
  counts and a healing log.
- `og.import_cad(path, spacing, max_turn_angle, domains, ...)` → `CadImport`: curves, surfaces,
  connectors, domains ("unstructured" / "auto" / "none"), skipped faces with reasons, free and
  non-manifold edges, `report()`. One undo step; the journal holds explicit `create_*` calls.
- `CreateEntities(optional_kinds=...)`: invalid domains are skipped and reported instead of
  failing the import.
- GUI: File → Import CAD (STEP / IGES), report in the log.

### Results (demo `09a_cad_import.py`)
- Bracket (block with four filleted edges and a hole, written by OCCT): 30 curves, 6 surfaces,
  12 domains (10 structured with "auto"), 0.05 s import.
- Solid: 12 441 tetrahedra; volume within +3.3e-4 of the CAD volume (faceting).
- Fluid around it: 12 layers at every wall node, 115 979 cells, 0 below the quality floor
  (worst scaled Jacobian 0.082, skewness 0.917), 6.8 s; CGNS passes cgnscheck.
- IGES of unsewn faces: 30 free edges sewn at 1e-5 to a closed shell (34 free edges without
  sewing).

### Verified
- `ctest` +4 (STEP topology and geometry to 1e-9, end poles = vertices, lengths; IGES sewing;
  no splitting; scale and errors) → 125.
- `pytest` +9 (model; entities on the CAD geometry; turning-angle dimensions; "auto" structured
  faces; option checks; IGES sewing; sphere faces reported; imported solid volume; GUI) → 263.
- Strict CI build, clang-tidy on the new files, ruff, clang-format, ternary check; all demos
  01–09a run and earlier outputs are unchanged.

## Phase 9b — scripting documentation, golden cases, topology editing (2026-10-12) — approved

### Built
- `docs/SCRIPTING.md` (guide) and `docs/API.md` (reference generated by `tools/gen_api_docs.py`
  from the Session docstrings; the 30 undocumented methods got docstrings); README updated.
- `examples/golden/g01`–`g10`: the ten golden cases of SPEC.md as standalone API scripts with
  invariant checks.
- Explicit distributions (C++ `DistributionType::Explicit`, Python `{"type": "explicit",
  "positions"}`), `CurveChain::locate_length`.
- `og.split_connector`, `join_connectors`, `split_domain`, `join_domains`, `split_block`,
  `join_blocks`, `ogrid_domain`, `ogrid_block` (`model/topology_edit.py`, command
  `TopologyEdit`); GUI Topology menu.
- `.gitignore`: `dist/` → `/dist/`; the distribution sources (never committed before) are now
  tracked; a copy of the tracked files builds and passes ctest.

### Results (demo `09b_topology_golden.py`)
- Pipe (radius 0.5, length 3): O-H section 0.05 → 0.71 min scaled Jacobian; 5 extruded blocks,
  153 600 hexahedra, 8 interfaces; CGNS 5 zones with 16 connections, cgnscheck clean.
- Solid sphere: one block −0.031 (inverted) → 7 blocks, 32 000 hexahedra, min 0.577; CGNS 7
  zones, 36 connections, cgnscheck clean.
- Duct split at k = 20 into two 21³ blocks and joined back to the same points.
- Golden cases: all ten pass (about 13 s), journals replay
  (`examples/output/09b_golden_report.txt`).

### Verified
- `ctest` +2 (explicit distribution, chain locate) → 127.
- `pytest` +24 (11 topology editing incl. GUI, 11 golden, 2 API docs) → 287.
- Strict CI build, clang-tidy on changed C++ files, ruff (incl. Markdown code blocks),
  clang-format, ternary check; all demos 01–09b run, earlier outputs unchanged.

## A tidier repository root, dependencies by need, Linux machines without internet (2026-10-21) — requested by the user
The user asked why `.clang-format`, `.clang-tidy`, `.gitignore` and `CLAUDE.md` were in the
repository, found the dependency list inconvenient for their SLES 15 machine, and could not
run `install_from_source.sh` there (no curl / internet).
- The clang configs moved to `tools/clang-format.yaml` and `tools/clang-tidy.yaml` (passed with
  `--style=file:` / `--config-file=`; editors no longer find them by themselves);
  `tools/lint.sh` runs every style check. `CLAUDE.md` (Claude Code's notes) is no longer
  tracked; its conventions are in `docs/ARCHITECTURE.md` §9. `.gitignore` stays: git reads
  it only at the root.
- Dependencies: the README lists them by need (building, the API, the application, each
  optional feature, tests, in-tree code). `h5py` is loaded on first use (`opengrid.io`
  re-exports the VTKHDF writers lazily): OpenGrid imports and runs without it; meshio left
  the conda package's run requirements (tests only); the pip extras are `gui`, `export`,
  `console`, `test`.
- `install_from_source.sh` downloads with curl, wget or Python, uses a Miniforge installer
  copied by hand (`MINIFORGE_INSTALLER`), and when a download or conda fails says that the
  machine needs internet (or `https_proxy`) and points to the bundle. `docs/INSTALL_LINUX.md`:
  a machine without internet takes the bundle (copied with scp or a USB disk).
- The bundle in `packaging/dist/linux` is still the one from 2026-10-19 (no Docker builds on
  this Mac any more, at the user's request: emulated x86_64 is too slow). 554 passed, 1
  expected failure; lint clean.

## Boundary types on 3D faces, the Script tab runs code, bold selection, points before meshing, conversion of sketch domains (2026-10-21) — requested by the user
The user asked for the Script tab to be editable with new lines run in the window, selected
curves in bold (colours already mean boundary types), a preview of where points land before
meshing, said converting a structured domain to unstructured did not work, and that a click
on a face of an extruded box selected nothing: how are boundary conditions set in 3D?
- **3D faces** (`app/scene.py`, `MainWindow._select_face`). Clicks found nothing on a block:
  when no curve or point was under the pointer, the second picking pass took every actor,
  and the reference grid's lines in front of the faces caught the ray. It now picks among
  the model's entities only. Blocks and volumes are drawn from their boundary faces with a
  "face" cell tag (`Scene._boundary_faces`: a block's six faces, a volume's mesh faces: ends,
  side walls, bounding domains), so a click selects that face as a patch member (`{"block":
  id, "face": "imin"}`, `{"volume": id, "face": "top"}`, a side wall `{"connector": id}`):
  filled and outlined in bold; Shift adds faces; a second click on the selected face selects
  the block (and the face's panel has a "Select block-1" button). The face's panel shows its
  points, patch and a boundary type (Apply: `set_boundary_type` / `clear_boundary_type`);
  the Boundary panel's Set type and the right-click menu take selected faces; the Boundary
  panel's 3D rows now select their faces (`BoundaryReport.members`, `api._face_member`);
  Colour by patch fills 3D faces with their patch colours, unassigned ones red
  (`Scene.show_patch_faces`, kept current from the patches after each edit). The selection
  keeps its faces across model changes (a block that goes takes its faces along).
- **Script tab** (`app/script_panel.py`: `ScriptEditor`; `MainWindow.run_script`;
  `SessionTrace.typed`). What already ran is read-only (typing over it writes at the end;
  Backspace / Delete / Cut stop at it; the editor's own undo is off), the lines after it are
  tinted and run with Ctrl+Enter, Shift+Enter or Run: statement by statement in the
  console's namespace (IPython's: names defined in either are seen by the other), inside one
  refresh; each statement that runs joins the trace as written (not as the calls it made:
  variables and loops survive, and the trace still replays); the first that fails stops,
  and it and what follows stay for editing, with the error under the editor (what ran before
  the error in a failing statement is recorded as its calls). Printed text goes to the Log.
  The window's own actions join the trace before the lines not run yet.
- **Bold selection** (`SELECTED_CURVE_WIDTH` 9 px, was 4 like the plain connectors): selected
  curves and points-curves are drawn bold, in the selection colour, or in their patch colour
  when colouring by patch.
- **Points before meshing** (`Session.preview_points`, `connectors.preview_points`: the C++
  node distribution on normalised trial parameters, nothing stored). Selected curves without
  points show pink dots where the "points" box would put them (they follow the box); in a
  curve's panel, editing the number of points, the mode, the distribution or the breakpoints
  previews the new points before Apply.
- **Conversion** (`Session.convert_domain`). The user's domain was a sketch's, which was
  refused ("list its block in the sketch's unstructured blocks instead"); it now converts
  through its sketch (`edit_sketch(unstructured=...)`, both ways, name kept). A domain in a
  smoothing group is converted with `leave_groups=True` (it leaves the group; a group left
  empty goes), which the window asks before (`_confirm_leave_groups`).
- Tests: `test_faces_3d.py` (a block face made an inlet from its panel, an outlet from the
  Boundary panel, the right-click menu, the Boundary panel's rows selecting faces, patch
  fills; an extruded volume's top and side wall), `test_convert_domain.py` (sketch domains
  both ways with undo and replay, a grouped domain, the window's confirmation), the Script
  tab's runs, failures, read-only trace and shared names (`test_console_trace.py`), the
  preview and the bold selection (`test_curves_view.py`); the 3D picking test now expects
  faces, then the block. macOS: 553 passed, 1 expected failure; lint clean.
- Follow-up (the user asked whether Pair as periodic works with faces): it did (their case,
  replayed from its trace: both pairs matched face to face, CGNS 4 periodic connections,
  Gmsh 2 links), but the Boundary panel still showed the 2D sketch's report ("no type yet"
  on the curves) under "press Refresh". The panel now reruns a 3D report after each change
  (it reuses the meshes drawn in the view) unless the last one took more than a second
  (`SLOW_REPORT_SECONDS`), and a stale panel clears its rows. Test: the user's workflow in
  the window (`test_faces_pair_as_periodic`). 554 passed, 1 expected failure.

## One "curve" in the window, the IPython console in the lower right, the Script tab (2026-10-21) — requested by the user
The user found the database curves and the connectors redundant (a sketch's edges appeared
twice in the tree, `sketch-1-e0-1` and `sketch-1-c0-1`) and asked for one concept called
curves; the console in the lower right corner instead of the whole bottom, IPython instead of
plain Python; and the commands done by hand traced as Python in a temporary file, reachable
from the console.
- **Curves** (`model/curves.py`, `CurveIndex`). The model keeps its two layers (the database
  curve is the shape that handles, templates and surfaces use; the connector holds the points
  that domains and patches use), and the window shows a curve: a connector covering one whole
  database curve, the only one on it, is one item with that curve. The tree has a Curves
  category instead of Database and Connectors (one row per curve, named after its points,
  "line, 21 pts"; a shape's own name is added, "line rectangle-1-bottom, 21 pts", unless it
  is automatic or a sketch's), and Surfaces and Points categories that show only when not
  empty. Picking or selecting the shape selects the curve (`MainWindow.set_selection`,
  `pick_at`); the property panel shows the shape's parameters (`shape:<field>` editors, a
  sketch's read-only) and the points together, one Apply; a curve without points has a "Give
  it points" button; deleting a curve takes its shape. Points along part of a curve or along
  several curves are curves of their own ("6 pts, part of a line", "21 pts, along 2
  curves"). The actions that give points (Points on curves, Chain of curves, Domains from
  curves) take only the curves without points and say so otherwise; Structured domain now
  works on curves without points too (`Session.structured_domain_from_curves`, new: the sides
  are found on trial points in a transaction that is rolled back, then everything is created
  as one undo step, as `unstructured_domain_from_curves` did). Automatic connector names are
  `curve-<n>` (`store.NAME_PREFIX`; counters still per kind, replays unchanged). The
  window's texts say curve(s) where they said connector(s). The API keeps its names.
- **Console** (`app/console.py`). `make_console` gives Jupyter's Qt console on an in-process
  IPython kernel (qtconsole + ipykernel, now in `environment.yml`, the conda recipe and the
  `gui` extra; without them, or with `OPENGRID_PLAIN_CONSOLE`, the plain console). One kernel
  per process (IPython's shell is a singleton; per-window kernels left asyncio tasks pending
  at exit), started when a console is first shown or used, so windows that are never shown
  (the tests) start none; each console binds `og`, `window`, `opengrid`, `trace_path` again
  when it gets the focus or runs code. Found by the keyboard test: qtconsole's text area does
  not claim plain keys, so `2` typed in it switched the view; the console now accepts the
  shortcut override of keys without Ctrl / Alt / Meta, as a line edit does. `run_line` returns
  the printed text, the value or the traceback (from the kernel's messages). The tip IPython
  prints at start is off.
- **Script tab** (`app/trace.py`, `app/script_panel.py`). `SessionTrace` writes the session's
  journal to `<temp>/opengrid-traces/session-<time>-<pid>.py` as it grows (a week's traces
  kept), continued across File > New (`og.new()`) and Open (`og.load(...)`; `Session.load`
  now puts that line in the new journal before telling the listeners), plus the exports and
  saves the window makes (`note`). The Script tab shows it with Python colours, Copy, Save
  as… and Open file; `trace_path` in the console. The trace replays (`opengrid.replay`) to
  the session it followed (test).
- **Layout.** Console, Script and Log are tabs of one dock in the right area, under Quality
  / Smoother / Boundary; the viewport keeps the window's height; default widths 360 / 460 px
  (`window/state_v3`: the saved layout of earlier versions is not restored). Three tabified
  docks were tried first: Qt 6.11 leaves a stale tab bar over the viewport when a second
  group of tabified docks is formed before the window is shown (reproduced in a bare
  QMainWindow, also on the native macOS platform), hence the tab widget.
- Tests: `test_curves_view.py` (pairing rules, texts, automatic names and replays,
  `structured_domain_from_curves` with undo / redo / replay / refusal, the one-row selection,
  the combined panel moving the shape and changing the points in one Apply, Give it points,
  deleting a curve, Structured domain on curves in the window), `test_console_trace.py` (the
  trace across New / Open / save and its replay, old traces removed, the Script tab after
  drawing, saving and typing, the dock layout, the IPython console, the plain fallback);
  existing tests follow the new names and texts. macOS: 544 passed, 1 expected failure; lint
  clean. The Linux bundle in `packaging/dist/linux` predates this round (rebuild it with
  `build_bundle.sh` to get the curves, the IPython console and the Script tab there).

## Repository cleanup, README, complete examples, Linux (x86_64) check (2026-10-19) — requested by the user

The user asked whether the clang files and `.github/workflows` were needed, for the README to
reflect the latest changes, whether the examples were representative (also for scripting),
and how to be sure OpenGrid compiles on openSUSE 15 on an Intel Xeon (x86_64).

### Done
- `.github/workflows/ci.yml` removed: it never ran (GitHub Actions is disabled for the private
  repository), its macOS job would cost 10x the Linux minutes, and its Linux job had no C++20
  compiler for distributions with an old system GCC. `.clang-format` and `.clang-tidy` stay at
  the root: clang-format and clang-tidy (and editors) find them by walking up from each source
  file, so they must sit above `core/`, `bindings/` and `backends/`. `.gitignore` regrouped (its
  comments had been sorted away from their patterns).
- `packaging/linux/check_build.sh [image]`: builds and tests OpenGrid in a clean container
  (openSUSE Leap 15.6 on x86_64 by default) as a user would: Mesa, Xvfb and fonts from the
  distribution, Miniforge, `environment.yml` plus conda-forge's `cxx-compiler`, the `ci` preset
  (warnings as errors), ctest, the Python package built from the same build directory, pytest
  under `xvfb-run`; Miniforge cached in a Docker volume; PASS / FAIL as the log's last line.
- Two complete examples on the newest API: `12a_airfoil_cases.py` (a NACA 2412 in a far field
  and a periodic NACA 6412 cascade, 40 inflation layers from 2e-5 chord, boundary types,
  `make_periodic_pair`, `extrude_model`, CGNS / Gmsh; 2 s) and `12b_bodies_in_a_box.py` (a
  sphere and a closed swept wing in boxes, `mesh_surfaces`, volumes with 12 / 8 prism layers
  and a hexahedral trailing-edge fan, wall / far-field patches; 10 s); figures saved as
  128-colour PNGs (a third of the size). `examples/README.md` lists every script by topic;
  `tests/python/test_example_cases.py` builds both cases so they cannot rot.
- Every example was run from a scratch copy: all 30 scripts and 10 golden cases exit cleanly
  (`00_launch.py` opens the application and waits, by design; the automatic topology of the
  two-element airfoil in `10e` fails as recorded in OPEN_ISSUES and its strict xfail).
- README: the newest features moved to the sections they belong to (they had been appended
  under Boundary conditions), a shapes-to-solver snippet (run as written), the gallery's new
  row, the examples by topic, test counts, a Linux / openSUSE build section with the packages
  and compiler minimums verified in the container, the status. SCRIPTING: the shape and
  surface-mesh sections moved under Typical workflows (they sat under Commands), a stale
  sentence about the deleted pole cleanup replaced, a pointer to the complete examples.

### Linux (x86_64): build, tests, bundle, installer
The user's target is openSUSE 15 on an Intel Xeon, with nothing to install by hand (Miniforge
fine if a script installs it). Everything below ran in Docker containers on the Mac, emulated
x86_64 (Docker Desktop raised to 12 GB and 8 CPUs for it).
- **The first build found real portability bugs** (conda-forge's GCC 15.3, the `ci` preset,
  warnings as errors; Apple Clang never showed them): three C++ tests used `std::sort`,
  `std::clamp`, `std::reverse`, `std::min_element` without `<algorithm>` (libc++ includes it
  transitively, libstdc++ does not); `__int128` under `-Wpedantic`; a sign conversion in
  `bindings/og_3d.cpp` (nanobind's `shape` takes signed sizes); `-Wshadow` warnings inside
  nanobind's headers (now system headers: `NB_SUPPRESS_WARNINGS`); GCC false positives:
  `-Warray-bounds` in `std::sort` on part of a `std::array` (`volume_quality.cpp`, now an
  insertion sort of the 3-4 face nodes) and `-Wstringop-overflow` on vector assignments in a
  test (now `std::copy`). CMake 4 scanned every file for C++20 modules (off:
  `CMAKE_CXX_SCAN_FOR_MODULES`). With conda-forge's `cxx-compiler` 2.0 (plain `c++`, no
  activation) CMake found no `ar` in a minimal system: `CMakeLists.txt` now stops with a
  message naming binutils or the prefixed compiler. After the fixes: 0 warnings, ctest 153 of
  154 (the GPU test skipped: no OpenCL device). macOS: 154 of 154, 0 warnings.
- **Python suite on Linux** (Leap 15.6, xvfb-run), first full run: 520 passed, 8 failed, 3
  skipped, 1 expected failure, in 3 h emulated. The failures: macOS tar's `._` files in the
  copied source (the ternary check read one; the scripts now set `COPYFILE_DISABLE=1`); five
  real-event tests whose keyboard shortcuts never arrived because, without a window manager,
  the window never became active (the harness now calls `activateWindow()` and waits; Qt then
  sets the X input focus itself; making Qt's offscreen platform the Linux default instead was
  tried and abandoned: pyvistaqt can abort when a window closes there); a boundary-panel test
  racing the panel's refresh timer (waits for the new count); the worker test counting timer
  ticks over 0.3 s, which a slow machine misses (the worker now waits for five ticks: no
  timing); and the automatic topology of the symmetric bracket case, which decomposes on macOS
  but not with glibc's `cos` / `sin` (last-bit differences in the boundary samples amplified
  by refinement and tracing; deterministic per platform: an expected failure on Linux,
  OPEN_ISSUES). Rerun with every fix: 527 passed, 2 expected failures (the two-element
  airfoil everywhere, the bracket on Linux), 3 skipped (no OpenCL device, the macOS launcher
  check, no `gmsh` executable), none failed, in 1 h 29 min emulated.
- **The conda package and the self-contained bundle**: rattler-build compiles the recipe with
  conda-forge's sysroot pinned to glibc 2.28 (`variants.yaml`; without it `_core` needed glibc
  2.34, which Leap 15.3-15.5 lack); its own tests (imports, OpenMP, the smoke test with a
  Netgen volume) pass. The bundle: 833 MB, 2.9 GB unpacked; its binaries need no glibc symbol
  newer than 2.17. `test_bundle.sh` in clean containers (no conda, no compiler, only Mesa,
  fonts and Xvfb): openSUSE Leap 15.6 (glibc 2.38), Leap 15.5 (2.31) and Rocky Linux 8 (2.28)
  all PASS: the imports, `12a` (2D), `12b` (Netgen volumes), `09a` (STEP / IGES), the
  application window, a project opened and coloured by quality (the same 0.36 / 0.92 scaled
  Jacobian as on macOS), a moved folder, `install.sh`.
- **`install_from_source.sh`** in a clean Leap 15.6: Miniforge downloaded into
  `~/.local/share/opengrid/miniforge`, the environment, conda-forge's GCC 15.3, the build,
  the launchers; `opengrid-python` then ran `12a`: PASS in 11 minutes emulated.

## Spheres mesh in a conformal plane (2026-10-18) — the user: "stuck here, too slow"

Meshing the sphere at the suggested spacing (about 4 000 triangles per half) hung: the
parameter-space mesh crowds the poles (3x the triangles needed), and the cleanup added the
day before is Python with quadratic cost.

### Built
- `frames.RevolutionFrame` (method in the module docstring): for a surface of revolution
  whose generating curve meets the axis, the isothermal coordinate sigma = integral ds / r
  tabulated along the curve, w = exp(sigma + i theta) and the Mobius map z = (w - 1) /
  (w + 1), scaled to 3D length units on average; exact both ways (the far pole handled by
  its 1 / (1 - t) asymptote). `make_frame(..., surface_entity, store)` picks it for such
  surfaces of less than a full turn; `domain_frame` passes the entity.
- The on-surface cleanup of the day before (collapses, smoothing, flips in Python) is
  removed: a mesh algorithm in Python breaks the project's rule (grid algorithms in C++),
  and the conformal plane makes it unnecessary; the pole opening (index bookkeeping around
  the C++ mesher) remains as the fallback for other frames, without cleanup.

### Verified
- A unit sphere at spacing 0.15: about 850 triangles per half (645 ideal; 2 870 crowded, 490
  after the slow cleanup), smallest angle 19 degrees, in under 0.1 s; the user's case
  (radius 0.5, spacing 0.033): 4 200 per half in 0.1 s instead of minutes; the sphere in a
  box meshes to a volume with no cell under the quality floor; the frame's round trip,
  poles and middle, flip, and the plain frame for a cylinder half (`test_surface_mesh.py`).

## Menu actions no longer pass Qt's "checked" flag into their slots (2026-10-18) — bug reported by the user

Create > Mesh surfaces did nothing: `QAction.triggered` emits a ``checked`` bool, and the
helper connected the slot directly, so ``mesh_selected_surfaces(spacing=False)`` ran and was
refused ("spacing: must be a positive length, got False"); the same hit Extrude 2D model
(``thickness=False``) and File > New (``confirm=False``: no "discard changes?" question).
`MainWindow._action` now connects through a wrapper that drops the flag. Test: the three
actions triggered from their menus (dialogs answered by patches) do what their slots do.

## Spheres mesh: poles opened in parameter space, cleaned up on the surface (2026-10-18) — the user

The user asked how to mesh a sphere; neither the shape nor a CAD sphere could be (its poles
leave no boundary, and the importer skips faces with degenerate edges).

### Built
- `unstructured.py`: `_open_poles` finds loop nodes where the surface's parameter line
  collapses (`_degenerate_direction`), opens each into two 2D nodes joined along the
  degenerate side (taking the neighbours' other coordinate) so the 2D polygon is the full
  parameter region; `_close_poles` folds the duplicates back after triangulation (the loop
  nodes stay first, in loop order; cells left with a repeated node go, quads with one become
  triangles). `_improve_on_surface`, for such meshes only: interior edges shorter than 0.55
  of the target collapsed (link condition, no turned triangles), interior nodes smoothed and
  projected back onto the surface, interior edges flipped where the smallest angle grows.
- `create_sphere` makes two halves ``<name>-a`` / ``-b`` of the meridian ``-meridian-a``
  with the opposite meridian ``-meridian-b``; `mesh_surfaces` meshes them (the pole arcs are
  dropped as zero-length sides, the two meridians shared). A full-turn surface closing on its
  axis says to build it as two halves.

### Verified
- pytest: a unit sphere at spacing 0.15 gives about 490 triangles per half (2870 before the
  cleanup), smallest angle above 21 degrees (1 degree before), all nodes on the sphere,
  manifold, meridian nodes first, fans at the poles; a sphere in a box meshes to a two-shell
  volume with no cell under the quality floor (473 before); replay; the shape and
  surface-mesh tests updated to the new sphere parts.

## Surface meshes of the 3D shapes (2026-10-18) — requested by the user

The user tried "Unstructured domain" on a cylinder's generating lines, which are no loop;
the 3D shapes had no way to be meshed.

### Built
- `model/surface_mesh.py` (method in its header): a surface's four parametric sides as curve
  specs (plane lines; ruled curves and their joining lines; an extruded curve, its translate
  and the extrusion lines; a revolved curve at the start / end angles as rotated copies and
  the arcs swept by its ends), zero-length sides dropped; sides matched to existing curves by
  sampling (a full circle by arcs tiling it) or made; connectors on them reused or made from
  a spacing / point count; flat surfaces as planar domains, curved ones as domains on the
  surface (loop reversed when clockwise). `Session.mesh_surfaces(surfaces, spacing=None,
  points=None, generator=None)` (one undo step), GUI Create > Mesh surfaces (spacing asked,
  a thirtieth of the selection's size suggested; the domains selected for Volume from faces).
- `create_cylinder` makes the side as two half surfaces between seam lines and rim arcs (a
  domain cannot wrap a full turn); the discs stay, meshed as planar domains on the rim arcs.
- Refusals with the reason: a full turn of a curved surface, a sphere (poles), NURBS
  surfaces, non-surface entities.

### Verified
- pytest +4 (`tests/python/test_surface_mesh.py`): a box (12 shared lines and connectors,
  planar domains, undo, reuse with ``points``, replay), a cylinder in a box (rims found, the
  disc's loop on the rim halves, the half sides on the surface, a volume with two shells of
  over a thousand tetrahedra, the full-turn refusal), a wing (leading / trailing edge lines
  shared, flat caps as two-curve loops, a volume of the wing alone), an extruded surface
  (the translated arc made), the sphere refusal, the GUI action feeding Volume from faces.

## Predefined shapes are closed bodies (2026-10-18) — requested by the user

The wing had an open tip. `create_wing(caps=True)` now adds the flat ``<name>-root`` and
``<name>-tip`` sections (ruled surfaces between each section's upper and lower curves), the
dialog has the switch; the box, cylinder (caps) and sphere were closed already, the 2D
shapes are closed loops. Test: the cap lies in the tip plane; ``caps=False`` leaves it open.

## Selecting inside a box: lines first, hide selected (2026-10-18) — requested by the user

The user could not select a wing inside its far-field box: the box's faces took every click.

### Built
- `Scene.pick_display` picks in two passes: the line and point layers and the sketches
  first (a curve seen through a surface is what the click means), then everything.
- `Scene.hidden` / `set_hidden`: hidden entities are neither drawn (layers and per-entity
  actors) nor picked, through rebuilds; View > Hide selected (H) and Show all (Shift+H),
  `MainWindow.hide_selected` / `show_all`, with the count in the status bar.

### Verified
- pytest +2 (`tests/python/test_gui_hide_pick.py`): a root spline of a wing picked and
  clicked through the box's faces, the box face picked beside it; hiding the box (actors
  gone, clicks land on nothing, hidden through a model change), showing all, the H key.

## Standard shapes (2026-10-18) — requested by the user

The user asked for default shapes: circles, NACA 0012, cylinder, sphere, wing.

### Built
- `model/shapes.py` (methods and references in its header): `naca4` (4-digit sections,
  cosine spaced, closed trailing edge), `placed` (chord, angle of attack, leading edge),
  `ellipse_nurbs` (the nine-point rational circle scaled), `rectangle_corners`, `frame`,
  `box_planes` (outward normals), `wing_sections` (taper, sweep, dihedral).
- `Session.create_rectangle / create_ellipse / create_airfoil / create_box /
  create_cylinder / create_sphere / create_wing`: one `CreateEntities` step each, parts
  named `<name>-upper`, `-side`, `-xmin`, ... (`<kind>-<n>` by default, chosen without
  touching the name counters so replays match).
- GUI: Create > Shapes, a dialog per shape (`app/shapes_dialog.py`: fields parsed as
  numbers, vectors, text, yes/no; invalid input kept in the dialog), the result selected.

### Verified
- pytest +4 (`tests/python/test_shapes.py`): the section formulas (thickness, symmetry,
  camber, closure), the ellipse / frame / box / wing helpers, every creator (names,
  geometry, an airfoil in a far-field rectangle meshed with the airfoil as a hole, undo as
  one step, replay), the GUI (create_shape, the menu, the dialog's parsing).

## Closed loops stay closed; the region between two loops (2026-10-17) — the user

The user could not mesh a freehand loop after moving its control points: dragging the
handle at the seam (the first or last control point of a closed curve) moved one end only
and opened the loop. They also asked how to mesh between an inner and an outer curve.

### Built / fixed
- `topology_edit.move_handle`: a closed curve's seam (first and last defining points
  coincide) moves as one point, so the loop stays closed; open curves unchanged.
- The loop-assembly message for a curve whose two loose ends are close (within a tenth of
  its length) now says "its two ends are 0.06 apart: close the loop (drag one end onto the
  other, or release a freehand stroke at its start)".
- Meshing between two loops already worked (select both, "Unstructured domain": the largest
  loop is the outer boundary, the others holes); a test locks it in.

### Verified
- pytest +3 (`tests/python/test_closed_loops.py`).

## The camera stays put while dragging (2026-10-17) — bug reported by the user

Every control-point drag reset the camera, too far out, and mirrored the view. Cause: the
drag preview takes the affected entities out of the merged layers and hides their actors, so
with a small model the scene looked empty while the drag's command ran; the window then
"fitted the first geometry", and `Scene.fit` went through the view preset whose default
toggles to the far side when the view already faces the plane - hence the mirror image.

### Fixed
- `Scene.has_geometry` counts a drag preview as geometry; `Scene.fit` asks for the preset's
  own side explicitly (only the view keys toggle).

### Verified
- pytest +1 (`tests/python/test_gui_camera_stability.py`): the camera state is identical
  after a control-point drag; fit twice gives the same view from +z; the view key flips.

## Handles for what moves the selection (2026-10-17) — requested by the user

The user could not move the corners of a meshed sketch block: a click in the mesh selects
the domain, which had no handles (only the sketch itself offered its vertices).

### Built
- `handles.handles_for`: a selected domain (structured or unstructured) offers the handles
  of the curves under its connectors, like a connector did; anything generated by a sketch
  (its domains, connectors, curves) offers the sketch's vertices instead, since the sketch
  is what moves them (its own generated lines are not offered: moving them is refused).
  Dragging a vertex regenerates the sketch and remeshes, inflation layers included.

### Verified
- pytest +2 (`tests/python/test_handles_selection.py`): the handle sets per selection, and
  a drag of a corner through the window's handlers with the domain selected (vertex moved,
  remeshed with layers, journaled as move_handle, undo).

## A sketch block's generator edited from its domain (2026-10-17) — requested by the user

The user could not change the smoothing, algorithm or element type of a sketch's
unstructured block: the domain's panel edits were refused as "generated by sketch".

### Built
- `topology_edit.set_unstructured_generator` / `Session.set_unstructured_generator(name,
  generator=None, **changes)`: a plain domain's generator as before (now one topology
  command, names in the journal); a sketch's block sends the triangulation settings to the
  sketch's `generator` and the inflation to walls on the block's edges (first height, layers,
  growth; the block's own walls not listed any more dropped, walls on edges shared with
  another block kept); a total thickness, layers past the core size or a connector of
  another block refused with the reason. The property panel routes an unstructured
  domain's generator through it and says so under the editor for sketch blocks.

### Verified
- pytest +2 (`tests/python/test_sketch_block_generator.py`): the sketch's generator and
  walls, the domain's resulting inflation and mesh layers, dropping, undo, the refusals,
  replay; a plain domain and the panel (smoothing, element type, a wall with the suggested
  first height).

## A connector's points follow through structured domains (2026-10-17) — requested by the user

The user could not change the number of points of a connector once a structured domain stood
on it: the edit was refused because the opposite edge kept its count.

### Built
- `model/dimension.py` (`propagate_dimension`): breadth-first through every structured domain
  using a changed connector, the opposite edge brought to the same count by its only free
  connector, else the longest (warning); contradictions (every opposite connector fixed and
  still unequal), counts below 2 and sketch-generated connectors refused with the reason.
- `topology_edit.set_dimension` / `Session.set_dimension(name, dimension)` (replacing the old
  one-connector version): applies the propagation, validates the dependents, and for a
  sketch's connector sets the edge family in the sketch (`edit_sketch(dimensions=...)`);
  returns the connectors changed. The property panel's "dimension" of a connector goes
  through it (the log says which connectors followed); `modify(dimension=...)` is unchanged.
  Follow-up (the user could not set tanh on a sketch's connector): a sketch connector's
  ``distribution`` goes to the sketch's per-edge ``distributions`` the same way, from the
  panel too; other parameters of a sketch's connector are refused saying what can be set.

### Verified
- pytest +4 (`tests/python/test_dimension.py`): stacked domains (the shared connector and the
  far edge follow, undo / redo, replay), a split edge (the longest connector, the warning,
  the reverse direction, the "at least 2" refusal), sketch families and the refusal for a
  plain domain on a sketch connector, the panel; `test_domains.py` updated (the propagation
  instead of the refusal).

## Deleting with dependents, membership never blocking, sketch blocks removable (2026-10-17) — requested by the user

The user could not delete connectors, sketches or domains: every deletion was refused
because something used the entity (a domain its connector, a patch its members, a group its
domain), with the reason in the status bar only.

### Built
- `commands/geometry.py`: `delete_plan(store, names)` (the entities plus the sketches'
  outputs, the transitive dependents, the refusals) and `DeleteEntities(names, dependents)`:
  with ``dependents`` the dependents go too, without it they are named ("delete it too, or
  delete with dependents=True"). `model/membership.py`: patches, zones and smoothing groups
  lose deleted members instead of blocking (a group left without domains goes; a periodic
  patch losing its partner is unpaired); sketch regeneration uses the same rule for outputs
  no longer generated. `Session.delete(names, dependents=False)`,
  `Session.deletion_dependents(names)`.
- `sketch.remove_blocks` / `Session.remove_sketch_blocks(name, blocks)`: blocks out of a 2D
  sketch in one regeneration (remaining blocks renumbered, edge-keyed items of edges left in
  no block dropped, unstructured indices remapped, vertices kept; the last block refused).
  `sketch.output_owner(store)`.
- GUI `delete_selected`: a sketch's domain removes its block (the whole sketch when every
  block is selected); a sketch's connector or curve gives a hint; other entities with
  dependents ask for confirmation listing them (`_confirm_deletion`), then go with them;
  a status line says what was deleted.

### Verified
- pytest +4 (`tests/python/test_delete.py`): dependents named then deleted (undo / redo,
  replay), membership removal and unpairing, sketch outputs and block removal (walls,
  pinned, unstructured, patches, renumbering), the GUI flow (hint, block removal, whole
  sketch, confirmation cancel / accept, a domain alone); the old refusal expectations in
  `test_groups.py`, `test_sketch.py` and `test_patches.py` updated.

## The Boundary panel made plain (2026-10-17) — requested by the user

The user found setting boundary conditions confusing: the panel's report was a list of
raw messages (cut off, "3 problem(s)"), the four periodic patches were named periodic,
periodic-2, periodic-3, periodic-4 with no sign of which paired with which, an "interface"
patch on a shared connector was reported as an error, and a mismatched pair said only
"closest target is 0.012 away".

### Built
- `app/boundary_dock.py` rewritten: a status line (green "All 6 boundaries have a type,
  2 periodic pair(s) match: ready to export", red "2 of 6 boundaries have no type; 1
  periodic pair(s) do not match"), then groups of rows - Boundaries (boundary, patch, type,
  status), Periodic pairs (a ↔ b, the transform, match or why not), Interfaces (shared by
  two domains: no type needed), Empty patches, Other problems; a row click selects the
  entity. Actions on the selection: Set type (type combo, optional patch name), Pair as
  periodic, Clear type; New zone; Refresh - in three short rows, so the panel fits in
  300 px (the user: "too many buttons horizontally"). The panel follows a 2D model at every change
  (`MainWindow._refresh_boundary`, coalesced by a timer); a 3D model marks it stale and
  Refresh recomputes. Setting a type or pairing raises the panel and switches the patch
  colouring on. `BoundaryReport` gained `assignment` (label -> patch) and `ids` (label ->
  entity to select), filled in 2D and 3D.
- `model/periodic.py`: a connector without an image says why: the point counts differ
  (both counts given), the ends land off the target ("not images of each other under the
  translation: snap to the reference grid or correct the transform"), or the ends coincide
  and the distributions differ.
- `patches.make_periodic_pair` names a pair `periodic-1a` / `periodic-1b` (`<name>-a` /
  `<name>-b` with a name; `pair_names`); the model tree shows "periodic ↔ periodic-1b, 1
  connector" or "periodic, NOT PAIRED" (`model_tree.patch_summary`, `patches.partner_of`).
- `patches.refuse_shared_connectors`: a connector used by two domains (an interface) is
  refused as a patch member by `set_boundary_type`, `make_periodic_pair`, `create_patch` and
  `add_to_patch`, with the two domains named; `Session.clear_boundary_type(members)`.
- Grid snapping: with the reference grid shown (G), 2D clicks snap to its minor lines when
  no entity point is near (`MainWindow.snap_point`; `ReferenceGrid.divisions`), so sketched
  blocks have exact, equal sides.

### Verified
- pytest +6 (`tests/python/test_boundary_panel.py`): the three mismatch reasons, pair
  names, the refusals (four entry points), clear + tree labels, the live panel (rows, named
  and typed patches, pair row, clicks, clear, mismatch, 3D stale + Refresh), grid snapping;
  the dock tests in `test_gui_events.py` / `test_gui_patches.py` moved to the new widgets.

## Periodic pairs in one step; the 2D model extruded one cell thick (2026-10-17) — requested by the user

The user asked how to set periodic pairs, found the patch-first route a detour in 2D (the
sides are connectors) and suggested the one-cell extrusion with periodic front and back that
3D solvers want for 2D cases.

### Built
- `Session.make_periodic_pair(source, target, ...)` (`patches.make_periodic_pair`): the two
  sides' members become two periodic patches (``name`` / ``name-2``) paired in one undo step,
  the translation from the sides' centroids unless given (`_members_centroid` handles
  connectors, domains, block faces and extruded ends). GUI: "Pair as periodic" in the Create
  and right-click menus - two connectors, two faces, one extruded volume (its ends), or two
  selected periodic patches; the source patch is selected afterwards so the panel shows the
  pairing.
- `Session.extrude_model(thickness, layers=1, ends="periodic")` (`topology_edit.extrude_model`):
  every domain extruded along z (blocks from structured domains with the edge-to-face map,
  mirrored blocks included; volumes of prisms and hexahedra from unstructured ones), the 2D
  patches carried to the sides (block faces by edge, side-wall connectors for volumes) with
  their pairings, zones extended, the fronts and backs in "front" / "back" as a periodic pair
  or symmetry planes; refusals: a structured edge in two patches, a 3D model, taken names.
  GUI: Create > "Extrude 2D model (one cell thick)…" asks the thickness (the mean boundary
  spacing suggested) and the ends, meshes in the worker and switches to 3D navigation.

### Verified
- pytest +4: the pair from two sides (names, types, translation, report, undo / redo,
  replay); the channel-with-dome model extruded (block and volume, every patch mapped to the
  right face, the side pairing kept, front / back paired with (0, 0, 0.1), the zone, the
  shared connector an interface, a clean boundary report with two pairs, CGNS connections,
  undo back to 2D, replay and project round trip, "3D already"); symmetry ends and the
  refusals; the GUI actions.

## Boundary types from the selection, suggested values, sketch attachment, domain conversion (2026-10-17) — requested by the user

The user asked which distribution types cluster and for suggested spacings, how to give a
selected line a boundary type (the Boundary dock was not found), why two sketch blocks did
not share their edge (they had become two sketches: the tool starts a new sketch when none
is selected, and a sketch never adopts another sketch's connectors), why a domain cannot
change kind, and for suggested first layer heights.

### Built
- `Session.set_boundary_type(members, patch_type, name=None)` (`model/patches.py`): the
  members join the patch called ``name`` (made when absent), else the first patch of that
  type, else a new patch named after the type, leaving the patches they were in. GUI: a
  "Boundary type" submenu in the right-click and Create menus (the patch colouring switches
  on), and a "boundary" row in a connector's property panel (the patch it is in and a combo).
- Distribution editor suggestions: the uniform spacing from the connector's length and point
  count (following the dim box), and on a type change the fields filled with a clustering to
  start from (tanh / sine: a quarter of the uniform spacing at the clustered ends; geometric:
  growth 1.2), placeholders naming the suggestion; extrusion editors in fractions.
- Inflation suggestions: the unstructured and volume editors take the mean boundary spacing
  and suggest a first height of a twentieth of it (placeholder and hint; used when the field
  is empty); the sketch panel says the same for wall first-cell heights.
- Sketch tool: a block drawn onto the vertices of an existing sketch extends that sketch
  (`attached_sketch`), and with one sketch in the model new blocks go to it; so the shared
  edge is one connector.
- The toolbar's "dim" box is now "points", with a tooltip listing what uses it (the user:
  "dim is misleading").
- `Session.convert_domain(name)` (`topology_edit.convert_domain`): structured <->
  unstructured on the same connectors, name kept, zones and 3D patches follow; refusals for
  sketch domains, group members, faces of blocks / volumes and holes. Topology menu and the
  context menu.

### Verified
- pytest +6: the panel's suggestions and boundary row (tanh / sine / uniform fills, the hint
  following the dim box, wall -> inlet -> none moves), the selection's boundary types (new,
  existing and named patches, a member moving, the unique name when "inlet" is a curve, the
  menus, replay), the sketch block attached without a selection (one sketch, seven
  connectors), conversion both ways with zone and patch following, undo / redo, replay and
  project round trip, and every refusal.

## Unstructured sketch blocks (2026-10-17) — requested by the user

The user asked why sketch blocks are limited to structured grids and whether wall layers
exist only there (they do not: unstructured domains have inflation layers; structured grids
add the orthogonal elliptic wall treatment), and asked for unstructured blocks in sketches.

### Built
- Sketch parameters ``unstructured`` (block indices) and ``generator`` (`model/sketch.py`):
  those blocks become unstructured domains on the sketch's connectors; a wall on one of their
  edges becomes inflation layers from that connector (`_block_generator`: the wall's
  first-cell height - required -, layers, growth); the sketch's generator (no inflation of
  its own) is shared by them. They are not group members; connectors between a structured
  and an unstructured block are pinned; with ``extrude`` they give volumes; a sketch with no
  structured block has no group. The vertex passes treat every block as structured. GUI:
  an "Unstructured sketch blocks" toggle beside the Sketch block tool (new blocks), the
  property panel's "unstructured blocks" list and the generator editor (walls greyed: the
  sketch's walls give them). API: `create_sketch(unstructured=..., generator=...)`.
- Reference grid made 3D (the user: "only 2D and only in the z plane"): an axis view shows
  the plane facing it (XY, XZ or YZ) with edge labels, a tilted view the three axis planes
  through the round point nearest the focal point, the facing one with minor lines, with
  billboard labels along the axes; the toolbar is back on two rows with 20 px icons (the
  single row overflowed) and the saved window layout key was bumped so the new default shows.
- Fixed on the way: a smoother result arriving after its group's members had changed kind
  (a block made unstructured while the smoother ran) crashed the completion slot; the dock
  now discards results whose members are gone and says so, and the window's completion
  handler skips entities that no longer exist.

### Verified
- pytest +2: a four-block square with blocks 1 and 3 unstructured (kinds, the group's
  members and pinned connectors, the shared connector's nodes in both meshes, the wall's
  inflation layers on the right connector with the given height / layers / growth, the
  generator applied, the refusals, extrusion to blocks and volumes, back to structured, all
  unstructured with no group, project round trip, replay); the GUI toggle through the tool
  and the panel.

## Sketch edges adopt existing connectors; reference grid (2026-10-17) — requested by the user

The user drew a sketch block under an unstructured dome and found the shared boundary not
shared (the sketch made its own 41-point connector beside the dome's 51-point one), asked how
to mesh a block unstructured (the sketch tool is structured by design: draw the boundary and
use Unstructured domain), where the inflation layers are (the generator's boundary-layer rows
of the property panel), and for a reference grid to read sizes off the view.

### Built
- Sketch parameter ``connectors`` (`{(a, b): connector}`, `model/sketch.py`): a 2D edge adopts
  a connector made by hand. Its vertices go to the connector's ends (`place_vertices`), the
  domain's loop uses the connector itself (so the other side shares its nodes), the edge
  family takes its point count (`_family_dimensions`; conflicts with given dimensions or
  another adopted count, and counts not fitting `interval_multiple`, reported by edge), the
  group keeps it `pinned`, and the vertex passes use a sketch-owned copy at half resolution.
  `sketch_references` holds it (deleting it is refused); editing it, or a curve under it,
  regenerates the adopting sketches (`regenerate_adopters`, from ModifyEntity and
  `_edit_curve`). `suggestions` offers connectors under edges (ends and middle), before
  curves; `assign_sketch_geometry` and the GUI's assign action take connectors; the Sketch
  block tool adopts the connectors its snapped clicks landed on (`matching_connectors`).
- Reference grid (`app/reference_grid.py`, View > Reference grid, `G`, persisted): lines of
  the z = 0 plane at a 1 / 2 / 5 x 10^k step from the visible extent with minor divisions,
  the x and y axes marked, coordinates of the major lines along the viewport edges in views
  along z; rebuilt after every camera move (interaction and wheel observers, the zoom window,
  the camera actions).

### Verified
- pytest +5: adoption (the domain uses the connector, the family follows, nodes shared with
  the dome and the connector an interface in the boundary report, the group leaves it in
  place, editing its dimension regenerates the sketch, conflicts and refusals by edge, replay),
  suggestions and assignment of a connector under a loosely drawn edge with the project round
  trip, the GUI tool adopting snapped connector ends; the grid's step table, actors, labels in
  z views only, zoom following, and persistence.

## GUI round: navigation, speed, 3D picking, discoverability, workers, settings (2026-10-17) — requested by the user

The user could not turn the camera in 3D after starting in 2D, missed a camera reset, asked
whether a C++ GUI would be faster, and asked for an assessment of what would make the
interface more intuitive, friendlier and faster; then asked for all of it.

### Found
- 3D rotation was bound to Ctrl + right-drag only, and "Fit view" always snapped back to the
  top view. Picking intersected the view ray with z = 0 and never picked blocks or volumes.
- After every command the scene removed and re-added the actors of every connector and
  domain through pyvista's named-actor lookup, which walks the renderer's prop list: 51 s to
  put 62 connectors on a 31 x 31 lattice of lines, 2.4 s per edit of one of them (the model
  itself: 1 ms). The connector evaluator's shared-end lookup was a Python loop over all
  connector ends (O(n^2)).
- No icons, tooltips or tool shortcuts; no hover, cursor coordinates, context menu; nothing
  persisted; only the smoother ran off the GUI thread.

### Built
- Navigation modes (`NAVIGATION`): 2D as before; 3D with left-drag rotate, Shift + drag pan,
  Alt + drag zoom window (Alt tracked in the window's Qt event filter: the Qt -> VTK adapter
  forwards Ctrl and Shift only). The mode follows the model (blocks, volumes, surfaces) until
  the user picks one. `F` resets the camera keeping its direction, `Z` zooms to the
  selection, `1 2 3 7` view presets, `V` perspective; a View toolbar. Click observers are
  passive VTK observers (the camera style grabs the interactor's focus on every press it
  sees; before, a plain left-drag mode never delivered the release).
- Scene (`app/scene.py`): curves, points and connectors are merged layers (plain / selected /
  nodes / ends, entity ids as cell data), rebuilt from numpy when any of them changes; domains,
  blocks, volumes, surfaces and sketches keep one actor each, redrawn only for the dependents
  closure of a change and recoloured in place on selection; actors are removed by object.
  Multi-command actions refresh once (`batch_updates`). The connector evaluator indexes all
  ends once per invalidation (vectorised lookups).
- 3D picking: `vtkCellPicker` on the drawn actors (`Scene.pick_display`) whenever the view is
  not along z or 3D navigation is on; the 2D picker (inside a domain counts) otherwise.
- Discoverability: painted toolbar icons (`app/icons.py`, palette colour, HiDPI), tooltips
  with shortcuts, tool keys, hover highlight after the pointer rests 30 ms, cursor coordinates
  in the status bar (2D views), a right-click menu in the viewport and the tree with the
  actions that apply to the selection.
- Workers (`app/worker.py`): `run_blocking` runs a function in a thread behind a window-modal
  busy dialog while the GUI thread processes events (the core and Netgen release the GIL);
  the window's listener only records changes meanwhile. Used for volume / block meshing after
  creation, templates, automatic topology, sketches, CAD import, the high-order overlay, the
  exports and project loading.
- Settings (QSettings; tests use a throw-away ini): window layout, zoom window / nodes /
  stretching options, the dim box, recent files (File > Open recent) and the last directory.
  Enter applies in the property panel; an edit plus a rename is one undo step
  (`ModifyEntity(new_name=...)`, `Session.modify(new_name=...)`).

### Verified
- Timings on this machine (GUI attached): 62 lattice connectors 51 s -> 0.68 s; one lattice
  edit 2.4 s -> 15 ms; a spline edit on the airfoil model 62 -> 27 ms; scene update of one
  entity 40 -> 4 ms; 500 x 200 domain respacing 79 -> 88 ms (dominated by the actor rebuild
  of 100k quads). The remaining per-command cost is the quality histogram redraw (~8 ms).
- pytest: new tests for the navigation modes (zoom window in 2D, rotation and Alt zoom in
  3D, the user's choice sticking), camera reset / presets / zoom to selection / perspective,
  ray picking in an isometric view (blocks, empty space, back to the 2D picker), hover and
  coordinates, the context menu, the worker (events processed meanwhile, exceptions, cancel),
  one refresh per long operation, settings round trip with recent files, and the one-step
  rename + edit with its journal lines; the full suite.
- Follow-up (the user did not find the way back to the XY / XZ / YZ planes in 3D): the view
  buttons read XY / XZ / YZ, a preset pressed again looks from the other side (`-xy`, `-xz`,
  `-yz` in the View menu too), the status bar names the view, the single-key shortcuts are
  window-wide (the model tree lets them through its type-ahead search) and the 3D hint lists
  the keys.

## 3D periodic pairing (2026-10-17) — requested by the user

Periodic patches could only be paired in 2D (connectors); 3D faces were refused ("not
supported yet"), so the ends of an extruded volume or the sides of a block sector could not
be periodic. The user asked for it before the extruded-volume work is pushed.

### Built
- The periodic transform takes a rotation axis (`periodic.rotation_matrix`, Rodrigues;
  `rotation_axis` in the pairing, z by default, so old projects load unchanged):
  x_target = R(axis, angle)(x - c) + c + t; `pair_periodic(..., rotation_axis=...)`. The
  inverse keeps the centre and turns the translation back (the CGNS reader sees the same
  centre on both sides).
- `model/periodic_faces.py`: the nodes of the source patch's faces (block faces, volume
  faces, the ends and side walls of extruded volumes) are transformed and matched to the
  target's by a uniform grid hash within 1e-8 of the model size (`match_points`, cells of
  twice the tolerance); the nodes must map one to one and every source face onto a target
  face with the same nodes (`match_faces`). `analyze_3d` runs the pairing checks shared with
  2D (`topology.check_periodic` with a matcher) and the report lists the pairs (`FaceMatch`:
  nodes and faces). Mismatches say how many nodes have no image and where the first one lands.
- Exports: the 3D assembly carries the links in its merged numbering (`assembly_links`), so
  Gmsh `$Periodic`, VTKHDF field data and the high-order node extension work unchanged;
  CGNS 3D writes 1-to-1 periodic connections between block faces that map corner to corner
  (`one_to_one_periodic`) and point-list connections with the periodic property for the rest,
  both ways; the C++ writer takes RotationCenter / RotationAngle (axis x angle) / Translation
  with three components (2D bases keep two), the high-order CGNS export too. Periodic faces
  are connectivity, not BCs (as in 2D). Curved high-order export leaves periodic faces
  straight-sided (`curving._periodic_sources`), as it did with periodic connectors.
- GUI: the property panel of a periodic patch has a pairing editor (partner, rotation axis /
  centre / angle in degrees, translation; empty = centroid to centroid); the target says
  where it is paired from.
- Fixed on the way: the 3D assembly stored an extruded volume's face key ("top" / "side") as
  the face's domain, which crashed the curved high-order export of a patched extruded volume.

### Verified
- pytest +5 (`tests/python/test_periodic_3d.py`): the matcher (images found within the
  tolerance only, cell boundaries); an extruded disc's ends paired by the automatic
  translation (match = the 2D mesh's points and cells; CGNS 2 point-list connections whose
  Translation maps the lists, boundaries exclude the periodic faces; Gmsh `$Periodic` affine
  check to 1e-12 and a Gmsh read without warnings when Gmsh is present; VTKHDF fields; order
  2 links = vertices + edges + quad faces, images exact to 1e-12, high-order Gmsh / CGNS
  exports; replay and project round trip); a block sector turned 60 degrees about y paired by
  rotation (CGNS 1-to-1 periodic both ways, Transform (1, 2, 3), RotationAngle (0, +-pi/3,
  0), the range corners mapped to 1e-6 in single precision; Gmsh affine = the transform
  matrix; a wrong angle reported with the node count); a rotated extrusion's ends paired by a
  rotation about x through (0, -2, 0) (CGNS point lists mapped by the stored property) and
  the refusals (rotation without a translation, zero axis, a patch with other faces: node and
  face counts); the property-panel pairing, re-apply = no change, unpairing.
- The full suite, ruff, clang-format, clang-tidy on the changed C++ and the ternary check.

## Extruded volumes: prisms and hexahedra from unstructured domains (2026-10-17) — requested by the user

"Extrude domain" only took structured domains (it builds a block). The user asked why; the
answer was that the sweep of an unstructured domain is a different kind of volume, which had
no entity, and the user asked for it.

### Built
- A volume may be given an ``extrusion`` instead of ``faces`` (`model/volumes.py`,
  `Session.create_extruded_volume(domain, mode, ...)`): the same schema as a block's
  extrusion (`blocks._normalize_extrusion`, now parametrised by the domain kinds and modes
  allowed) with translate / rotate / normal (no hyperbolic marching) and a connector-style
  distribution of the layers. `model/extrusion.py` sweeps the domain's mesh points layer by
  layer (Rodrigues' formula for rotations; the plane normal of a planar domain for "normal")
  and stacks a prism on every triangle and a hexahedron on every quad; the base is mirrored
  when the sweep runs against its normal, so cells are positive; folded sweeps (an axis
  through the domain) are reported, never hidden. Quality, the quality floor, the cut plane
  and the exports (VTKHDF, CGNS, Gmsh) work unchanged on the prisms and hexahedra.
- Boundary for patches: the domain is the bottom, ``{"volume": name, "face": "top"}`` (or
  "bottom") the ends, and the domain's connectors the side walls (one quad strip each, built
  from the mesh's boundary nodes, which come first in loop order). `analyze_3d` reports them
  (`volume.top`, `volume.connector`), treats a side shared by two extruded volumes as an
  interface, and now accepts connector members in 3D patches when they are such side walls.
  So a 2D model's connector patches carry over to the one-cell-thick 3D mesh some solvers want
  (an open issue since Phase 4). CGNS boundary sections are named by the face (`face_name`).
- GUI: "Extrude domain" accepts an unstructured domain (a quarter of its size along its mesh
  normal, "dim" box layers); the property panel edits the extrusion; "Create patch" / "Add
  to patch" take the selected extruded volume as both its ends and connectors as its sides.

### Verified
- pytest +3: a quad-dominant disc extruded (5 x the 2D cells, volume = area x height to
  1e-12, faces and side strips on the connector's nodes in order, outward orientation of
  bottom / top / sides, patches on all three with a clean boundary report, CGNS 3
  boundaries, Gmsh wedge / hexahedron / quad / triangle counts, VTKHDF cell count, replay and
  project round trip); mirrored, rotated (Pappus with straight sides, 1e-9) and normal
  extrusions, a folded rotation reported, refusals (hyperbolic, a structured domain, faces and
  extrusion together, inflation); the GUI action, the property-panel edit and the patches
  from the selection.
## Curve points added and removed in the GUI (2026-10-17) — requested by the user

### Built
- Double-click (Select tool): on a polyline, spline or NURBS curve - or on a connector lying
  on one - adds a point there; on a handle of a selected curve, removes that point. Qt
  delivers the double click (VTK has no such event): an event filter on the viewport, which
  consumes it so that the click before it does the selecting. The status bar reports the
  point count and that undo restores the curve. A click on a handle without moving it no
  longer records an empty move.
- `Session.insert_handle(name, xyz)` / `remove_handle(name, index)` (journaled, one undo
  step; `model/topology_edit.py`): polylines and splines get the projected point between its
  neighbours (chord-length parametrisation) and lose point ``index``; NURBS curves get a
  control point by knot insertion at the projected parameter (`model/knots.py`, Piegl and
  Tiller A5.1: the curve is unchanged) and lose one by knot removal (A5.8: exact when the
  curve allows it, otherwise the warning says how far the curve moved; control points without
  an interior knot - the ends - are refused, with the removable ones listed).
- `move_handle`'s network bookkeeping became `_edit_curve`, shared by the three edits; with
  ``remap`` the connector segment parameters on a reparametrised spline are re-projected
  first, so connectors and domains stay where they were (crossings are solved again as for a
  drag).

### Verified
- pytest +4 -> 461 passed, 1 xfailed: polyline / spline insert and remove (on the curve, in
  order, one undo step, replay), NURBS insertion exact to 1e-12 and its removal exact again,
  a genuine removal warned, the refusals, a lattice of splines keeping its 7 interfaces and
  its domains (grids move < 0.05) after an insertion, and a real double-click on the spline
  under a connector adding, then removing, a point with the domain following.
- ruff, format, ternary check, API docs regenerated; `tools/check_launcher.py` current.

## Freehand loops meshed from the GUI (2026-10-17) — requested by the user

The user drew a freehand curve, clicked "Unstructured domain" and nothing happened: the
stroke is a curve, the action needs connectors (its only answer was a five-second status
message), and a hand-drawn loop never closes exactly (connector ends must coincide to 1e-9),
so there was no route from a stroke to a mesh.

### Built
- A stroke released within the pick tolerance (8 pixels) of its start is closed: the start
  point is appended to the samples before the fit (`fit_stroke(close_tolerance)`,
  `create_freehand(close_tolerance)`), so the loop's curves meet exactly; strokes shorter than
  four tolerances are left alone. The new curves are the selection after a stroke.
- `Session.unstructured_domain_from_curves(curves, dimension, connectors=())`: one connector
  per curve plus the given connectors, assembled into loops (`assemble_polylines`, now shared
  with `create_unstructured_domain(connectors=...)`), and the domain, in one `CreateEntities`
  (one undo step, journaled as explicit calls). "Unstructured domain" uses it when curves are
  selected (points from the "dim" box). Loop assembly errors name the connector and its open
  end ("its end at (x, y, z) meets 0 other end(s), not exactly one").
- The status bar says what is selected and what is needed when an action does not fit the
  selection (`_selection_hint`: "d is a domain: select the connectors or curves of the
  boundary loops"), for the unstructured and structured domain, connector and network actions.
- Curve networks (`model/network.py`): self-crossings of a curve are cut (segment tests of
  the tessellation against itself, Newton refinement); dangling pieces (an end that meets
  nothing: open strokes, the tails beyond a crossing) are pruned before the faces are traced.
  They used to be walked on both sides, which made a zero-area "face" that passed the
  `area > 0` test by rounding and, with `unstructured=True`, a domain with a repeated
  connector that failed to triangulate ("duplicate points"). A curve of which nothing remains
  is reported; a face whose boundary runs along an edge twice (a curve joining two loops) is
  reported and not meshed.
- `tests/python/conftest.py` forces Qt's offscreen platform on macOS whatever `DISPLAY` says
  (XQuartz exports one in every shell; on-screen test windows get resized by the window
  manager and the real-event tests' display coordinates drift).

### Verified
- pytest +6 -> 457 passed, 1 xfailed: closed strokes, self-crossing / open / loose / bridged
  networks, `unstructured_domain_from_curves` (one step, replay, error text), a loop drawn
  with real mouse events and meshed with one click, the selection hints. No C++ change.
- ruff, format, ternary check, API docs regenerated, `tools/check_launcher.py` current.

## Phase 12d — 3D speed (2026-10-17) — requested by the user, committed when its checks passed

### Measured first
- `OG_PROFILE=1` (a `PhaseTimer`, `util/phase_timer.hpp`, in the 3D group smoother's loop;
  the report goes with the result's warnings). The 1M-cell ball (10b: 7 blocks, 1.06 M nodes,
  328 iterations), M1 Max: 10 threads 21.0 s = line sweeps 8.5 s (48 %), Anderson mixing
  6.2 s (35 %), sliding projections 1.4 s (8 %), coefficient pass 1.3 s (8 %); one thread
  110.6 s with the same split (sweeps 51 s, Anderson 37 s). (macOS `sample` attached but
  recorded nothing for this process.)
- The knobs are at their optimum: Anderson depth 10 -> 439 iterations, 24.5 s; depth 30 ->
  274 iterations but 76 ms each, 20.9 s; over-relaxation omega 1.2 -> 354 iterations,
  22.9 s; 1.4 -> 387.

### Built
- Line sweeps read a compact per-node record (`SweepCoefficients`: the six line
  coefficients, the diagonal and the frozen cross term, 80 bytes instead of 104) and compact
  line tables (32-bit node ids, byte directions: 22 bytes per line node instead of 80). Both
  exact: the grids are bit-identical to before.
- Anderson acceleration (`util/anderson.cpp`): the history rings in single precision (they
  hold differences between consecutive iterates and residuals and only enter linearly);
  the dot products as four fixed-order partial sums (one accumulator was a latency-bound
  chain under strict floating point); the extrapolation column by column (each unknown
  still subtracts the columns in the same order: the same values, but streaming and
  vectorised). The chunked, index-ordered sums are unchanged: results stay independent of
  the thread count (identical digests on 1, 4 and 10 threads). Smoothed grids differ from
  before at the level of the float history (the 1M ball: same 328 iterations, same kink).

### Results (1M-cell ball, 10 threads)
| Step | Total | Sweeps | Anderson |
|---|---|---|---|
| Before | 21.0 s | 8.5 s | 6.2 s |
| Compact coefficients (exact) | 20.5 s | 8.1 s | 6.1 s |
| + float history | 20.0 s | 8.1 s | 5.7 s |
| + Anderson loops | 19.4 s | 8.5 s | 4.6 s |
| + compact line tables (exact) | 18.5 s | 7.6 s | 4.6 s |

- 12 % in all; the 10 s target of Phase 9c stays missed. What remains is memory latency:
  every j- or k-line node gathers its coefficient record, its four side neighbours and
  itself from scattered addresses (about 1.5 KB of cache-line traffic per node per
  iteration, 1.7 GB per iteration); the Thomas solves are a small part. A per-direction
  layout of the coefficients and positions (ADI style) or a structured fast path for
  interior nodes would attack that; the solver-level routes (coarse-to-fine, GPU) were
  measured in 11c / 11d and do not pay to the same accuracy.

### Verified
- `ctest` 154; `pytest` (the smoothed grids changed at the 1e-7 level: no test needed a
  change); thread independence on the ball; lint, clang-tidy on the changed C++, ternary
  check, launcher check.

## Phase 12c — repair of curved elements and rigorous validity (2026-10-17) — requested by the user, committed when its checks passed

### Built
- `model/repair.py`:
  - Validity by Bernstein-Bezier bounds (Johnen, Remacle, Geuzaine 2013): the Jacobian
    determinant of an order-p element is a polynomial (degree 2(p-1) on triangles, 3(p-1) on
    tetrahedra, 2p-1 per direction on quadrilaterals, 3p-1 on hexahedra, (3p-2) x (3p-1) on
    prisms), interpolated exactly on the lattice of that degree and converted to Bernstein
    coefficients, whose minimum bounds J from below; inconclusive bounds are tightened by
    degree elevation (three rounds). Lattices: Gauss-Lobatto points per tensor direction and a
    GLL-warped simplex lattice (the equispaced Bernstein Vandermonde is too ill-conditioned:
    an exact rational inverse of the 220 x 220 tetrahedral matrix took minutes, the
    equispaced degree-11 tensor conversion lost 1e-4). Accurate to 1e-6 through order 4 for
    every kind; hexahedra above order 4 (degree 17: 0.1 error) and pyramids (rational
    Jacobian) are sampled instead (`rigorous`).
  - Repair: elements that fold, fall below a scaled Jacobian of 0.25 or lose 30 % of their
    straight-sided quality are optimised with their neighbours as guards - the free
    high-order nodes (not on the geometry or the boundary, never the linear vertices)
    minimise sum_q (r-1)^2 + alpha / h(r), r = J / J_straight, h the regularised Jacobian
    of Escobar et al. (rounds of decreasing delta only while something is folded), with
    L-BFGS, analytic gradients (cofactors) and BLAS-shaped evaluations; capped at 400
    evaluations; kept when the worst element improved.
- `curve_boundary(..., repair=True)` (default; `Session.high_order(order, curved, repair)`,
  the exports): assess, repair the flagged elements, assess again, then check every
  element whose sampled scaled Jacobian is below 0.3 rigorously; folded and undetermined
  elements are reported with their Bezier bound. Report: `curving["repair"]` (flagged
  elements, free nodes, iterations, patch, before / after, kept) and `curving["validity"]`.

### Results
| Case | Without repair | With repair | Time |
|---|---|---|---|
| Triangulated sphere, prism layers, order 3 | 0.638 (straight 0.966) | 0.705 (patch 0.638 -> 0.744) | 1.4 s |
| Structured sphere faces, hexahedral layers, order 3 | 0.623 (0.712) | 0.663 (0.623 -> 0.687) | 2.2 s |
| Coarse ring over 1e-5 layers, order 4 | 0.823 (0.915) | not flagged | - |
| One block on a sphere, order 4 | folds | still folds (reported, 10 s) | - |

- Bounds on real meshes are tight: ring 0.826 against a sampled 0.823, sphere tetrahedra
  0.617 against 0.638; 6000 tetrahedra checked in 0.08 s.
- The block on the sphere cannot be untangled by moving free nodes: its curved face nodes
  are pinned to the sphere and the straight corner cells are nearly flat.

### Verified
- `pytest` +9 (`tests/python/test_repair.py`: coefficients of affine elements, pushed /
  folded quads, curved rings, sampled kinds, the repair, the hopeless case); lint, ternary
  check, API docs, launcher check.

## Phase 12b — end-to-end tests of feature combinations (2026-10-17) — requested by the user, committed when its checks passed

- `tests/python/test_combinations.py`: one checker for every model - the merged mesh is
  consistent (every point used by a cell, every patch line / face a side of a cell, the
  patches cover exactly the mesh boundary, no side shared by more than two cells), every
  export read back by its reference reader with the right counts (VTKHDF through VTK, Gmsh
  through meshio by element type and dimension, CGNS through `cgnscheck`), curved high-order
  meshes at orders 2 (and 3 in 2D) valid with boundary elements matching the patches and
  written to Gmsh / CGNS and read back, the journal replaying to the same model, the project
  round-tripping to the same assembly.
- Nine models: C-grid template + 12 wall layers + smoothing group + patches; the periodic
  channel with zones; the hybrid ring; a cylinder with 2D inflation layers; the circle O-grid
  smoothed as a group; a sketch; the pipe template (block faces without domains patched as
  block faces); the extruded bump channel with 3D wall layers; the sphere in a box with
  hexahedral layers. 20 s in all.
- Found: nothing new in the code (the fixtures needed the 3D patch-member and face-table
  forms). The connector-row export bug of 11g is now covered by the first and fifth cases.

## Phase 12a — GUI reliability (2026-10-17) — requested by the user, committed when its checks passed

### Built
- `tests/python/test_gui_events.py`: a `Live` harness that gives the off-screen window a real
  interactor on its widget (as the on-screen window has) so `QTest` mouse and key events on
  the viewport run through VTK's interactor, the window's observers and the camera style;
  17 tests: the chain itself, click / Shift-click selection, the line tool with its preview,
  a drag in a drawing tool (zooms, adds no point), Esc and Return in the viewport, handle
  drags, freehand strokes, the sketch tool, tool switches, the boundary dock's buttons, the
  Delete key, model-tree clicks, wheel zoom and F, Shift-drag panning, a property edit typed
  into the panel and applied.
- Found and changed: a click in the viewport picks the connector lying on a curve, never the
  curve, so its control points could not be dragged from the viewport (only after selecting
  the curve in the tree). A selected connector now offers the handles of its curves too
  (`handles_for`); a drag moves the curve and the connector and domains follow.
- View -> "Show high-order elements…" (`app/high_order_view.py`, `Scene.show_high_order`,
  `MainWindow.show_high_order`): the element edges sampled along their curves and the
  high-order nodes drawn over the mesh (2D: all edges; 3D: the patches' faces), refreshed on
  every model change, with the count, node count and the curved elements' scaled Jacobian in
  the status bar; `curved=False` shows the straight-sided elements.

### Verified
- `pytest` +18; no other interaction bug was found by the chain tests (the earlier
  world-space tests remain). Lint, ternary check, API docs, launcher check.

## GUI navigation, docks and a launcher app (2026-10-16) — requested by the user

- Zoom to a window (`app/zoom_box.py`): a left-drag (a three-finger drag on a trackpad) draws a
  rectangle that follows the pointer and zooms onto it on release (box centre to the view
  centre, zoom by the box's size); a click still selects, handle drags and freehand strokes
  keep the button first. Panning moved to Shift + left-drag and middle-drag; View -> "Left-drag
  zooms to a window" (on by default) gives left-drag back to panning.
- The Properties, Smoother and Boundary forms scroll inside their docks: a connector's form
  needed 484 px, which kept the Model tree above it at 107 px.
- `packaging/macos/install_launcher.sh`: OpenGrid.app in /Applications that starts the
  development environment (current code, no rebuild; log in ~/Library/Logs/OpenGrid.log),
  with the O-grid icon; installed.
- Tests `tests/python/test_gui_navigation.py` (zoom maths against the view's world corners, the
  rectangle and the click threshold, docks).
- Fixed (reported by the user: the zoom window stayed and could not be cancelled): the camera
  style saw the press, started a pan and grabbed the interactor's focus (VTK's trackball
  style does in OnLeftButtonDown), after which only it got the moves and the release, so the
  box never finished and then followed the pointer. The zoom observer now takes the press
  (as handle drags do); Esc and tool changes cancel a box; a move without the button ends it
  (a missed release). `test_zoom_window_through_the_event_chain` drives a real interactor on
  the off-screen view (it fails without the fix; the earlier tests did not reach the chain).
- .ogp projects open from the Finder (double-click; also while OpenGrid runs): the launcher
  declares the type `org.opengrid.project` and is made its default handler
  (`packaging/macos/set_default_handler.py`; the stale self-contained bundle in packaging/dist
  is unregistered); `python -m opengrid project.ogp` opens it too (`__main__.py`: file-open
  events and arguments, unsaved changes confirmed). Round icon (transparent outside the disk).
- Staying current: the window title shows the commit it runs ("+": uncommitted changes); a
  status-bar notice appears when OpenGrid's code changes on disk while it runs (checked every
  30 s, Python files and extension modules); the bundle records what it was made from
  (`Contents/Resources/source.txt`), and `tools/check_launcher.py` / `tests/python/
  test_launcher.py` fail when the installer or icon script changed since. The app's log goes
  to ~/Library/Logs/OpenGrid.log. (Found on the way: an app started by the Finder can sit in a
  folder it may not read, and git then fails even with -C: git now runs in the checkout.)

## Phase 11g — curved boundary elements (2026-10-16) — requested by the user, committed when its checks passed

### Built
- `model/curving.py` + `Session.high_order(order, curved=True)` (now the default;
  `export_high_order(..., curved=)`, GUI choice): the boundary's high-order nodes moved onto
  the geometry - 2D: the connectors' curves; 3D: border edges of each face domain onto the
  nearest of its connectors' curves, other nodes onto its surface, or the degree-p
  interpolation of a structured face's grid without one; planar faces and periodic patches
  stay - then the displacement extended into the mesh and every changed element checked
  (folds: ValidationError with element, kind and location). Report `mesh.curving`: moved and
  solved nodes, iterations, curved elements, their scaled Jacobian curved and straight-sided.
- `core/src/highorder/curving.cpp` (`og::extend_displacement`, binding
  `_core.extend_displacement`): weighted graph Laplacian on the element lattices' neighbour
  edges (1 / length^3), Dirichlet boundary, region of normalised distance <= 1 from the moved
  nodes (3 boundary-element lengths, multi-source Dijkstra); Eigen sparse Cholesky up to 20k
  unknowns, conjugate gradients with incomplete Cholesky above. Eigen (already in the
  environment) is now linked by the core (private, header-only).
- `model/shape.py`: Lagrange shape functions and gradients of all seven kinds at any order
  (tensor products, Silvester's simplex formula, a rational pyramid basis in collapsed
  coordinates), Jacobians at the order-(p + 2) lattice.
- `ElementMesh.edge_nodes`, cached key lookups; assemblies record every boundary connector /
  face source; without patches `high_order` works (no boundary elements).
- Fixed: the 2D assembly took a connector's nodes by merging the connector's own nodes with the
  grids by coordinates, which fails for group-smoothed domains (their boundary nodes slide):
  326 of the 392 patch lines of the C-grid template's export pointed at stray points (Gmsh,
  CGNS, VTKHDF exports of such models). Connector nodes now come from the grid row they bound.
- Demo `examples/11g_curved_elements.py`; tests `core/tests/test_curving.cpp`,
  `tests/python/test_curved_elements.py`.

### Results
| Cylinder ring, element edges vs the circle | 16 wall edges | 32 | rate |
|---|---|---|---|
| straight-sided, any p | 1.9e-2 | 4.8e-3 | 2.0 |
| p = 2 | 4.6e-5 | 2.9e-6 | 4.0 |
| p = 3 | 1.1e-4 | 6.9e-6 | 4.0 |
| p = 4 | 3.4e-7 | 5.3e-9 | 6.0 |
| p = 5 | 1.3e-6 | 1.9e-8 | 6.0 |

- Thin layers: 16 wall edges over a 1e-5 first cell (wall nodes moved 1.9e-2 = 1921 first
  cells), order 4: valid, scaled Jacobian 0.823 (straight-sided 0.915). The weight exponent:
  1 / length^2 lost 0.20 there, ^3 0.10, ^4 0.08 (but more on the C-grid): 3.
- NACA 0012 C-grid (1e-5 wall, 24 layers), order 4 in 0.55 s: wall nodes 2.2e-4 -> 3e-17 off
  the airfoil, quality unchanged (0.480, the far field's).
- Sphere in a box, order 3 (1.1 s each): hexahedral layers 8.8e-3 -> 2.1e-5 off the sphere
  (scaled Jacobian 0.712 -> 0.623); prism layers on triangulated faces 1.3e-2 -> 8.5e-5
  (0.966 -> 0.638: a tetrahedron on top of the layer stack, its bottom face curved).
- O-grid ball (32000 hexahedra), order 4: 8.6e-4 -> 1.7e-11, 567k unknowns in 92 CG
  iterations, 10 s (37 s before caching the key lookups and BLAS Jacobians).
- CAD bracket (STEP), fluid volume, order 2 (116k elements, 7.8 s): nodes on its NURBS faces
  4.7e-3 -> 7e-16.
- Pipe template (block faces without surfaces, grid interpolation), orders 2-4: wall nodes
  3.3e-4 -> 7.7e-6 / 1.2e-6 / 2.6e-7; quality unchanged; planar end caps untouched.

### Verified
- `ctest` +3 -> 154; `pytest` +17 (and the 11f tests on curved defaults); strict CI build,
  clang-tidy on the new C++ (`-DEIGEN_DONT_VECTORIZE`, noted in CLAUDE.md), ruff, clang-format,
  ternary check, API docs; demos 01-11g run.

## Phase 11f — straight-sided high-order elements of every cell type (2026-10-16) — requested by the user, committed when its checks passed

### Built
- `io/lagrange.py`: the equispaced node lists of order-p lines, triangles, quadrilaterals,
  tetrahedra, pyramids, prisms and hexahedra as exact fractions in one canonical reference
  element per kind, in the node orders of Gmsh, VTK and CGNS, each generated from a table
  of edges (with directions), faces (with vertex orders) and interior modes; shape-function
  weights of the linear element (the pyramid's collapsed map), permutations between formats,
  Gmsh element type numbers (orders 1-6), SIDS names (orders 1-4), VTK Lagrange types.
- `model/elements.py`: `build_elements(assembly, order)` turns every cell of the linear mesh
  into an order-p element with the linear cell's geometry (straight edges, planar triangles,
  bilinear quadrilateral faces); edge and face nodes merged exactly through integer keys,
  coordinates of a shared node computed once from its key; boundary lines / triangles /
  quadrilaterals per patch; periodic node correspondences extended to the high-order nodes.
  `ElementMesh`: `points` (the linear vertices first), `cells` / `faces` per kind (Gmsh
  order), `connectivity(kind, fmt)`, `element_nodes(kind, "gll" | "gauss")` for
  quadrilaterals and hexahedra, `counts()`.
- Writers: Gmsh MSH 4.1 high-order types with the physical groups, entities and `$Periodic`
  links of the linear export (`write_gmsh_elements`; one writer for both); CGNS
  (`io/cgns_elements.py` + C++ `write_cgns_elements`: one unstructured zone per zone, SIDS
  sections up to order 4, boundary sections with BCs per patch (EdgeCenter / FaceCenter),
  vertex connectivity between zones, translational periodic links, 2D bases with X, Y);
  VTKHDF Lagrange cells (pyramids linear: VTK has no Lagrange pyramid); HDF5 element nodes
  (equispaced, or GLL / Gauss for quadrilaterals and hexahedra).
- `Session.high_order(order)` / `export_high_order(path, order, nodes)` (by extension:
  `.msh`, `.cgns`, `.vtkhdf`, `.h5`); the 11e curved elements are now
  `curved_high_order` / `export_curved_high_order`. GUI: File -> Export high-order elements
  (Gmsh, CGNS, VTKHDF, HDF5 node sets) and Export curved high-order elements (structured).
- Demo `examples/11f_high_order_elements.py`; tests `tests/python/test_high_order_elements.py`
  with reference data in `tests/python/data/` (generated with Gmsh outside the repository).

### Results
- Node orders: Gmsh 42 of 42 (7 kinds x orders 1-6), VTK 36 of 36 (no Lagrange pyramid in
  VTK 9.7), CGNS 26 of 26 written by Gmsh's CGNS writer, and all 46 SIDS node descriptions
  (TETRA_35, PYRA_55, PENTA_40, PENTA_75). Findings on the way: CGNS quadrilateral faces of
  order 4 list their interior ring by ring (not as a nested QUAD_9), HEXA_125 and PENTA_75
  interiors are layer by layer, Gmsh's prism interior is column by column with the two ends
  first, VTK's wedge orders its triangle faces by rows (its triangles recursively); VTK's CGNS
  reader puts TETRA_35's interior node at vertex 0 (its table is wrong).
- Read-back: Gmsh 4.15 interpolates every element of our files (all kinds, orders 2-4, prisms
  order 2: Gmsh cannot evaluate higher-order prisms) to within 1e-14 of the straight-sided
  map; `cgnscheck` accepts the CGNS files; VTK evaluates the Lagrange cells to 1e-12.
- 2D hybrid cylinder (1600 quadrilaterals + 1484 triangles): order 4 in 24 ms, 37752 nodes.
  The straight-sided wall is the linear polygon at every order (edges off the cylinder by the
  sagitta 7.7e-4 of an 80-edge circle); the curved order-4 elements of the structured ring
  follow it to 1e-8.
- 3D body in a box (27090 tetrahedra, 1010 pyramids, 820 prisms, 5650 hexahedra): order 4,
  706630 nodes in 1.6 s; order 3 written as Gmsh 1.4 s, CGNS 0.9 s, VTKHDF 0.8 s.

### Verified
- `pytest` +28 -> 392 passed + 1 xfailed (node orders against the three references, meshes
  of every kind, every export read back element by element, periodic and multi-zone CGNS,
  GUI); `ctest` 151 (strict CI build); clang-tidy on the changed C++, ruff, clang-format,
  ternary check, API docs; demos 01-11f run.

## Phase 11e — high-order elements (2026-10-16) — requested by the user, committed when its checks passed

### Built
- `core/src/highorder/nodes.cpp`: equispaced, Gauss-Lobatto-Legendre and Gauss-Legendre
  points and weights (Newton on Legendre polynomials), barycentric interpolation and
  differentiation matrices.
- `core/src/highorder/elements.cpp`: order-p quadrilaterals / hexahedra from a structured grid
  whose interval counts are multiples of p (p x p (x p) cells per element; the grid's nodes
  are the equispaced node set), other node sets by tensor-product interpolation (sum
  factorisation), Jacobian min / max at (p + 2)^d GLL points per element.
- `model/highorder.py` + `Session.high_order(names, order, nodes, snap)`: domains (2D) or
  blocks (3D); GLL nodes snapped onto the curves of boundary sides (not interfaces) and the
  surfaces of block faces (face interiors). `Session.export_high_order`: VTKHDF
  (VTK_LAGRANGE_QUADRILATERAL / HEXAHEDRON, equispaced, shared nodes merged; point order as
  vtkHigherOrder*::PointIndexFromIJK, checked for p = 1..8) or `.h5` (any node set, merged
  numbering). GUI: File -> Export high-order elements.
- Sketches and templates: `interval_multiple` rounds every edge family's interval count (and
  the extrusion's) up to a multiple.
- Demo `examples/11e_high_order.py`; tests `core/tests/test_highorder.cpp`,
  `tests/python/test_high_order.py`.

### Results
| Quarter annulus, inner arc | 4 x 4 elements | 8 x 8 | rate |
|---|---|---|---|
| p = 1 | 1.9e-2 | 4.8e-3 | 2.0 |
| p = 2 | 4.6e-5 | 2.9e-6 | 4.0 |
| p = 3 | 1.2e-5 | 7.6e-7 | 4.0 |
| p = 4 | 3.8e-8 | 5.9e-10 | 6.0 |
| p = 6 | 1.5e-11 | 5.9e-14 | 8.0 |

- (h^(p+1), one order more for even p on a circle.) NACA 0012 C-grid template, 1e-5 wall,
  24 wall layers, order 4: 828 curved quadrilaterals, scaled Jacobian 0.064 - 0.981 (the
  smallest at the trailing edge), GLL nodes on the airfoil moved at most 3.6e-5 chords by the
  snapping; pipe (3D), order 2: 1536 hexahedra, scaled Jacobian 0.164 - 0.999.

### Verified
- `ctest` +4 -> 151; `pytest` +6 (VTK reads the Lagrange cells back; one cell evaluated by VTK
  equals our mapping to 1e-12); strict CI build, clang-tidy on the new C++, ruff, format,
  ternary check, API docs.

## Phase 11d — GPU backend, experimental (2026-10-16) — requested by the user, committed when its checks passed

### Built
- `core/include/opengrid/core/gpu/gpu.hpp`, `core/src/gpu/` (optional `OG_WITH_OPENCL`, the
  macOS framework or an ICD loader): device selection (the GPU with the most compute units),
  kernels built once, a float-float self-test (a driver that reassociates fails it), modes
  off / auto (>= 200k nodes) / on; `opengrid.gpu_info()`, `opengrid.set_gpu()`,
  `opengrid.gpu_mode()`; the GUI logs the compute devices at start.
- `Multiblock3DSweeps`: the coefficient pass and the line relaxation of the 3D group smoother
  on the device (float-float positions, float32 delta form, `FP_CONTRACT OFF`); the device
  colours lines that share nodes (8 colours on the 1M-cell ball instead of 21) and double
  buffers positions; mapped buffers. Anderson mixing, sliding, irregular and wall updates stay
  on the CPU; GPU runs stop at max(tolerance, 1e-5).
- Demo `examples/11d_gpu.py`; tests `core/tests/test_gpu.cpp`, `tests/python/test_gpu.py`
  (skipped without a device).

### Results (M1 Max, 1M-cell ball)
| | CPU | GPU line sweeps |
|---|---|---|
| sweeps per iteration | 36 ms | 16.6 ms (first version with the CPU colours: 22 ms) |
| whole iteration | 64-66 ms | 54 ms |
| run (tolerance 1e-6 / 1e-5) | 21.5 s, 328 iterations | 15.6 s, 290 iterations |
| grid vs the CPU's | - | up to 7 % of the local cell size (stopped at 1e-5) |
| to the CPU's accuracy | 328 iterations | floor 1.3e-6 (~570 iterations): slower than the CPU |

- 10d wing stub (159k nodes, lines of ~9 nodes): GPU sweeps 4.1 ms per iteration = the CPU's
  4.0 ms (latency-bound: ~1200 short lines per launch).
- Decision: kept as an experimental option, off by default (it does not accelerate to the
  same accuracy). What it would take: Anderson mixing on the device, a line solver with more
  parallelism (parallel cyclic reduction, or structured per-block kernels with interface
  exchange - where the benchmark's 7-9x comes from), and faster CAD projections of sliding
  nodes (CPU-bound).

### Verified
- `ctest` +2 -> 147 (release and strict CI presets, OpenCL deprecation warnings silenced);
  `pytest` +2; clang-tidy on the new and changed C++, ruff, format, ternary check; GPU vs CPU
  agreement and repeatability tested on a sheared, wall-clustered box.

## Phase 11d (assessment) — GPU precision and speed (2026-10-16)

Tools in `tools/gpu/` (not part of the package; `opencl_bench.py` needs pyopencl):
- `precision_study.py`: an O-grid with a 1e-6 first cell at radius 1 (outer radius 20),
  inverse control functions so that the exact polar grid is the discrete solution, a sheared
  start, radial line relaxation in delta form, as the CPU smoother does:

| Variant | Error at the wall row (local cell heights) | First cell (exact 1e-6) | Residual floor |
|---|---|---|---|
| float64 | 1.6e-8 | exact | (6000 iterations: 8e-9) |
| float32 arithmetic, positions in double (two float32 on a GPU) | 2.0e-7 | 0.9999998e-6 - 1e-6 | 2.6e-6 |
| float32 throughout | 5.7e-2 | 0.954e-6 - 1.029e-6 | 2.3e-5 |

- `opencl_bench.py` / `cpu_bench.cpp`: 3D Winslow sweeps on a 128 x 128 x 64 block, the same
  kernels on both sides (the CPU line kernel eliminates plane by plane, i innermost):

| Kernel | M1 Max GPU (OpenCL 1.2) | 10 CPU threads, float64 |
|---|---|---|
| point Jacobi, float32 | 0.29 ns/node | - |
| point Jacobi, float-float positions | 0.45 ns/node | 4.1 ns/node (9x) |
| line relaxation (Thomas along k), float-float | 0.71 ns/node | 5.2 ns/node (7.4x) |

  GPU float-float and CPU float64 agree to 2.8e-7 of the local cell height after 200 line
  sweeps; on a 40k-node block the GPU loses its lead (launch overhead, 0.2 ms per sweep). The
  current general multi-block smoother spends ~36 ns per node and iteration on its sweeps and
  coefficients (indirect neighbour access, per-node frames).
- Decision: the backend pays for large 3D groups. Plan: OpenCL (Apple, NVIDIA, AMD, Intel;
  OpenCL is deprecated but working on macOS), positions as float-float, float32 arithmetic in
  delta form, coefficients + line relaxation + Anderson mixing on the GPU, sliding,
  irregular and wall-control updates on the CPU through shared buffers, a float64 CPU finish
  below the float32 residual floor; runtime detection, a self-test against the CPU and a size
  threshold, CPU fallback. Expected 2-7x end to end (sliding projections stay on the CPU).

## Phase 11c — coarse-to-fine smoothing (2026-10-16) — assessed, not built

- Prototype (scratch, not in the code): the coarse grid is every other node of each block of
  the merged layout (all interval counts even), smoothed with the same settings; its
  displacement is interpolated (trilinear in index space) onto the fine starting grid.
- 1M-cell ball (10b, 65 x 33 points per block): cold 328 iterations / 21.0 s; coarse 134k
  nodes, not converged in 2000 iterations (17.0 s, stall guard fired twice); fine 218
  iterations / 13.5 s with the cold run's absolute stopping level; total 30.5 s. The result
  differs from the cold one by 9.6e-3 (domain size 1): sliding nodes' limits are relative to
  their starting positions. 33-point ball: cold and sequenced both run to 2000 iterations at a
  residual floor (relative 1.5e-6).
- Decision: not built (target 2x). The iteration counts are set by a residual floor and
  the slide limits rather than by smooth error modes; see OPEN_ISSUES.

## Phase 11b — CPU performance (2026-10-16) — requested by the user, committed when its checks passed

### Built
- `parallel_for`: dynamic schedule in chunks (about 16 per thread) instead of a static
  partition (69 % of the thread time was idle on the wing stub).
- `util/scratch.hpp`: stack scratch for the NURBS curve and surface basis tables (no heap
  allocation per evaluation; half the busy time had been malloc/free); `Nurbs::eval` position
  only; `Surface::values` (position + derivatives in one NURBS evaluation) used by
  `project_from`.
- 3D group smoother: state packing and point updates in parallel.
- Demo `examples/11b_performance.py`; tests: NURBS curves (high degree, eval vs d1) and
  surfaces (values vs separate calls, high degree).
- Tried and dropped: a closed-form projection onto planes (3 % faster, but results changed at
  1e-7 and the iteration count went 908 -> 954).

### Results (10 threads, M1 Max)
| Case | Before | After |
|---|---|---|
| 10d wing stub group smoothing | 39.3 s | 9.2 s (one thread 43.7 s) |
| 10d wing stub sketch passes | 27.5 s | 7.5 s |
| 9c 2D elliptic 500 x 200 | 1.10 s | 0.68 s |
| 9c 3D elliptic block 81 x 41 x 41 | 0.96 s | 0.58 s |
| 10b 1M-cell ball | 26.2 s | 22.3 s |

- Every grid is bit-identical to 11a (wing stub digest; demos unchanged except timings).

### Verified
- `ctest` +2 -> 145 (release and strict CI presets); `pytest` 356 passed, 1 xfailed;
  clang-tidy on the changed C++ files, ruff, format, ternary check; all demos run.

## Phase 11a — wall layers (2026-10-16) — requested by the user, committed when its checks passed

### Built
- `structured/wall_layers.cpp` (`march_wall_layers`): the first N layers off a wall marched
  along the front normals (facet-normal bisectors, 1/10 neighbour blending, symmetry planes
  kept), heights spacing * growth^k, folds / overlaps / short lines reported; the grid beyond
  shifted along its lines. Used by the 2D and 3D group smoothers before their line setup; layer
  nodes are node class 4 (`layer`, kept fixed); the wall control moves to each line's last layer
  node (next height spacing * growth^N).
- Python: `layers` / `growth` on group walls, sketch walls (also `{edge: {"spacing", "layers",
  ...}}`) and templates (all walls; a C-grid's wake cut then takes the wall spacing). GUI wall
  text `name [spacing [layers [growth]]]` (groups and sketches); group status counts the layer
  nodes.
- Anderson guard: a non-finite residual or 200 iterations without a 10 % drop also restart
  from the best iterate.
- Demo `examples/11a_wall_layers.py`; tests `core/tests/test_wall_layers.cpp`,
  `tests/python/test_wall_layers.py`.

### Results (first cell 1e-5 chords, growth 1.2; deviation of grid lines from the wall normal)
| Case | Layer 8 | Layer 16 | Layer 29 | Iterations | min SJ |
|---|---|---|---|---|---|
| 2D C-grid, wall control only | 2.4° (max 12°) | 10.2° | 18.8° | 630 | 0.576 |
| 2D C-grid, 30 layers | 0.09° (max 1.1°) | 0.09° | 0.07° | 142 | 0.468 |
| 3D tapered + swept wing, wall control only | 3.2° | 11.7° | 20.8° | not converged (4000) | 0.434 |
| 3D tapered + swept wing, 30 layers | 0.07° | 0.07° | 0.05° | 1770 (15.7 s) | 0.407 |

- First cells exactly the requested height; the cell beyond the layers is controlled
  (2.5-2.8° mean deviation at layer 32 in the outer grid).
- 1e-6 with 30 layers: 0.24° (±2-node tangents) / converged in 147 iterations; 40 layers:
  converged after the stall guard (min SJ 0.29: only 21 lines left for the outer grid).

### Verified
- `ctest` +4 (curved and sheared 2D/3D marching, one-sided ends, reports) -> 143; `pytest` +4
  -> 356 passed, 1 xfailed; strict CI build, clang-tidy on the new and changed C++, ruff,
  format, ternary check, API docs; all demos 01-11a run.

## Phase 10e follow-up — extruded templates, symmetry planes, wall control fixes (2026-10-16) — approved (pre-authorised)

### Built
- Sketch `extrude` (2D sketches): the smoothed 2D grid swept into owned extruded blocks
  (`vector`, `points`, `distribution` with spacings as lengths). The extruded templates (pipe,
  cylinder_in_channel, bump_channel, wing_in_box) use it; `three_d_blocks=True` keeps six-face
  3D sketches.
- `wing_in_box` `taper` / `sweep`: a 3D sketch whose tip section is a similarity of the root
  (about the leading edge), the wing on ruled surfaces, the leading-edge span edge on a line.
- Symmetry planes (`symmetry` on 3D groups and sketches; C++ `Multiblock3DInput::symmetric`):
  planar faces whose nodes move as interior nodes mirrored across the plane (no slide limit),
  with edge rules (two planes; a plane and another surface) and wall lines through them.
  Templates set their end planes; the 10d wing stub its side walls. Sketch face corners in
  one surface or plane now slide on it (`joints="slide"`).
- Fix, 2D internal walls: the wall-normal direction of a fixed interior wall node comes from
  its wall neighbours (a turned local frame applied the control along the C-grid template's
  wake cut: lines at ~45°, kink 120-130°).
- Fix, Anderson divergence: `AndersonAccelerator::guarded_step` (all four elliptic smoothers):
  restart from the best iterate when the residual exceeds 100x its lowest value, plain
  iteration after two restarts, with a warning. The C-grid at wall spacings 1e-5 / 1e-6 had
  diverged and silently fallen back to its TFI grid.
- C-grid default `wake_points` stays 41 (61 was better only while the wake wall was broken).
- Demos: 10d adds the wing stub as a 2D sketch extruded on the CAD root curves
  (`10d_wing_extruded.ogp`); 10e adds the tapered + swept wing and the C-grid wall control
  at 1e-3 / 1e-5 / 1e-6.

### Results
| Case | Before | Now |
|---|---|---|
| Wing in box template (extruded) | 0.249 (missed) | 0.619 (met) |
| Wing in box, six-face 3D / taper 0.6 / sweep 20° / both | 0.24 / - / - / - | 0.48 / 0.55 / 0.48 / 0.60 |
| C-grid template | 0.430, wake kink 130° | 0.619, kink 22° |
| 10d wing stub (3D sketch), wake blocks | 0.12 | 0.616 (0.653 elsewhere) |
| 10d wing stub, 2D sketch extruded | - | 0.679, nodes on the CAD faces within 4.8e-15 |
| C-grid wall spacing 1e-5 / 1e-6 | diverged (TFI returned) | converged in 600 / 3125 iterations |

- Wall control at convergence (C-grid, 1e-3 / 1e-5 / 1e-6): first cell exactly the requested
  height on the body and 2x on the wake cut, first grid line within 0.005° of the normal.
  Further layers (measured, see OPEN_ISSUES): 1e-3 - 4° at layer 4, 9° at layer 8; 1e-5 -
  1e-6 - below 0.5° to layer 4, 1-2° at layer 8 (means).
- Earlier demos: grids unchanged (projects and journals gain `"symmetry": []`).

### Verified
- `ctest` +1 (Anderson guard) -> 139 (release and strict CI presets); `pytest` +3 (C-grid
  wall and wake orthogonality at 1e-3 and 1e-5, extruded / tapered + swept wing) -> 352 passed,
  1 xfailed; clang-tidy on the changed C++ files, ruff, format, ternary check, API docs;
  all demos 01-10e run.

## Phase 10e — templates and automatic topology (2026-10-16) — approved (pre-authorised)

### Built
- `model/templates.py` + `Session.create_template`: 2D ogrid, butterfly, cgrid, hgrid,
  body_in_channel; 3D pipe, cylinder_in_channel, bump_channel, wing_in_box (extruded) and
  sphere_in_box. Sketch: edges across closed-curve seams, `slide_limit`.
- `model/crossfield.py` + `Session.auto_topology` / `cross_field`: cross field, singularities,
  index theorem, separatrices, block decomposition into a sketch.
- GUI: Topology -> Template: O-grid / C-grid / butterfly from the selected curves, Automatic
  topology (2D) from the selected curves.
- Demo `examples/10e_templates_auto.py`; tests `tests/python/test_templates.py`.
- Test fix: `test_cancel_smoothing` cancels at once (the small domain could reach round-off
  before a delayed cancel: a timing race, failing about every second full run).

### Results
| Template | Blocks | min scaled Jacobian | Target 0.4 |
|---|---|---|---|
| O-grid (ellipse) | 4 | 0.953 | met |
| Butterfly (disk) | 5 | 0.864 | met |
| C-grid (NACA 0012) | 12 | 0.430 | met |
| H-grid (bump channel) | 3 | 0.989 | met |
| Body in channel (cylinder) | 12 | 0.782 | met |
| Pipe (3D) | 5 | 0.864 | met |
| Cylinder in channel (3D) | 12 | 0.673 | met |
| Bump channel (3D) | 3 | 0.993 | met |
| Wing in box, C-H (3D) | 12 | 0.249 | **missed** |
| Sphere in box (3D) | 32 | 0.525 | met |

| Automatic topology | Singularities | Index sum = Euler | Blocks (min SJ) |
|---|---|---|---|
| Airfoil in a box | 2 x -1/4 | 0 = 0 | 6 (0.995) |
| Two-element airfoil | 4 x -1/4 | -1 = -1 | **failed**: a face in the slot has 5 corners |
| Plate with three holes | 12 x -1/4 | -2 = -2 | 32 (0.706) |
| Bracket section | 4 x -1/4, 4 x +1/4 | 0 = 0 | 24 (0.131) |

### Verified
- `pytest` +21 (`test_templates.py`; the two-element case is a strict xfail) -> 349 passed;
  ruff, format, ternary check, API docs; demos 10d and 10e run (no C++ changes).

## Phase 10d — sketch topology (2026-10-16) — approved (pre-authorised)

### Built
- `model/sketch.py`: the `sketch` entity (vertices + blocks, 2D and 3D), edge families and
  dimension propagation (conflicts reported by edge), geometry assignment and placement,
  generation of owned connectors / domains / blocks / smoothing group, vertex passes (2D
  joints on curves, 3D vertices on planes and surfaces), suggestions, ownership checks.
- Session: `create_sketch`, `edit_sketch`, `regenerate_sketch`, `sketch_report`,
  `sketch_suggestions`, `accept_sketch_suggestions`, `assign_sketch_geometry`;
  `move_handle` on sketch vertices; `modify` of a sketch regenerates it.
- Smoothing groups: `joints="slide"`, `pinned`, `slide_limit`, internal walls (2D and 3D,
  C++), `absolute_tolerance` (C++ option); one plane patch for coplanar sliding faces; wall
  faces projected onto their surface.
- Fixes: `Store.rollback` restores ids and name counters; the smoother dock replaces a
  running job instead of mixing two (results tagged by job).
- GUI: Sketch block tool (four clicks, snapping to sketch vertices), sketch wireframe and
  vertex handles (drag with merged-actor preview), Topology -> Assign selected geometry /
  Regenerate sketch / Sketch report, sketch rows in the property panel (report, dimensions,
  walls, relax), sketches in the model tree.
- Demo `examples/10d_sketch_topology.py`: airfoil (2D) and wing stub (3D) around STEP models
  written and imported through OCCT, sphere in a box.

### Results
| Check | Result | Target |
|---|---|---|
| Airfoil, C-topology (12 blocks): min scaled Jacobian | 0.679 | >= 0.5 met |
| Airfoil: floating vertices moved by 20 % of the block size, 3 trials | grid changes <= 2.1e-4 local spacings | < 1e-2 met |
| Airfoil: nodes on the assigned CAD curves | 1.8e-15 | 1e-9 met |
| Airfoil: generation (3 edits, 11 vertex passes) | 2.7 s | — |
| Airfoil: drag preview of a sketch vertex | 95 ms (11 frames/s) | — |
| Sphere in a box (32 blocks): min scaled Jacobian | 0.688 | >= 0.5 met |
| Sphere: inner vertices moved by 20 %, 2 trials | <= 1.2e-8 local spacings | — |
| Sphere: nodes on the sphere | 2.2e-16 | 1e-9 met |
| Wing stub (12 blocks): min scaled Jacobian | 0.12 (0.68 outside the 4 wake blocks) | — |
| Undo per edit, journal replay, .ogp round trip | identical grids | met |

### Verified
- `ctest` +3 (internal walls 2D and 3D, absolute tolerance) -> 138; `pytest` +9
  (`test_sketch.py`) -> 329; CI build, clang-tidy on the changed C++, ruff, clang-format,
  ternary check, API docs; all demos run. Changed outputs: groups in projects and journals
  now carry `joints`, `pinned`, `slide_limit`.
- Not verified: dragging sketch vertices with a real mouse; 3D vertices have no axis handles
  (they move in the view plane).

## Phase 10c — live mouse editing, freehand curves, curve networks (2026-10-16) — approved (pre-authorised)

### Built
- `model/fitting.py` + `Session.create_freehand`: strokes -> NURBS curves (corners kept).
- `model/network.py` + `Session.domains_from_curves`: crossing curves -> conforming structured
  domains (or unstructured for faces without four corners), one undo step.
- `Session.move_handle` (topology edit): move a point / curve defining point; crossings of the
  moved curve with other curves are re-solved and the connectors follow.
- GUI: drag handles of selected entities with live preview (`app/handles.py`), Freehand and
  NURBS tools, Create -> Domains from curves; drags rerun the affected smoothing on release.
- Demo `examples/10c_freehand_live.py`.

### Results
| Check | Result |
|---|---|
| Quarter-circle stroke, 0.5 % tolerance | 4 control points (<= 12 required) |
| NACA / noisy S strokes | within 0.5 % of the stroke length |
| 3 x 2 lattice of crossing strokes | 6 structured domains, 7 conforming interfaces |
| Five-sided face | reported by name (or unstructured on request) |
| Drag preview, spline with a 200 x 100 domain | ~12 ms (>= 10 frames/s required) |
| Drag preview, six-domain freehand network | 83 ms (12 frames/s) |
| Drag = one undo step, journal replays | yes (`og.move_handle(...)`) |

### Verified
- `pytest` +10 (`test_live_editing.py`) -> 320; ruff, ternary check, API docs; all demos run.
- Not verified: real mouse and trackpad dragging (offscreen tests drive the world-space
  handlers; the VTK observer path ran only in the screenshot session).

## Phase 10b — multi-block smoothing groups, 3D (2026-10-15) — approved (pre-authorised)

### Built
- C++ `structured/multiblock_elliptic3d.cpp` + binding `_core.smooth_multiblock_3d`: one
  elliptic smoothing over the hexahedra of several blocks (floating interfaces, edges and
  corners; irregular edges; sliding on surfaces and curves; wall control along grid lines;
  line relaxation across blocks; thread-independent).
- Python `model/groups3d.py`: groups of six-face blocks (`create_smoothing_group(blocks=...)`,
  face walls and sliding faces), face roles, planes for planar sliding faces, routing of block
  grids and face displays, block topology edits keep groups; project section `groups`.
- GUI: Smoothing group from selected blocks (selected faces become walls), blocks in the
  property panel, a selected block runs its group in the dock.
- Anderson accelerator: fused passes (bitwise identical results, ~2x faster on big states).
- Group default Anderson depth 20 (was 10): 2D 68k O-grid 575 -> 436 iterations, 1M-cell
  ball ~800 -> 320.
- Demo `examples/10b_multiblock_3d.py`: sphere in a box (CGNS export), the solid ball, a
  1M-cell timing.

### Results
| Case | Per block | Smoothing group |
|---|---|---|
| Split sphere block (vs unsplit) | — | equal to 1e-9 (point sets) |
| Ball (core + 6 shells): min scaled Jacobian | 0.577 | 0.750 |
| Ball: interface kink | — | 0.52° (away from irregular edges and face edges) |
| Ball: nodes on the sphere | — | 2.2e-16 |
| Sphere in a box: first-cell height on the sphere | — | 0.01000 .. 0.01000 |
| Sphere in a box: min scaled Jacobian | 0.577 | 0.327 (box corners: three blocks in one 90° corner) |
| Sphere in a box: interface kink | — | 31° with wall control (inherited from the cubed sphere's own 36° grid turns at its seams); 1.5° without |
| 1M-cell ball | — | 328 iterations, 27 s (target 10 s **missed**) |

### Verified
- `ctest` +2 → 135 (split block = single block, bad input); `pytest` +5 (`test_groups3d.py`).
- Strict CI build, clang-tidy, ruff, clang-format, ternary check, API docs; all demos run.

## Phase 10a — multi-block smoothing groups, 2D (2026-10-15) — approved (pre-authorised)

### Built
- C++ `structured/multiblock_elliptic.cpp`: one elliptic smoothing over the quads of several
  planar structured domains — floating interfaces and block corners, grid lines continued
  across domains and relaxed implicitly (greedily coloured, parallel, thread-independent),
  irregular points, sliding boundary nodes on curve chains (`CurveChain::project_param`; the
  curve's local Newton projection is now public), wall control carried along grid lines,
  Anderson acceleration. Binding `_core.smooth_multiblock`.
- Python `model/groups.py`: entity `smoothing_group` (validation, journal names, project
  section "domains"), connector roles (interface / wall / sliding / fixed / outside), global
  numbering, `GroupEvaluator` (cache, prepare / run / finish for the GUI thread), interface
  kink metric. Session: `create_smoothing_group`, `smooth_group`, `group_status`,
  `display_connector_nodes`; members' `domain_grid` / `domain_status` / `smooth_domain` go
  through their group. Topology edits update groups.
- GUI: Create → Smoothing group, property panel editor, smoother dock runs groups, model tree
  entry, moved connectors drawn in place.
- Demo `examples/10a_multiblock_smoothing.py`: cylinder in a channel and circle O-grid.
- Fix in the 9b and 10a demo pictures: quality colours were assigned to the wrong cells
  (PyVista orders structured cells with the first array axis fastest); `09b_topology.png`
  regenerated.

### Results
| Case | Per-domain smoothing | Smoothing group |
|---|---|---|
| Split bump channel with wall (vs unsplit) | — | equal to 1e-11 |
| Circle O-grid: interface kink | 76° (34° away from corners) | 0.38° |
| Circle O-grid: min scaled Jacobian | 0.71 | 0.86 |
| Cylinder in channel: interface kink | 90° | 0.99° |
| Cylinder in channel: min scaled Jacobian | 0.45 | 0.36 (box corners, see OPEN_ISSUES) |
| Cylinder: first-cell height on the wall | 0.01 | 0.01000 .. 0.01000 |

- Sliding nodes stay on their curves to 1e-16; on a straight sliding wall (bump channel) the
  orthogonality error is 2.6° (5.5° with fixed nodes).
- 68k-node O-grid: 456 iterations, 1.5 s (10 threads; depth 10). Same digest on 1, 4 and 10
  threads.

### Verified
- `ctest` +3 (split = single block with wall, O-grid irregular corners + sliding + thread
  independence, input checks) → 133; `pytest` +10 (`test_groups.py`).
- Strict CI build, clang-tidy, ruff, clang-format, ternary check, API docs; all demos run.

## Phase 9c — performance and packaging (2026-10-13) — approved

### Built
- OpenMP (`OG_WITH_OPENMP`, `util/parallel.hpp`): deterministic `parallel_for` (static
  schedule, outputs per index, sums in index order). `opengrid.threads()` /
  `opengrid.set_threads(n)`; `llvm-openmp` added to `environment.yml`.
- Elliptic smoothing, 2D and 3D: zebra ordering relaxed in parallel, precomputed wall-source
  weights, Anderson acceleration (`util/anderson.hpp`; new setting `acceleration`, default 5,
  0 = plain; GUI spin box).
- Volume quality metrics in C++ (`quality/volume_quality.cpp`), computed once per volume mesh
  and reused by the floor and volume checks, `volume_quality` and the exports (the NumPy
  version is the test reference).
- 3D optimiser: whole-mesh scans in parallel; swap candidates from the bad tetrahedra only
  (same swaps, identical meshes on five reference cases). Netgen backend: bulk NumPy transfer
  and Netgen's thread pool. Assembled shells cached by face-mesh digest.
- Benchmark `examples/09c_performance.py` (report, JSON and chart in `examples/output`).
- Packaging: conda recipe (`packaging/conda`, rattler-build; recipe tests: imports, OpenMP, a
  smoke test meshing 2D + a Netgen volume + exports) and `packaging/macos/build_app.sh`
  (`OpenGrid.app` + DMG via conda-pack; binaries re-signed ad hoc after prefix rewriting).
- Fix found by the conda compiler (Clang 23): `geom/surface.cpp` lacked `#include <string>`.

### Results (Release, M1 Max, 10 threads)
| Case | Cells | Time | Target |
|---|---|---|---|
| 2D structured TFI 1001² (grid + quality + VTKHDF) | 1.00 M | 1.8 s | < 10 s met |
| 2D Delaunay | 1.23 M | 3.4 s | < 10 s met |
| 2D advancing front | 1.25 M | 4.1 s | < 10 s met |
| 2D hybrid (layers + triangles) | 1.54 M | 4.9 s | < 10 s met |
| 500 x 200 elliptic, wall control | 0.1 M | 1.1 s, converged (402 it.) | < 2 s met |
| 3D tetrahedra (Netgen) | 4.96 M | 176–186 s | < 120 s **missed** |

- The elliptic case on one thread: 3.2 s, identical grid. Without acceleration: not converged
  after 4000 iterations (residual 9e-3 of the first). Before this phase: 2000 iterations in
  11.7 s, not converged.
- Demo smoothers: NACA O-grid 891 → 129 iterations (0.8 → 0.1 s), C-grid 2110 → 249 (3.2 →
  0.3 s), spacing and orthogonality unchanged. Golden wing stub 8.6 → 3.0 s; pytest 194 → 63 s.
- 5M tetrahedra: Netgen ~157 s with its thread pool; our own part ~19 s (optimiser 13 s, quality and checks ~5 s). At 278k tetrahedra the whole
  volume went from 29 s to 8.8 s; before this phase 2.2M tetrahedra took 193 s for the volume,
  19 s for quality and 21 s for the VTKHDF export.
- Packages: `opengrid-0.0.1-py312_0.conda` (947 KiB) builds and passes its tests;
  `OpenGrid.app` 1.9 GB, DMG 722 MB; the app copied elsewhere runs the smoke test and the GUI
  with an empty environment; the phase screenshot was taken with the bundle's Python.

### Verified
- `ctest` +3 (parallel_for, Anderson, volume quality) → 130; `pytest` +9 → 296.
- Thread independence: identical digests for 1–10 threads (2D / 3D smoothing, volume meshes,
  quality). C++ volume quality = NumPy reference (1e-9; angles 1e-5 degrees).
- Strict CI build, clang-tidy on changed C++ files, ruff, clang-format, ternary check, API docs
  check; all demos 01–09b run. Changed outputs: projects / journals now carry
  `acceleration: 5`; smoothed grids differ at tolerance level (pictures 03, 05a, 05b, 07a).
