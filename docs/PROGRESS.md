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
