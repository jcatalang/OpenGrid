# OpenGrid — Progress log

## Phase 0 — Kickoff and scaffold (2026-10-11) — awaiting review
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

## Phase 1 — 2D geometry kernel + project I/O (2026-10-11) — awaiting review
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

## Phase 2 — Connectors and distributions (2026-10-11) — awaiting review
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

## Phase 3 — Structured 2D domains + first export (2026-10-11) — awaiting review
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

## Phase 4 — Patches, zones, solver-ready export (2026-10-11) — awaiting review
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

## Phase 5a — Unstructured 2D: constrained Delaunay + refinement (2026-10-11) — awaiting review
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
