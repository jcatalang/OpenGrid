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
