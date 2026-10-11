# OpenGrid — Architecture (Phase 0, v1)

This document records the design decisions taken before coding. `docs/SPEC.md` is the
requirements source; any deviation from it is justified here (see *Deviations and risks*).

## 1. Module layout

```
gauss/                         repository root (product name: OpenGrid)
├── CMakeLists.txt             top-level build (core lib, gtest, nanobind module)
├── CMakePresets.json          debug / release presets (pins the C++ compiler on macOS)
├── pyproject.toml             scikit-build-core packaging + ruff config + pytest config
├── environment.yml            primary dev environment (conda-forge)
├── core/                      C++20 library `og_core` — NO GUI / Python dependencies
│   ├── include/opengrid/core/ public headers, one subdirectory per module:
│   │     geom/     database entities: curves (later surfaces), evaluation, projection
│   │     topo/     connectors, domains, blocks, patches, zones, entity store, validators
│   │     dist/     1D distribution functions (uniform, geometric, tanh/sine Vinokur, ...)
│   │     structured/   TFI, elliptic (Winslow / Poisson-TTM / Steger-Sorenson), hyperbolic
│   │     unstructured/ CDT + refinement, advancing front, quad merge
│   │     inflation/    advancing-layer boundary layers (2D, then 3D)
│   │     smoothing/    Laplacian / angle-based / optimisation-based untangling
│   │     quality/      per-cell metrics
│   │     io/           VTKHDF, CGNS, Gmsh 4.1, Plot3D writers/readers
│   │     util/         robust predicates, progress/cancel token, small math
│   ├── src/               implementation, same subdivision
│   └── tests/             GoogleTest
├── bindings/              nanobind module `opengrid._core`
├── backends/              optional, separately switched plugins (OG_WITH_NETGEN, OG_WITH_MMG,
│                          OG_WITH_OCCT); each wraps an external LGPL library behind an interface
├── python/opengrid/       Python package
│   ├── _core.*.so          (built extension)
│   ├── model/              thin Python facade over core entities (Phase 1+)
│   ├── commands/           Command base class, CommandStack (undo/redo), Journal
│   ├── project/            .ogp read/write and schema migration
│   ├── api.py              public scripting API (= the command layer)
│   └── app/                PySide6 GUI (main window, viewport, tree, panels, console)
├── tests/python/          pytest (bindings, commands, GUI smoke with pytest-qt)
├── examples/              runnable phase demos; `examples/output/` holds their results
├── tools/                 repo tooling (e.g. check_no_ternary.py)
└── docs/                  SPEC, ARCHITECTURE, ROADMAP, PROGRESS
```

Rule: the dependency direction is strictly `app → api/commands → model → _core → og_core`.
`og_core` never includes Python, Qt or VTK headers. VTK is used in Python only for rendering
and picking, and in tests for read-back validation.

*Revised in Phase 3:* the VTKHDF and Plot3D writers are written in Python (h5py + NumPy, in
`opengrid.io`), not in C++. They are pure array I/O, so speed isn't a concern, and headless
scripted runs still need no VTK. CGNS (Phase 4) still needs the CGNS C library; where it goes
will be decided then.

## 2. Dependencies and licenses

| Dependency | Use | License | Linking | Status |
|---|---|---|---|---|
| C++ std lib (libc++ / libstdc++) | core | Apache-2.0 w/ LLVM exc. / GPL w/ runtime exc. | dynamic | required |
| Eigen 3 | small dense/sparse linear algebra | MPL-2.0 | header-only | required |
| nanobind | Python bindings | BSD-3 | static into extension | required |
| scikit-build-core, CMake, Ninja | build | Apache-2.0 / BSD-3 | build-time | required |
| GoogleTest | C++ tests | BSD-3 | test-only | required |
| HDF5 | VTKHDF + CGNS writing | BSD-style | dynamic | required (Phase 3) |
| CGNS (CGNS/MLL) | CGNS export | zlib | dynamic | required (Phase 4) |
| Shewchuk robust predicates | orient2d/3d, incircle/insphere | public domain | vendored source | required (Phase 5) |
| Python, NumPy | runtime | PSF / BSD-3 | — | required |
| PySide6 / Qt 6 | GUI | LGPL-3.0 | dynamic (relinkable) | required for GUI |
| VTK, pyvista, pyvistaqt | rendering, picking | BSD-3 / MIT / MIT | dynamic | required for GUI |
| meshio, h5py | test read-back, import | MIT / BSD-3 | — | test/import |
| pytest, pytest-qt, ruff, clang-format/tidy | dev tools | MIT / MIT / MIT / Apache-2.0 | — | dev only |
| **Netgen** | 3D volume tet generation | LGPL-2.1 | dynamic, optional plugin | Phase 8 |
| **mmg** (mmg3d/mmgs) | 3D tet optimisation / adaptation | LGPL-3.0 | dynamic, optional plugin | Phase 8 |
| OpenCASCADE | STEP/IGES import | LGPL-2.1 + exception | dynamic, optional plugin | Phase 9 |

**Not used by default (license incompatible with closed in-house use):** TetGen (AGPL-3),
Shewchuk's *Triangle* program (non-free for commercial use — note: only his *predicates* are
public domain and are used), GPL parts of CGAL, and the **Gmsh library** (GPL-2+). Gmsh's
`.msh` 4.1 format is written by our own code; Gmsh/meshio are used only to *read back* files in
tests. Any of the excluded tools may later be added only as an isolated, clearly flagged
plugin that is off by default and never linked into `og_core`.

LGPL compliance: all LGPL libraries are linked dynamically and isolated in `backends/`, so an
end user can replace them. This is noted again for Phase 9 packaging (app bundle/AppImage
must keep Qt and backend libraries as separate replaceable shared libraries).

## 3. Build system

- CMake ≥ 3.28, Ninja, C++20, double precision, `std::int64_t` entity/mesh indices.
- `CMakePresets.json` provides `debug` and `release` presets (build dirs `build/debug`,
  `build/release`) used for C++ development and `ctest`.
- Python packaging via **scikit-build-core** (`pyproject.toml`); the same top-level
  `CMakeLists.txt` builds the extension. Developer install:
  `pip install --no-build-isolation -e .` (editable, rebuilds on import is disabled by default
  for predictability; rerun the install after C++ changes).
- Compilers: Apple Clang ≥ 15 (tested with 16 from the Command Line Tools) and GCC ≥ 11.
  On this development machine `xcode-select` points to an old Xcode (Apple Clang 13, whose
  C++20 library is incomplete). Rather than hard-coding a machine-specific compiler path in
  the presets, the `opengrid` conda env sets `DEVELOPER_DIR=/Library/Developer/CommandLineTools`
  (Apple Clang 16), and `CMakeLists.txt` refuses Apple Clang < 15 with an explanatory message.
- Compile flags: `-Wall -Wextra -Wpedantic -Wshadow -Wconversion` (warnings as errors in CI),
  and **`-ffp-contract=off`** everywhere to keep floating-point results identical between
  debug/release and independent of FMA contraction decisions (determinism requirement).
- Options: `OG_BUILD_TESTS` (ON), `OG_BUILD_PYTHON` (ON), `OG_WITH_NETGEN`, `OG_WITH_MMG`,
  `OG_WITH_OCCT` (OFF).

## 4. C++ / Python boundary

- **Ownership (revised in Phase 1):** the authoritative entity store lives in **Python**
  (`opengrid.model.Store`). It holds small immutable parameter records, with stable integer IDs
  that are assigned monotonically and never reused, plus unique user-visible names. It also
  caches the C++ objects built from them (`_core.Line`, `_core.CubicSpline`, ...). The algorithms
  stay in C++ and take explicit inputs (curves, distributions, and later meshes) instead of a
  store reference.
  *Why the change:* the store is not a hot path, keeping it in Python avoids a C++ JSON
  dependency and a large binding surface, and it makes `.ogp` determinism trivial (Python's
  `repr` floats round-trip exactly). Heavy data (point arrays, meshes) still lives in C++/NumPy.
  If later phases need topology queries on the C++ side, a C++ topology view will be built from
  the store per operation.
- **Data exchange:** point/cell arrays cross the boundary as NumPy arrays via
  `nb::ndarray` (zero-copy views of core buffers where lifetimes allow, copies otherwise).
- **Long operations** (smoothing, triangulation, inflation, export) are C++ functions that
  take a `ProgressToken` (atomics: fraction done, residual history ring buffer, cancel flag).
  Bindings release the GIL (`nb::call_guard<nb::gil_scoped_release>`). The GUI runs them in a
  `QThread` and polls the token with a `QTimer` — C++ never calls back into Python from a
  worker thread, which avoids GIL re-entry and keeps cancellation simple.
- **Errors:** validators return structured reports (entity name, index, message) rather than
  throwing; true programming/IO errors throw C++ exceptions translated to Python exceptions.

## 5. Command / undo / journal design

- Every model mutation is a `Command` (Python, `opengrid.commands`). The GUI and the scripting
  API both call commands; there is no other mutation path (enforced by convention + tests).
- `Command.apply(store)` performs the change. The default `SnapshotCommand` records a JSON
  snapshot of every entity it touches **before** and **after** (created/modified/deleted sets).
  - **undo** restores the *before* snapshot (delete created entities, restore modified ones,
    recreate deleted ones with their original IDs);
  - **redo** restores the *after* snapshot — it does not re-run the algorithm, so redo is fast
    and exact.
  Generated meshes are **not** part of snapshots: they are derived data, cached and
  invalidated by entity revision counters, and regenerated deterministically on demand.
- `CommandStack` holds the undo/redo lists; commands may be grouped into a macro (one undo
  step for e.g. "create domain from selection + auto-assemble edges").
- **Journal:** each command renders itself as one line of Python against the public API
  (`og.create_line(p0=(0, 0, 0), p1=(1, 0, 0), name="line-1")`). Undo/redo are journaled as
  `og.undo()` / `og.redo()`. Replaying the journal in a fresh session must reproduce the
  project exactly (a test asserts JSON equality of the resulting `.ogp`). Entities are
  referenced by **name** in the journal for readability.

## 6. Native project format `.ogp` v1

UTF-8 JSON, describing intent, not meshes. Keys are sorted and floats written with
`repr`-exact round-tripping so that save → load → save is byte-identical.

```json
{
  "format": "opengrid-project",
  "version": 1,
  "units": "m",
  "tolerance": 1e-9,
  "next_id": 42,
  "name_counters": {"line": 3, "spline": 1},
  "database": [
    {"id": 1, "name": "le", "type": "point", "xyz": [0.0, 0.0, 0.0]},
    {"id": 2, "name": "chord", "type": "line", "p0": [0.0, 0.0, 0.0], "p1": [1.0, 0.0, 0.0]},
    {"id": 3, "name": "arc-1", "type": "arc", "center": [0, 0, 0], "radius": 10.0,
     "normal": [0, 0, 1], "start_angle": 0.0, "end_angle": 3.14159},
    {"id": 4, "name": "upper", "type": "spline", "method": "natural|clamped|catmull_rom",
     "points": [[0, 0, 0], [0.5, 0.06, 0], [1, 0, 0]],
     "start_tangent": [0, 1, 0], "end_tangent": null},
    {"id": 5, "name": "n1", "type": "nurbs", "degree": 3, "knots": [], "control_points": [],
     "weights": []}
  ],
  "connectors": [
    {"id": 10, "name": "con-1",
     "segments": [{"curve": 4, "t0": 0.0, "t1": 1.0}],
     "dimension": 101,
     "distribution": {"type": "tanh", "mode": "arc_length",
                      "begin_spacing": 1e-5, "end_spacing": null,
                      "breakpoints": []},
     "propagate_spacing": true}
  ],
  "domains": [
    {"id": 20, "name": "dom-1", "type": "structured",
     "edges": [[{"connector": 10, "reversed": false}], [], [], []],
     "generator": {"method": "tfi_linear", "elliptic": null}},
    {"id": 21, "name": "dom-2", "type": "unstructured",
     "loops": [[{"connector": 11, "reversed": false}]], "elements": "tri",
     "size_field": {}}
  ],
  "blocks": [],
  "patches": [{"name": "wall", "type": "wall",
               "members": [{"kind": "connector", "id": 10}]}],
  "zones": [{"name": "fluid", "members": [20, 21]}],
  "settings": {}
}
```

As implemented in Phase 1: one entity record per line, keys sorted, and floats written with
`repr` (exact). Explicit names of the form `<kind>-<n>` advance the automatic-name counter, so
replaying a journal (which always passes explicit names) reproduces `name_counters` too.

Fields are added per phase; every format change bumps `version` and adds a migration function
`migrate_v{n}_to_v{n+1}` in `opengrid.project`, with a test loading a stored v{n} file.

### Curve parametrisation (Phase 1)
All curves are parametrised on t ∈ [0, 1]. Polylines and splines use normalised cumulative chord
length, arcs use the angle and NURBS use their knot domain mapped linearly. Each curve
reports *breakpoints* (knots) that bound its smooth segments. Arc length (adaptive
Gauss–Kronrod 7/15 with a round-off floor), arc-length inversion (safeguarded Newton,
1e-14·L) and projection (sampling, then golden-section search, then Newton on (x−p)·x′) all work segment by
segment, so derivative jumps at knots never cost accuracy.

### Connectors and distributions (Phase 2)
- **Distributions** (`core/…/dist`) work on the unit interval:
  - uniform;
  - geometric, given exactly one of begin spacing, end spacing or growth ratio;
  - tanh and sine, each one- or two-sided.
- **tanh** is Vinokur's (1983) tanh/tan family.
- **sine** is the companion sinh/sin family. The spec only says "sine", which I've read as the
  Hyperbolic-sine distribution: one-sided uses sinh to cluster and sin to
  de-cluster, two-sided uses sin to cluster and sinh to de-cluster. The two-sided sine family
  cannot cluster an end below about (π/2n)²; such requests fall back to tanh with a warning.
- **Discrete spacing matching.** Vinokur prescribes end *derivatives*. An outer Brent
  (one-sided) or damped Newton (two-sided) iteration makes the *discrete* first and last
  intervals match the request, to about 1e-13 relative.
- **Far-end precision.** All formulas are evaluated without cancellation, and the far end
  comes from the complement 1 − s. Even so, a spacing δ next to the far end is the difference
  of two positions near 1, which limits its relative accuracy to about ε/δ. Spacings of
  1e-7·L or more meet the 1e-8 target; smaller far-end spacings are still finite and monotone,
  just less exact.
- **Fallbacks.** Infeasible requests always yield a finite, strictly increasing distribution
  plus a warning, never an exception or a NaN. Examples: spacing ≥ length, end spacings that
  sum to ≥ 1, too few points, or ratios outside floating-point range.
- **Connectors** reference curves by entity id, so renames are safe. They are chains of
  (curve, t0, t1) pieces (`og::CurveChain`), with optional breakpoints that create
  sub-connectors.
  - Spacings are absolute lengths.
  - In arc-length mode the distribution maps to arc length. In parametric mode it maps to the
    chain parameter, so the achieved physical spacing differs on curves with non-uniform
    speed.
- **Dependencies.** A curve used by a connector cannot be deleted, and it cannot be modified
  in a way that breaks a connector's chain. Both are rejected with a validation error naming
  the connector.
- **Spacing propagation** is pull-based and one level deep. A tanh/sine end with no spacing,
  on a connector with `propagate_spacing`, adopts the smallest explicit spacing that other
  connectors prescribe at the same node. It is evaluated on demand, so undo/redo and edits
  stay consistent automatically.
- **Derived data.** Nodes, chains and diagnostics are cached by the Session's
  `ConnectorEvaluator` and invalidated on every model change.
- **Stretching threshold.** The warning threshold for the stretching ratio (default 1.3) is a
  session/GUI preference. It is not stored in the project.
- **Quadrature changed in Phase 2** from recursive to globally adaptive (QAG-style, bounded
  number of subintervals). The recursive version had exponential worst-case cost on integrands
  carrying round-off noise. That case was real: the Phase 1 spline derivative cancelled on
  short segments, and building a finely spaced spline's length table took 7.5 s. Splines are
  now evaluated in difference form.

### Structured 2D domains (Phase 3)
- **Edges.** A domain is four connector chains listed as a counter-clockwise loop (edge 0 =
  j=0 with i increasing, then i=imax, then j=jmax with i decreasing, then i=0).
- **Validation.** The validator checks chain continuity, loop closure, equal point counts on
  opposite edges (naming the connectors and counts), a CCW orientation, and that the domain is
  planar (z = const).
- **Automatic assembly** handles two cases. The first is a single loop: corners are chosen
  among the junctions so that opposite counts match, preferring the sharpest turns. The second
  is the O-type: two closed loops joined by one cut connector, with the inner loop placed at
  j = 0. C-grids, which use a connector on both sides of the cut, are built from explicit edges.
- **Shared nodes are bitwise identical.** Coincident connector ends are snapped to one
  canonical coordinate (the lowest-id connector's), so domains sharing connectors or cuts have
  identical points. VTKHDF export then merges points exactly (np.unique), which is
  deterministic.
- **TFI.** Linear TFI uses Soni's arc-length blending parameters, which makes it exact on an
  annular sector. Hermite TFI uses wall-normal derivatives scaled by |top − bottom|; near the
  end columns the direction blends into the side boundaries'. *Limitation:* on C-grids
  (|T − B| ≫ wall size, cusp at the TE) Hermite TFI can fold. Use plain TFI plus elliptic wall
  control there.
- **Elliptic smoothing** (C++, see the header of `structured/elliptic.hpp` for references):
  - Controls are written in the α/γ-scaled (Thomas–Middlecoff) form, so P and Q stay O(1) even
    at a 1e-5 wall spacing.
  - The Thomas–Middlecoff control preserves boundary clustering.
  - TTM attractions are applied in the same scaled form.
  - **Wall control** uses Steger–Sorenson exponential-decay sources, with p and q set by
    *direct feedback* on the achieved spacing and angle (Hilgenstock). The classic GRAPE
    update left about 6% spacing and 4° angle errors. Defaults: decay 0.1/line, gain 0.1,
    clamp 100. When a source stays at the clamp, the run is reported as a warning, because the
    iteration can still "converge" without meeting the targets. That case was found on the
    C-grid with an earlier clamp of 4.
  - **Solvers.** `line_sor` is alternating-direction line relaxation (η-lines then ξ-lines),
    with source terms implicit along the lines. With explicit sources, line SOR diverged on the
    C-grid; alternating directions sped up convergence by about 500× there. `point_sor` is
    also available.
  - **Convergence** is relative (tolerance × first residual), with an absolute floor at
    1e-13·L for grids that already solve the equations.
  - **Threading.** The solver releases the GIL and reports to a `SmoothingProgress` (atomics
    plus a mutex-protected residual history, and a cancel flag). The GUI runs it on a Python
    thread and stores the result on the GUI thread.
- **Derived data.** Smoothed grids are cached by a key over the generator settings and the
  boundary node bytes. Unrelated edits keep a smoothing result valid, and undo brings a cached
  result back.
- **VTKHDF.** The layout was checked against files written by the installed VTK 9.7
  `vtkHDFWriter` (format 2.8). We declare version 1.0, because only 1.0 features are used.
  String field data must be variable-length *ASCII*; VTK fails to read UTF-8 there. **meshio
  5.3.5 has no VTKHDF (or Plot3D) reader**, so read-back is validated with VTK's
  `vtkHDFReader` plus h5py layout checks rather than with meshio as the spec asks.

## 7. 3D tet backend decision

`TetBackend` is a C++ interface: input = closed, watertight, consistently oriented surface
mesh (tris + quads) with region seeds and a size field; output = tets (+ pyramids for quad
faces), with boundary faces preserved bit-for-bit.

- **First backend: Netgen (LGPL-2.1)** — it generates volume meshes from a closed surface mesh
  with boundary preservation, has a C++ API (`nglib`), and is permissively usable when linked
  dynamically.
- **Second backend: mmg (LGPL-3.0)** for *optimisation and adaptation* (mmg3d with
  boundary frozen). mmg cannot create a volume mesh from a bare surface on its own, so it cannot
  be the only backend; it complements Netgen.
- **Later: own 3D Delaunay** (Bowyer–Watson + boundary recovery) behind the same interface,
  removing the external dependency once mature.
- TetGen (AGPL) is excluded; it may only ever appear as an off-by-default isolated plugin.
- Availability of Netgen/nglib and mmg as conda-forge packages for osx-arm64 and
  linux-64 will be verified at the start of Phase 8; fallback is building them from source as
  CMake ExternalProjects.

## 8. GUI architecture

- `MainWindow` (QMainWindow): central `pyvistaqt.QtInteractor` viewport; docks for the model
  tree (Database / Connectors / Domains / Blocks / Patches / Zones), property panel, log panel
  and an embedded Python console (namespace pre-populated with `opengrid` and, from Phase 1,
  the session's `og` API object).
- The GUI observes the model through a `ModelChanged` signal emitted by the CommandStack
  after each command/undo/redo; views redraw from the store. Views never mutate the model.
- High-DPI: Qt6 is HiDPI-aware by default; VTK render window DPI is set from the screen's
  device pixel ratio.

## 9. Coding conventions (summary — full list in CLAUDE.md)

- clang-format / clang-tidy for C++, ruff for Python.
- **No inline ternary conditionals** (`?:` in C++, `a if c else b` in Python) — explicit
  if/else. Enforced by `tools/check_no_ternary.py` in CI (neither clang-tidy nor ruff has a
  rule for this).
- One generator per file with a header comment giving method and reference.
- Validators never silently fix input.

## 10. Deviations from SPEC and identified risks

1. **Determinism.** "Bit-identical" is achievable for *same platform + same compiler + same
   build type*. Across arm64/x86_64 or Apple Clang/GCC, `libm` transcendental results differ in
   the last ulp, so we guarantee bit-identity per platform/build and tolerance-level agreement
   (1e-12 relative) across platforms. Measures: `-ffp-contract=off`, no `-ffast-math`, fixed
   reduction orders, no unordered-container iteration in output paths, fixed RNG seeds, and
   OpenMP/TBB reductions (Phase 9) done in deterministic chunk order.
2. **Ternary ban** needs custom tooling (see §9); the C++ check is a lexical heuristic, with
   `clang-query` used for an exact AST check when available.
3. **Toolchain on macOS.** Apple Clang < 15 is not supported (C++20 library gaps). We also avoid
   `std::format` and floating-point `std::to_chars` in core to keep the macOS deployment
   target low; formatting uses streams.
4. **VTKHDF patches.** The partitioned-collection layout in VTKHDF depends on the VTK version
   actually installed; it will be checked against that VTK's reader in Phase 3 before
   implementing, as the spec requires. The fallback (boundary grid + `patch_id` + name table)
   is always written as well, so meshio users get the patch information.
5. **Performance target "3D 5M tets via backend < 2 min"** is at risk with Netgen, whose
   optimisation passes are largely serial; we will measure it in Phase 8 and may meet it only
   with optimisation reduced or delegated to mmg.
6. **Phase size.** Phase 5 (CDT + refinement + advancing front + quads + smoothing) and
   Phase 8 are large. Proposal: split Phase 5 into 5a CDT/refinement/size field, 5b advancing
   front, 5c quads + smoothing; split Phase 8 into 8a surface meshing + shell + backend,
   8b 3D inflation, 8c optimisation + mixed blocks. Not applied to ROADMAP until approved.
7. **Ruppert refinement with small input angles** does not terminate in general; we will use
   concentric-shell splitting plus a "no refinement of skinny triangles at small input angles"
   rule (Shewchuk's Terminator) and document the achievable angle bound.
8. **Plot3D 2D** has several variants (with/without iblank, z-coordinate present or not); we
   will write the most common `multi-grid, 2D, no iblank` form and document it, plus an
   option for the 3D-with-one-k-plane form many solvers expect.
9. **Linux CI GUI tests** need an X server: we run them under `xvfb-run` with VTK's
   offscreen-capable build; on failure they are skipped with a clear marker, not silently.
10. **Trackpad interaction.** VTK interactor styles handle macOS pinch/two-finger pan poorly;
    a custom interactor style mapping Qt gesture events will likely be needed in Phase 1.
