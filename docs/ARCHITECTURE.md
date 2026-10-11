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
  "patches": [{"id": 30, "name": "wall", "type": "patch", "patch_type": "wall",
               "members": [{"connector": 10}]}],
  "zones": [{"id": 31, "name": "fluid", "type": "zone", "members": [{"domain": 20}]}],
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

### Patches, zones and solver-ready export (Phase 4)
- **Patch** entity: `patch_type` (wall, inlet, outlet, symmetry, farfield, periodic,
  interface, user) and member connectors. **Zone** entity: member domains. Both are stored in
  the `.ogp` `patches`/`zones` sections and journaled by name. *Bug found and fixed:* a patch
  parameter was first called `type` and overwrote the record's entity `type` on save. The
  writer now refuses any parameter that collides with `id`/`name`/`type`.
- **Topology.** Every use of a connector in a domain edge is recorded with its index range.
  - One use: a boundary connector, which must be in exactly one patch.
  - Two uses: an interface, detected automatically, including cuts within one domain (O- and
    C-grids). Interfaces must not be in a patch.
  - More than two uses: non-manifold (an error).
  - Each domain may be in at most one zone. Unzoned domains form an implicit zone `default`
    in every export (the name is then reserved).
- **1-to-1 connectivity** comes from the matched connector end nodes. PointRange is sorted,
  PointRangeDonor holds the images of its ends, and the transform maps the running axes
  (signed) and the face-normal axes (+ for min↔max faces, − for faces of the same kind).
  Examples: O-cut i_min↔i_max gives [1, 2]; the C-grid wake cut on j_min↔j_min gives [−1, −2].
- **Periodic pairing** (added after the Phase 4 review): stored once, on the source patch,
  as `periodic = {partner, rotation_center, rotation_angle, translation}` with the CGNS
  convention x_target = R(angle)(x − c) + c + t. The validator requires every periodic patch
  to be paired exactly once and every source connector to map node by node (either direction)
  onto a distinct target connector. Exports: CGNS periodic GridConnectivity1to1 in both
  directions (inverse transform on the reverse side; translations only in 2D), Gmsh
  `$Periodic` (affine matrix + node correspondences), VTKHDF field data on the patch blocks.
  Paired periodic faces are connectivity, never BCs. Projects load pairings in a second pass
  because they reference other patches of the same section.
- **Exports refuse to run** while the boundary does not validate (CGNS, Gmsh, and VTKHDF when
  patches or zones exist). The issue list names each connector, domain and patch involved.
- **CGNS** (C++ on the CGNS mid-level library 4.5, optional via `OG_WITH_CGNS`; zlib-like
  license):
  - One base (CellDimension 2, PhysicalDimension 2) with DataClass Dimensional and SI units.
  - Coordinates carry their dimensional exponents.
  - One structured zone per domain, with a zone FamilyName for its OpenGrid zone.
  - One Family_t per patch with FamilyBC (wall → BCWall, inlet → BCInflow, outlet →
    BCOutflow, symmetry → BCSymmetryPlane, farfield → BCFarfield, others →
    BCTypeUserDefined). One BC_t (FamilySpecified + FamilyName) per boundary connector.
  - GridConnectivity1to1_t in both directions for every interface.
  - Names are limited to 32 characters; longer ones are reported, not truncated.
  - `cgnscheck` reports no warnings on the demo files.
  - Build note: the CGNS CMake package refers to HDF5's imported target, so HDF5's config is
    loaded first; the conda `zlib` package is needed for the libz link.
- **Gmsh MSH 4.1** (Python): one surface entity and 2D physical group per zone, one curve
  entity and 1D physical group per patch. All nodes are merged and listed in one block. Read
  back with meshio, the physical groups and face counts are verified.
- **VTKHDF with patches**: a PartitionedDataSetCollection (format 2.1).
  - One quad block per zone (`zone_id`, `block_id`, cell indices, quality) and one line block
    per patch (`patch_id`), grouped in the Assembly as Volume/<zone> and Patches/<patch>
    (soft links).
  - *Found while testing:* VTK's reader iterates links in creation order and hangs on files
    whose groups do not track it. All VTKHDF groups are now created with `track_order=True`,
    as VTK's own writer does.
  - Without any patch or zone, the Phase 3 single UnstructuredGrid is still written.

### Unstructured 2D domains (Phase 5a)
- **Entity.** `unstructured_domain`: an outer loop (CCW) plus holes (CW), each a closed chain of
  connectors, together with generator settings (min angle, growth rate, max size, sources, max
  points). Loops can be assembled automatically from a set of connectors (largest loop =
  outer). Boundary nodes are exactly the connector nodes, so structured and unstructured
  neighbours conform and merge bitwise in the exports. Zones and patches accept both domain
  kinds. Topology counts loop uses, so interfaces between the kinds are detected.
- **Mesher** (C++, `unstructured/triangulate.cpp`, own implementation, no third-party code):
  - Predicates (`util/predicates.cpp`): orient2d/incircle with Shewchuk's stage-A filters and
    an exact expansion-arithmetic fallback, written from the 1997 paper.
  - Bowyer–Watson insertion inside a super triangle; Sloan's flip-based segment recovery;
    Lawson restoration of the constrained Delaunay property; flood-fill classification from
    the oriented segments.
  - Refinement: circumcenters for oversized triangles and Üngör off-centers for bad shapes.
  - *Boundary segments are never split*: points across a segment or in its diametral lens
    (120°) are rejected. The centroid is tried only for oversized interior triangles. A
    triangle's target size is at least its longest boundary edge.
  - The size field is the minimum of boundary spacing + (g − 1)·distance, the maximum size,
    and point/line sources, evaluated through bucketed ring search (exact).
- **Found while testing:**
  - A centroid fallback for badly shaped triangles cascades into ever thinner triangles.
  - Triangles standing on a long boundary segment next to a much finer boundary can never meet
    the size field. Both made refinement loop until the point cap; fixed as described above.
- **Limitation.** Where the boundary spacing jumps abruptly (tenfold at a corner), triangles
  next to the coarse segments cannot be well shaped without splitting the boundary.
  Quality ≥ 20.7° holds for consistently spaced boundaries and away from such jumps.
- **Performance:** about 1.4 M triangles/s (1.39 M triangles in 1.0 s on the M1).
- **Exports.** The mesh assembly now holds mixed cells (quads and triangles; offsets +
  connectivity + VTK types). Boundary lines and periodic correspondences come straight from
  connector nodes. Gmsh writes triangles (type 2); VTKHDF writes mixed cells. CGNS unstructured
  zones are deferred to Phase 5c (export refuses with a message), and Plot3D refuses
  unstructured domains.

### Advancing front (Phase 5b)
- **Where.** `unstructured/advancing_front.cpp` (own implementation; Lo 1985, Peraire et al.
  1987), selected per domain with `generator.algorithm = "advancing_front"`. The default stays
  `"delaunay"`, and projects saved before the option existed load as Delaunay.
- **Front.** Shortest edge first. The ideal apex sits at equilateral height for the local
  target size, which is clamped to [0.55, 1.3]× the edge length.
- **Candidates.** Existing front nodes within 0.9 h of the ideal apex come first, then the
  ideal point, then other nearby nodes; the first valid candidate is taken. All validity tests
  use the exact predicates.
- **Reuse radius.** It must exceed 1/1.3 ≈ 0.77. Neighbouring apexes lie about one edge length
  from the ideal point, and the front only coarsens if it reuses them. Measured on a graded
  cylinder: radius 0.8 gives median edge/target 0.48; 0.9 gives 1.01. A best-quality choice
  among candidates (instead of first valid) did not coarsen either.
- **Leftover edges.** Edges still open after retries go to the Delaunay mesher, with a warning.
  None occurred in the tests or demos.
- **Post-processing:**
  - Lawson flips on a triangle-neighbour array (segments are never flipped).
  - Three passes of guarded Laplacian smoothing: a move is undone if it lowers the smallest
    angle around the node; sines are compared instead of angles.
  - Flips again.
- **Bug found during development.** A triangle closing against an existing front edge was
  first rejected as overlapping that very edge, which stalled the front everywhere; fixed by
  exempting the edges being closed.
- **Shared speed-ups (both meshers):**
  - The size field moved from a bucket ring search (O(r³) cell visits for points deep inside a
    boundary loop) to a k-d tree branch and bound. The result is still the exact minimum, now
    1–8 µs per query.
  - `orient2d` has an allocation-free exact path for exact coordinate differences (Shewchuk
    stage B); it handles the 16% of tests the filter cannot decide in near-equilateral meshes.
  - Delaunay on a 6000-node disc: 2.8 s → 0.3 s, same mesh.

### Quads, smoothing, CGNS unstructured zones (Phase 5c)
- **Pipeline** (`unstructured/domain_mesher.cpp`): triangulation (Delaunay or advancing front),
  then by `generator.element_type`:
  - `triangles`: kept as is;
  - `quad_dominant`: recombined (`recombine.cpp`);
  - `quads`: coarse mesh + recombination + subdivision.

  Then the requested smoothing (with edge swapping), then untangling. With `smoothing = "none"`
  and triangles, the output is bit-identical to Phase 5a/5b (tested).
- **Recombination.** Greedy matching of triangle pairs by Remacle's quad quality η (only strictly
  convex quads with η ≥ 0.3), then bounded augmenting paths (≤ 4 matched pairs) to cut the
  leftover triangles. Typically 2–8% of cells stay triangles.
  - Blossom-Quad's minimum-cost perfect matching was not used: it needs Edmonds' weighted
    blossom algorithm, a large piece of code, and a perfect matching does not always exist on a
    given triangulation.
- **All-quad** ("recombine, then subdivide", as Gmsh's full-quad option):
  1. Mesh a coarse boundary made of every second node of each loop, at twice the target size.
  2. Recombine.
  3. Split every triangle into 3 quads and every quad into 4. The midpoint of a coarse boundary
     edge is the skipped connector node, so the boundary is exactly the connector nodes.
  4. Smooth the new nodes with three guarded Laplacian sweeps.

  Requirements and choices:
  - Each loop must have an even number of segments. This is stricter than the necessary
    condition (an even total); two odd loops are not handled.
  - Loops touching at a point are rejected for `quads`.
  - The quad-dominant leftovers become 3 quads each, so no triangles remain.
- **Smoothing** (`smoothing.cpp`, Gauss–Seidel over free nodes):
  - Methods: Laplacian, length-weighted, and angle-based (Zhou–Shimada 2000).
  - Every move is guarded on the worst signed scaled Jacobian around the node, so smoothing
    never inverts a cell or worsens the worst cell.
  - Edge swapping (Delaunay flips between triangle pairs; boundary segments fixed) after each
    sweep.
  - Untangling: compass search on the worst signed scaled Jacobian, started from the
    neighbours' centre. A cell that stays bad produces a warning, never a silent result.
- **Quality.** `quality_cells` computes all 2D metrics for mixed meshes (triangles first, then
  quads; the same definitions as structured quads). `quality_triangles` delegates to it.
- **CGNS** (`io/cgns_structured.cpp` → `write_cgns_2d`; the structured function is a wrapper):
  - Unstructured zones get `TRI_3`/`QUAD_4` sections numbered 1.., then one `BAR_2` section
    per BC. The BC_t is a PointRange of those bars with GridLocation EdgeCenter (SIDS face BCs
    for CellDim 2).
  - Interfaces and periodic pairs involving an unstructured zone use GridConnectivity_t
    (Abutting1to1, Vertex, PointList/PointListDonor, (i, j) pairs on the structured side;
    periodic translation as for 1-to-1).
  - Vertex numbers come from an exact coordinate lookup (boundary points are the connector
    nodes bitwise).
  - cgnscheck reports no warnings or errors on the tests and demos.
- **Measured at 1 M cells** (disc, uniform size):

  | Triangulation | Element type | Cells | Time | Min SJ |
  |---|---|---|---|---|
  | Delaunay | triangles | 1.41 M | 1.2 s | 0.58 |
  | Delaunay | quad-dominant | 0.74 M (8% triangles) | 2.9 s | 0.45 |
  | Delaunay | all-quad | 0.72 M | 2.3 s | 0.46 |
  | Advancing front | all-quad | 0.57 M | 2.0 s | 0.62 (mean 0.89) |

### 2D inflation layers (Phase 6)
- **Placement.** Inflation is part of the unstructured domain generator
  (`generator.inflation`: wall connectors from the domain's loops, `first_height`,
  `growth_rate`, `layers`, `total_thickness`, `transition`).
  - `unstructured/inflation.cpp` marches the layers.
  - `domain_mesher.cpp` meshes the core, bounded by the final fronts and the non-wall
    boundary, with any Phase 5 generator. Core boundary nodes are the front nodes, so the mesh
    is conforming.
  - Boundary nodes are still the connector nodes. Patches, zones, periodicity and all exports
    work unchanged.
- **Marching:**
  - Layer heights are h₀·gᵏ.
  - Normals are front bisectors, Jacobi-smoothed with a visibility constraint, and scaled by
    1/min(n·n_seg). On a straight wall layer k is exactly h₀gᵏ (tested to 1e-6); on curved
    walls the first height is within 0.25% (demo; tested < 2%).
  - Sharp convex corners (normals turning > 60°) get fans in steps of ≤ 20°.
  - Converging columns are merged.
- **Stopping rules:**
  - The layer limit is reached.
  - **Transition**: the layer height reaches the column's isotropic size. That size is its
    width along the front, at least the wall spacing it started from, at most the core size
    field.
  - **Slope**: a column may not stand above a stopped neighbour by more than one layer and
    more than `max_slope` (0.5) times their distance along the front.
  - **Crossing**: exact intersection tests against the old fronts, the boundary and the other
    new segments.
  - **Proximity**: a facing segment ahead closer than the isotropic size; this handles
    multi-element gaps.
  - **Cell validity**: inverted cells, or an angle above 165°.

  Checks repeat until no column has to stop, so every committed layer is valid. Open wall
  chains keep their fixed end nodes, and the layers taper to them.
- **Lessons from the NACA demo** (each is fixed and covered by a test or the demo):
  - The slope rule must measure the gap along the front, not between the column tops; that
    distance already contains the height difference, so the rule never fired.
  - It must run to a fixpoint after all transition stops; a single pass was order-dependent
    and left cliffs.
  - It must allow single-layer steps (with a tolerance), or one stop sweeps the whole front
    once layers are thicker than the column spacing.
  - Transition on the live front spacing alone stops fan columns after one layer. Transition
    on the core size alone lets the layers grow to about a chord when the layer growth rate is
    close to the size-field growth rate.
- **All-quad cores.** The coarse boundary keeps sharp corners (turn > 35°, `loop_corners.hpp`)
  and pairs the nodes between them; inflation makes each stretch even by merging a pair of
  layer nodes on a smooth part (with fallbacks). This also changed Phase 5c all-quad meshes
  slightly (NACA demo: 10756 instead of 10760 quads, min SJ 0.665 instead of 0.683).
- **Optimisation.** After smoothing, nodes of cells with scaled Jacobian < 0.3 are optimised
  (worst-cell compass search, the untangler with a threshold). Well-shaped meshes (25°
  triangles score 0.49) are untouched.
- **Limitations:**
  - A layer front ending at an open wall end tapers to zero height there, because the
    neighbouring connector's nodes are fixed.
  - Layer cells are not smoothed.
  - `max_slope`, the fan angles and the proximity factor are C++ options, not yet exposed in
    Python.

### 3D structured foundation (Phase 7a)
- **Surfaces** (`geom/surface.cpp`): plane, ruled, extruded, revolved and NURBS
  (Piegl–Tiller basis), with analytic first derivatives.
  - Projection: sample grid, then Newton on the stationarity conditions with a
    finite-difference Hessian (Gauss–Newton where it is not convex), damped and clamped.
  - Accuracy: on-surface points are recovered within 1e-10·size (tested for all five types).
  - Model entities are database kinds (`plane`, `ruled_surface`, `extruded_surface`,
    `revolved_surface`, `nurbs_surface`). Curve references are entity ids. Cached C++ surfaces
    are dropped when a referenced curve changes.
- **Structured domains** take an optional `surface`:
  - Parametric TFI on it (`structured/surface_tfi.cpp`): project the boundary nodes, linear TFI
    of (u, v), evaluate. The domain must not cross a parameter seam.
  - Without a surface, domains outside a plane z = const use linear TFI in space.
  - Hermite TFI and elliptic smoothing stay restricted to z = const planes.
- **Blocks** (`model/blocks.py`, entity `block`):
  - From six domains in any order and orientation: the 8 corners are identified by exact
    coordinates, a right-handed frame is chosen at the first face's first corner, and each face
    is transposed/flipped onto its block face. Then 3D TFI (`structured/tfi3d.cpp`: Boolean sum
    of three linear projectors; parameters are face arc-length fractions blended in index
    space; an annular sector is reproduced exactly).
  - Or by extrusion of a domain (`structured/extrude.cpp`): translate, rotate (Rodrigues) or
    along node normals, with a connector-style distribution. The base is mirrored in i if
    needed for a right-handed block.
  - Grids are (nk, nj, ni, 3). The block evaluator caches by a hash of the parameters and the
    face grids.
- **3D topology:**
  - Block faces with identical nodes (8 orientations checked) are 1-to-1 interfaces, including
    a block-internal cut such as an extruded O-grid wake. Other faces are boundaries.
  - CGNS 1-to-1 data: begin/end donor indices by coordinate lookup. Transform of the normal
    axis: `−side·donor_side`; in-face axes from a unit step.
  - Patch members are connectors (2D), structured domains (faces of six-face blocks, the base
    of extrusions) or block faces `{"block", "face"}`.
  - `boundary_report()` switches to the 3D rules when the model has blocks: every boundary face
    in exactly one patch; connector members and periodic patches are reported.
- **Exports** with blocks:
  - VTKHDF: hexahedra plus boundary quads (outward) per patch, quality as cell data.
  - Gmsh: hexahedra and quads, physical volumes and surfaces.
  - CGNS 3D structured (`write_cgns_structured_3d`): zone families, BC point ranges, 1-to-1
    connectivity.
  - Plot3D 3D multi-block.
  - Projects gain the `blocks` section. Entities load one by one in id order, so references
    within a section resolve.
- **Quality** (`quality/quality3d.cpp`):
  - scaled Jacobian (8 corners), face equiangle skewness, edge-family aspect ratio, face
    angles, FV non-orthogonality, volume ratio;
  - exact trilinear volume (2×2×2 Gauss).

  In the GUI's quality views, area and area ratio stand for volume and volume ratio for blocks.
- **GUI:**
  - Blocks in the tree and viewport (boundary wireframe; cells when cut or shrunk);
    translucent surfaces.
  - "Block from 6 faces", "Extrude domain", block property panel with an extrusion editor.
  - View → cut plane (x = model centre, whole cells kept) and shrink.
  - Quality colouring and the worst-cell locator for blocks.
- **Limitations:**
  - Periodic pairing of 3D faces is not supported.
  - Degenerate (collapsed) block faces are not supported.
  - 2D quality of non-planar face domains uses the planar formula.
  - The cut plane is fixed at x = model centre.

### 3D smoothing and hyperbolic extrusion (Phase 7b)
- **3D elliptic smoothing** (`structured/elliptic3d.cpp`):
  - Poisson system with cofactor metric coefficients, point SOR, faces fixed.
  - Control: Laplace, or Thomas–Middlecoff with P_m from the four faces containing the m-lines,
    blended by transfinite interpolation in index space.
  - Wall control on any face: decaying sources with Hilgenstock feedback, using the same signs
    as the 2D smoother. The normal source comes from the first-cell height error, the two
    tangential sources from the lean of the first grid line.
  - Bump-channel demo: first-spacing error 15% → 0.0%, wall angle 19.5° → 0.00° (852
    iterations, 0.4 s).
  - Blocks get it through the optional `smoothing` setting (faces or extrusion blocks alike).
- **Surface smoothing** (`structured/surface_elliptic.cpp`):
  - The 2D Poisson system written with 3D tangent vectors, Laplace or Thomas–Middlecoff; every
    updated node is projected back onto the surface (Newton from its previous parameters).
  - Used by structured domains with a `surface` and an elliptic generator; walls and
    attractions are rejected.
- **Hyperbolic extrusion** (`structured/hyperbolic.cpp`, extrusion mode `hyperbolic`):
  - An explicit marching variant of the hyperbolic conditions: layer normals smoothed with a
    visibility guard, steps divided by the cosine to the raw normal (corners keep the height),
    and tangential dissipation.
  - Edge conditions: free, constant x/y/z, periodic.
  - This is not Chan–Steger's implicit block-tridiagonal solver; it is cheaper and adequate for
    smooth convex surfaces. Strongly concave regions may need more dissipation, and inverted
    layers are reported as warnings.
  - Cylinder demo: radii exact to 1e-9, first cell exactly 1e-3, worst SJ 0.9997.

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
