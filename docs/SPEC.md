# Claude Code prompt pack — interactive block-structured / unstructured mesher

## How to use this file

1. Create an empty repo, copy this file in as `docs/SPEC.md`.
2. Start Claude Code in the repo and paste **Prompt 0 (Kickoff)**.
3. After each phase is merged and you have checked the demo, paste the next **Phase N** prompt.
4. Keep `CLAUDE.md` (Claude Code creates it in Phase 0) short: conventions, build commands, current phase. The long-lived spec stays in `docs/SPEC.md`.

Tip: start each phase in plan mode (Shift+Tab twice), approve the plan, then let it run. Ask it to stop at the end of the phase rather than continuing on its own.

---

## MASTER SPEC (paste once with Prompt 0, keep as docs/SPEC.md)

### Goal
Build **"OpenGrid"** (placeholder name), a desktop mesh-generation application: interactive, bottom-up (geometry → connectors → domains → blocks), supporting structured multiblock, unstructured and hybrid grids in 2D and 3D, with boundary-layer inflation, smoothing, patch (boundary condition) definition, and export to standard formats. Target users are CFD engineers; it must produce meshes directly usable by finite-volume and finite-difference solvers.

### Platforms
- macOS on Apple Silicon (arm64, native, no Rosetta) and Linux x86_64 (Ubuntu 22.04+). Optional Linux aarch64.
- CI on both (GitHub Actions: `macos-14` arm64 runner + `ubuntu-latest`).

### Architecture (proposed; justify any deviation in docs/ARCHITECTURE.md before coding)
- **Core library in C++20** (`core/`): geometry kernel, topology/data model, distribution functions, structured generators, unstructured generators, inflation, smoothing, quality metrics, I/O. No GUI dependencies. Double precision everywhere, 64-bit indices for mesh entities.
- **Python bindings** via nanobind (or pybind11) built with CMake + scikit-build-core.
- **GUI in Python**: PySide6 (Qt6) + VTK via `pyvista`/`pyvistaqt` for rendering and picking. Both have native arm64 macOS wheels.
- **Scripting API** in Python exposing the same command layer the GUI uses. Every GUI action is a command; the session can be replayed as a Python journal script.
- **Command pattern** for all model mutations → undo/redo and journaling come for free.
- **Native project format**: versioned JSON (`.ogp`) describing geometry, connectors, distributions, domains, blocks, patches and generator settings — not the mesh itself (meshes are regenerated deterministically).
- Environment: `uv` or conda-forge; document both in README. Build must work with Apple Clang and GCC ≥ 11.

### Licensing (decide in Phase 0, record in docs/ARCHITECTURE.md)
The project may be used inside an aerospace company. Prefer permissive or LGPL dependencies. Explicitly evaluate and document:
- Own implementations (preferred for 2D: CDT, refinement, advancing front, structured generators).
- Netgen (LGPL), mmg (LGPL) for 3D tet generation/optimization behind a pluggable backend interface.
- Avoid by default: TetGen (AGPL), Shewchuk's Triangle (non-free for commercial use), CGAL GPL parts. If used at all, only as optional, isolated plugins clearly flagged.
- OpenCASCADE (LGPL) optional, for STEP/IGES import in a later phase.

### Data model (hierarchy)
- **Database entities** (geometry): points, line segments, circular arcs, polylines, cubic splines (natural / clamped / Catmull-Rom), NURBS curves; later surfaces (planar, ruled, revolved, extruded, NURBS, imported CAD). Each with evaluation `x(t)`, derivative, arc-length parametrization, closest-point projection.
- **Connector**: a curve (or chain of curves) with a number of points and a **distribution function**:
  - uniform; geometric (growth ratio); hyperbolic tangent / sine with prescribed begin and/or end spacing (Vinokur 1983 two-sided stretching); explicit spacing at breakpoints (sub-connectors); arc-length vs. parametric mode.
  - Spacing constraints can be propagated to neighbouring connectors (shared nodes).
  - Display: node markers, spacing ratio coloring, max stretching ratio warning.
- **Domain** (2D entity, or surface patch in 3D):
  - *Structured*: 4 edges, each a chain of connectors; opposite edges must have matching point counts (validate, report clearly).
  - *Unstructured*: one outer loop + any number of inner loops of connectors; element type triangles, quads or quad-dominant.
- **Block** (3D):
  - *Structured*: 6 faces (structured domains), consistent i/j/k orientation, point matching between faces.
  - *Unstructured*: closed shell of domains (any mix of structured/unstructured faces) → tets, with prisms/hexes in inflation layers and pyramids at quad-to-tet transitions.
- **Patch** (boundary condition group): a named set of connectors (2D) or domain faces (3D) with a type tag (wall, inlet, outlet, symmetry, farfield, periodic, interface, user-defined). Every boundary entity must belong to exactly one patch before export (validator).
- **Zone / volume tag**: named groups of domains (2D) or blocks (3D) for multi-region meshes.
- Topology queries: shared connectors between domains, non-manifold detection, orientation consistency, multiblock connectivity (1-to-1 interfaces) for structured export.

### Algorithms required
Structured:
- 2D/3D transfinite interpolation (linear Lagrange blending; Hermite blending option for wall orthogonality).
- Elliptic smoothing: Winslow (Laplace) and Poisson with Thompson–Thames–Mastin control functions; Steger–Sorenson (or Thomas–Middlecoff) source terms to enforce wall spacing and orthogonality. SOR / line-relaxation solver with residual history shown in the GUI.
- Hyperbolic grid generation (Steger–Chaussee) for O/C-grid extrusion from a curve (2D) and surface (3D), with dissipation and marching-step control.
- Extrusion: translate, rotate (axisymmetric), normal extrusion.
- Domain/block splitting and joining at given indices; O-grid / butterfly topology helpers (later phase).

Unstructured:
- 2D constrained Delaunay triangulation (Bowyer–Watson + edge recovery, robust predicates — Shewchuk's orient2d/incircle code is public domain and acceptable) with Ruppert/Chew refinement and a size field (from boundary spacing, background function, sources, max size, growth rate).
- 2D advancing-front triangulation as an alternative.
- Quads: triangle-merge (Blossom-Quad or greedy merge) for quad-dominant; full-quad via split. Document the choice.
- 3D: surface triangulation of domains on surfaces (parametric-space meshing with metric correction), volume tetrahedralization via pluggable backend (own Delaunay later; Netgen/mmg first), pyramids at quad faces.
- **Inflation / boundary layers**: advancing-layer extrusion from selected wall patches — first-cell height, growth rate, number of layers or total thickness, with: normal smoothing at convex/concave corners, multiple normals at sharp corners (fan/wedge in 2D), local layer termination on collision/proximity, quality-based stopping, and smooth transition to the isotropic core (last-layer size matching the core size field). 2D → quads; 3D → prisms (from tri walls) / hexes (from quad walls), then tets in the core and pyramids where needed.

Smoothing / optimization:
- Laplacian, angle-based and length-weighted smoothing for unstructured meshes with boundary constraints (nodes on curves/surfaces stay on geometry by projection).
- Edge swapping and optimization-based untangling (minimize inverse mean-ratio / maximize min scaled Jacobian).

Quality metrics (per cell, histogram + colormap in GUI, worst-cell locator):
- Scaled Jacobian / negative-Jacobian detection, skewness (equiangle), aspect ratio, max/min included angle, non-orthogonality (FV sense: face normal vs. cell-centre vector), area/volume ratio between neighbours, wall y-spacing.

### I/O
Export (each validated by reading back with VTK and meshio in tests):
- **VTKHDF** (primary): check the current VTKHDF file-format spec in the installed VTK version before implementing — do not rely on memory. UnstructuredGrid for all meshes (structured blocks exported as hexes/quads). Patches as a multiblock / PartitionedDataSetCollection if the installed VTK supports it, otherwise boundary faces as a separate UnstructuredGrid with `patch_id` cell data plus a patch-name table in field data. Also write `zone_id` and quality fields.
- **CGNS** (aero standard, via the CGNS library — zlib-like license): structured zones with 1-to-1 connectivity and BC_t for structured; unstructured zones with element sections and BCs for unstructured/hybrid.
- **Gmsh `.msh` v4.1** (physical groups = patches) and **Plot3D** (multiblock structured, whole/multi-grid, binary and ASCII).
- Optional later: OpenFOAM polyMesh, SU2.
Import: native project, Plot3D, Gmsh, simple geometry (DAT airfoil coordinates, CSV point lists); STEP/IGES via OpenCASCADE in a later phase.

### GUI requirements
- Main window: 3D viewport (VTK), model tree (database / connectors / domains / blocks / patches / zones), property panel for the selected entity, command console (Python), log/message panel, undo/redo.
- Picking and multi-selection in the viewport; snapping to points and curve ends.
- Creation tools: point, line, arc, spline through points, connector on curve, domain from selected connectors (auto-assemble edges), block from selected domains.
- Edit tools: dimension (number of points), distribution type and spacings, split/join connectors, switch domain type structured ↔ unstructured, generator settings, run smoother with live residual plot.
- Visualization: show/hide by entity type, mesh edges, cell quality colormap, cut planes for 3D blocks, shrink view, patch colouring.
- Long operations run in a worker thread with progress and cancel; GUI must stay responsive.
- Works with trackpad on macOS and mouse on Linux; HiDPI aware.

### Engineering rules
- Each phase ends with: passing tests, a runnable demo script under `examples/`, a screenshot or exported file under `examples/output/`, and an update to `docs/PROGRESS.md`. Then **stop and wait for review**.
- Unit tests: GoogleTest (C++) and pytest (Python). Regression cases with checked invariants (cell count, min Jacobian > 0, total area/volume vs analytic, patch face counts, watertightness).
- Golden validation cases, introduced as phases allow: unit square (uniform/stretched), quarter annulus (TFI exactness), NACA 0012 C-grid and O-grid, cylinder O-grid with wall spacing, backward-facing step multiblock, airfoil hybrid (inflation + triangles), 3D channel/duct, annular sector, sphere in box (prisms + tets), wing stub.
- Code style: clang-format + clang-tidy for C++, ruff for Python. **Do not use inline ternary conditionals (`? :` in C++, `a if c else b` in Python); write explicit if/else blocks.**
- Prefer clarity of algorithms over cleverness: each generator in its own file with a header comment giving the method and reference.
- Never silently fix invalid input: validators report what is wrong and where (entity name, index).
- Determinism: same project + same settings → bit-identical mesh (fixed seeds, stable ordering).
- Performance targets (Release, M1): 2D 1M cells < 10 s; structured elliptic smoothing on 500×200 block converges interactively; 3D 5M tets via backend < 2 min.

---

## PROMPT 0 — Kickoff (paste this first)

```
Read docs/SPEC.md completely. We will build this application progressively, one phase at a time; do not implement beyond the current phase.

Phase 0 tasks:
1. Write docs/ARCHITECTURE.md: module layout, dependency list with licenses, the build system, the C++/Python boundary, the command/undo/journal design, the native project-file schema (v1), and the decision on 3D tet backends (justify against the licensing section). Flag anything in SPEC.md you think is wrong or risky and propose alternatives.
2. Write docs/ROADMAP.md with the phases below, each with acceptance criteria you can test automatically.
3. Scaffold the repo: CMake + scikit-build-core, nanobind module with a trivial function, PySide6 + pyvistaqt main window with an empty viewport, model tree and Python console. pytest + GoogleTest wired. clang-format, clang-tidy, ruff configs. GitHub Actions for macos-14 (arm64) and ubuntu-latest.
4. Create CLAUDE.md with build/test commands, conventions (including: no inline ternary operators, explicit if/else), and "current phase".
5. Verify: builds and launches on this machine, tests pass. Report what you could not verify.

Start in plan mode, show me the plan, and stop after Phase 0 for review.
```

## PHASE 1 — 2D geometry kernel + project I/O

```
Phase 1 (see docs/ROADMAP.md). Implement in core/: points, lines, arcs, polylines, cubic splines, NURBS curves with evaluation, derivatives, arc-length reparametrization (accurate to 1e-10 relative) and closest-point projection. Command layer with undo/redo and Python journaling. Native project save/load (.ogp v1). GUI: create/edit/delete these entities, picking, snapping, model tree, property panel. Import of airfoil .dat and CSV point lists.
Tests: arc-length of circle and spline vs reference, projection accuracy, project round-trip, undo/redo of every command.
Demo: examples/01_geometry.py builds a NACA 0012 + farfield circle by script and saves a project. Stop for review.
```

## PHASE 2 — Connectors and distributions

```
Phase 2. Connectors on curves and curve chains; dimension; distribution functions: uniform, geometric, tanh and sine with one- and two-sided spacing (Vinokur), breakpoints with sub-connector spacings, arc-length vs parametric mode. Spacing at shared nodes propagates to neighbours (optional per connector). GUI: node display, spacing colouring, warnings for stretching ratio > user threshold.
Tests: achieved end spacings vs requested (rel. error < 1e-8), monotonicity, max stretching ratio, behaviour at extreme ratios (fall back with a warning, never NaN).
Demo: airfoil connectors with LE/TE clustering. Stop for review.
```

## PHASE 3 — Structured 2D domains + first export

```
Phase 3. Structured domains from 4 edges (each a connector chain), automatic edge assembly from a selection, orientation handling, validation of opposite point counts. Generators: linear TFI, Hermite TFI, elliptic smoothing (Winslow; Poisson with Thompson-Thames-Mastin control; Steger-Sorenson wall control for spacing/orthogonality) with SOR and line relaxation, residual plot in GUI, run in a worker thread with cancel. Quality metrics (2D set from SPEC) with histogram and colormap. Multi-domain point matching on shared connectors.
Export: VTKHDF (check the spec in the installed VTK version first) and Plot3D 2D multiblock.
Tests: TFI exactness on quarter annulus, positive Jacobians, elliptic convergence, VTKHDF read back with vtk and meshio, Plot3D round trip.
Demos: NACA 0012 C-grid and O-grid with prescribed wall spacing; backward-facing step multiblock. Stop for review.
```

## PHASE 4 — Patches, zones and solver-ready export

```
Phase 4. Patches (named BC groups with types) and zones; validator ensuring every boundary connector belongs to exactly one patch and interfaces are detected automatically. Structured 1-to-1 connectivity detection between domains. Export with boundary info: VTKHDF (patches as multiblock/partitioned collection if supported, otherwise boundary grid + patch_id + name table), CGNS (structured zones, 1-to-1 connectivity, BC_t), Gmsh msh 4.1 with physical groups. GUI patch editor and patch colouring.
Tests: patch face counts, CGNS file validated with cgnscheck if available, Gmsh file read by meshio with correct physical groups. Stop for review.
```

## PHASE 5 — Unstructured 2D

```
Phase 5. Unstructured domains from outer + inner loops. Own constrained Delaunay (Bowyer-Watson, robust predicates, edge recovery) with Ruppert/Chew refinement driven by a size field (boundary spacing interpolation, max size, growth rate, point/line sources). Alternative advancing-front generator. Quad-dominant and all-quad options via triangle merging (document method). Laplacian/angle-based smoothing with boundary nodes fixed, edge swapping. Mixed structured + unstructured domains sharing connectors must be conforming.
Tests: boundary conformity (every boundary segment present exactly once), min angle guarantees where theory allows, total area vs geometry, no inverted cells, determinism.
Demo: airfoil in a circular farfield, triangles and quad-dominant. Stop for review.
```

## PHASE 6 — 2D inflation layers (hybrid)

```
Phase 6. Advancing-layer inflation from selected wall patches: first height, growth rate, number of layers or total thickness, normal smoothing, multiple normals/fans at sharp convex corners, layer collapse at concave corners, local termination on proximity/collision and on quality, transition so the outer layer matches the core size field. Layers are quads; the remaining region is filled by the Phase 5 generator.
Tests: first-cell height and growth achieved, no overlaps (segment intersection checks), positive Jacobians, trailing-edge (sharp and blunt) handling.
Demos: NACA 0012 hybrid mesh; two-element airfoil (gap proximity). Stop for review.
```

## PHASE 7 — 3D structured

```
Phase 7. Surfaces (planar, ruled, revolved, extruded, NURBS) with evaluation and projection. Structured domains on surfaces (parametric TFI + projection, surface elliptic smoothing). Structured blocks from 6 faces, orientation and point matching, 3D TFI, 3D elliptic smoothing with wall control, hyperbolic extrusion from a surface, translate/rotate/normal extrusion of 2D meshes. Cut planes and shrink view in the GUI. Export: VTKHDF (hexes), CGNS 3D structured with connectivity and BCs, Plot3D 3D.
Tests: block volume vs analytic, positive Jacobians, 1-to-1 interface consistency.
Demos: annular sector duct; extruded airfoil wing section; O-grid around a cylinder in 3D. Stop for review.
```

## PHASE 8 — 3D unstructured and hybrid

```
Phase 8. Unstructured surface meshing of domains on surfaces (parametric-space meshing with metric correction), closed shell assembly and watertightness check. Volume generation behind a TetBackend interface; first backend as chosen in ARCHITECTURE.md (e.g. Netgen or mmg), keep room for an own Delaunay later. Prism/hex inflation from wall patches (3D version of Phase 6: normal smoothing, termination, ridge/corner handling). Pyramids at quad faces adjoining tets. Optimization: smoothing + untangling with quality-based objective. Mixed blocks (structured + unstructured) conforming at shared faces.
Tests: watertightness, element-type counts, no negative volumes, total volume vs analytic, inflation height/growth achieved.
Demos: sphere in a box with prism layers; wing stub in a box hybrid mesh. Export to VTKHDF, CGNS (mixed element sections + BCs), Gmsh. Stop for review.
```

## PHASE 9 — CAD import, scripting polish, performance

```
Phase 9. Optional STEP/IGES import via OpenCASCADE (isolated module), healing tolerances, automatic connector/domain creation from CAD edges/faces. Python scripting API documentation and examples reproducing all golden cases from scripts only. Profiling against the performance targets in SPEC.md; parallelize hot loops (OpenMP or TBB). Packaging: app bundle for macOS arm64 and AppImage or conda package for Linux. Stop for review.
```

---

## Notes for you (not for Claude Code)

- Phases 1–4 alone already give a useful structured 2D tool; you can stop there and use it before investing in 3D.
- If a phase gets too big in practice, ask Claude Code to split it ("split Phase 5 into 5a CDT, 5b refinement/size field, 5c quads") and update ROADMAP.md.
- Phase 8 is by far the hardest; expect several iterations on inflation robustness at ridges and concave corners.
