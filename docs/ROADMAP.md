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
  quality at sharp edges.
- Surface meshing on parametric surfaces, a watertight shell check (each edge shared by exactly 2 faces, consistent orientation).
- `TetBackend` with Netgen first and mmg for optimisation.
- Prism/hex inflation with ridge/corner handling, pyramids at quad faces.
- Element-type counts as expected, no negative volumes, volume vs analytic within 1e-6, inflation height/growth achieved.
- Demos: sphere in a box with prism layers; wing stub hybrid. Export VTKHDF, CGNS (mixed sections + BCs), Gmsh.

## Phase 9 — CAD import, scripting polish, performance
- Optional STEP/IGES import via OCCT, with healing tolerances and automatic connector/domain creation.
- All golden cases are reproduced from scripts alone (tested).
- Performance targets from SPEC.md measured by a benchmark script in Release on M1. OpenMP/TBB for hot loops, with determinism preserved.
- Packaging: macOS arm64 app bundle; AppImage or conda package for Linux.
