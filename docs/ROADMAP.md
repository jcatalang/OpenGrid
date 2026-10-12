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
- Python scripting API documentation and examples reproducing all golden cases from scripts only
  (tested).
- Domain/block splitting and joining at given indices; O-grid / butterfly topology helpers for
  multi-block hexahedral meshes.

### 9c — Performance and packaging
- Performance targets from SPEC.md measured by a benchmark script in Release on M1. OpenMP/TBB
  for hot loops, with determinism preserved.
- Packaging: macOS arm64 app bundle; AppImage or conda package for Linux.
