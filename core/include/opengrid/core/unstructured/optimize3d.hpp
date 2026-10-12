// Quality optimisation of mixed 3D meshes (tetrahedra, pyramids, prisms, hexahedra) whose
// boundary and layer nodes are fixed.
//
// Quality: the signed scaled Jacobian of a cell, the minimum over its corners of the triple
// product of the three corner edges divided by their lengths, normalised so the ideal cell scores
// 1 (regular tetrahedron, pyramid with equilateral triangles, right prism on an equilateral
// triangle with square sides, cube) and capped at 1; <= 0 means an inverted corner. Corners as in
// opengrid.model.volumes.CELL_CORNERS (a pyramid's apex is not a corner; it enters through the
// base corners). Cell conventions as in inflation3d.hpp.
//
// Method (L. A. Freitag, C. Ollivier-Gooch, "Tetrahedral mesh improvement using swapping and
// smoothing", Int. J. Numer. Meth. Engng 40 (1997) 3979-4002; worst-cell objective as in
// L. A. Freitag, P. Plassmann, "Local optimization-based simplicial mesh untangling and
// improvement", Int. J. Numer. Meth. Engng 49 (2000) 109-125):
//  * Node optimisation: the movable nodes of cells below the threshold t, worst first, move to
//    lower the sum of squared deficits sum(max(0, t - q)^2) of their cells without lowering
//    their worst quality: to the centre of the other nodes of their cells if that is better,
//    then by a compass search (14 directions, the step halved when no direction improves). A
//    pure max-min objective was tried first; it pulls good neighbours of an unimprovable cell
//    down to its level (more cells below t on a wing's anisotropic layers). Since no move lowers
//    the worst quality around its node, inverted cells are untangled when possible and no valid
//    cell is inverted.
//  * Swapping (tetrahedra from index n_fixed_tets on; others never change): a face shared by two
//    tetrahedra is replaced by three tetrahedra around the edge joining the opposite nodes (2-3),
//    and an interior edge with exactly three tetrahedra and no other cell around it is removed
//    (3-2), when the worst quality of the tetrahedra involved improves.
//  * Passes of both steps alternate until a pass changes nothing or `passes` is reached.
// Deterministic: nodes and swaps are processed in a fixed order (worst quality, then index).
#pragma once

#include <array>
#include <cstdint>
#include <vector>

namespace og {

struct Optimize3DOptions {
  double threshold = 0.3;  // cells below this quality are improved
  int passes = 4;
  bool move_nodes = true;
  bool swap = true;
};

struct Optimize3DInput {
  std::vector<std::array<double, 3>> points;
  std::vector<std::array<std::int64_t, 4>> tets;
  std::vector<std::array<std::int64_t, 5>> pyramids;
  std::vector<std::array<std::int64_t, 6>> prisms;
  std::vector<std::array<std::int64_t, 8>> hexes;
  std::vector<bool> movable;      // per point
  std::int64_t n_fixed_tets = 0;  // the first tetrahedra are never swapped (e.g. layer cells)
};

struct Optimize3DResult {
  std::vector<std::array<double, 3>> points;
  std::vector<std::array<std::int64_t, 4>> tets;  // the fixed ones first, unchanged
  std::int64_t moved = 0;                         // accepted node moves
  std::int64_t swaps_23 = 0;
  std::int64_t swaps_32 = 0;
  double worst_before = 1.0;
  double worst_after = 1.0;
  std::int64_t below_before = 0;  // cells below the threshold
  std::int64_t below_after = 0;
  std::int64_t inverted_after = 0;  // cells with quality <= 0
};

// Signed scaled Jacobian of a cell with n = 4 (tetrahedron), 5 (pyramid), 6 (prism) or 8
// (hexahedron) nodes.
double cell_quality_3d(const std::vector<std::array<double, 3>>& points, const std::int64_t* nodes,
                       int n);

// Throws std::invalid_argument for out-of-range node ids or sizes that do not match.
Optimize3DResult optimize_3d(const Optimize3DInput& input, const Optimize3DOptions& options);

}  // namespace og
