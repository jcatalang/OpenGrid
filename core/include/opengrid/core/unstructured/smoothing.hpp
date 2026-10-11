// Node smoothing and untangling of planar triangle / quadrilateral meshes with fixed boundary
// nodes (the first n_fixed points).
//
// Smoothing methods (one Gauss-Seidel sweep over the free nodes per iteration, in index order):
//  * laplacian        node -> mean of its edge neighbours;
//  * length_weighted  node -> mean of its edge neighbours weighted by the edge lengths, which
//                     pulls nodes towards distant neighbours and evens out edge lengths;
//  * angle_based      T. Zhou, K. Shimada, "An angle-based approach to two-dimensional mesh
//                     smoothing", 9th Int. Meshing Roundtable (2000): every edge neighbour N_j
//                     proposes the point on the bisector of the angle formed at N_j by its two
//                     ring neighbours, at the current distance |V - N_j|; the node moves to the
//                     mean of the proposals.
// Every move is guarded: it is undone if it lowers the worst scaled Jacobian of the cells around
// the node (so smoothing never inverts a cell or makes the worst cell worse). With edge_swap,
// triangle pairs are flipped to the Delaunay configuration after every sweep (boundary
// segments are never flipped; quads are left alone).
//
// Untangling / optimisation: nodes of cells with a scaled Jacobian at or below a threshold (0:
// inverted or non-convex) are moved by a compass search that maximises the worst signed scaled
// Jacobian of their cells (optimisation-based untangling, cf. P. M. Knupp, Eng. Comput. 15 (1999)
// 263-268, with the worst-cell objective instead of a sum).
#pragma once

#include <array>
#include <cstdint>
#include <vector>

namespace og {

enum class SmoothingMethod { none, laplacian, length_weighted, angle_based };

struct SmoothingOptions {
  SmoothingMethod method = SmoothingMethod::none;
  int iterations = 5;
  bool edge_swap = true;
};

// Signed scaled Jacobian of a cell: min over corners of sin(angle) (divided by sin 60 deg for
// triangles); 1 is ideal, <= 0 means inverted or non-convex.
double cell_quality(const std::vector<std::array<double, 2>>& points,
                    const std::array<std::int64_t, 4>& corners, int n);

// Returns the number of accepted node moves.
std::int64_t smooth_mesh(std::vector<std::array<double, 2>>& points,
                         std::vector<std::array<std::int64_t, 3>>& triangles,
                         const std::vector<std::array<std::int64_t, 4>>& quads,
                         const std::vector<std::array<std::int64_t, 2>>& segments,
                         std::size_t n_fixed, const SmoothingOptions& options);

// Nodes of cells whose signed scaled Jacobian is <= threshold are optimised (threshold 0:
// untangling; a positive threshold also improves poor cells). Returns the number of cells that
// are still inverted or non-convex afterwards.
std::int64_t untangle(std::vector<std::array<double, 2>>& points,
                      const std::vector<std::array<std::int64_t, 3>>& triangles,
                      const std::vector<std::array<std::int64_t, 4>>& quads, std::size_t n_fixed,
                      double threshold = 0.0);

}  // namespace og
