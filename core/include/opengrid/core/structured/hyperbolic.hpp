// Hyperbolic extrusion: a block grown layer by layer from a structured surface grid (the k = 0
// face).
//
// Method: explicit marching on the hyperbolic grid conditions (grid lines orthogonal to the
// layers, prescribed layer heights) - the marching direction of a node is the normal of the
// current layer (central differences of r_xi x r_eta), smoothed a few Jacobi sweeps along the
// layer, and the step is the layer height divided by the cosine between that direction and the
// local face normals (so convex and concave corners keep the height); explicit second-order
// dissipation then moves the new nodes towards the mean of their neighbours in the layer, with
// the component along the marching direction removed (it smooths the spacing along the layer
// without changing heights), as in the dissipation of the hyperbolic schemes of Steger &
// Chaussee (SIAM J. Sci. Stat. Comput. 1 (1980) 431-437) and Chan & Steger (Appl. Math. Comput.
// 51 (1992) 181-205). This is an explicit variant, not their implicit block-tridiagonal solver.
// Edges of the base grid: free (one-sided differences), constant x / y / z (the edge nodes stay
// in the plane of their base node, e.g. symmetry planes), or periodic (i or j closed: the first
// and last lines coincide).
// The block is right-handed when marching to the side of x_i x x_j (negative distances march
// to the other side with i mirrored).
#pragma once

#include <array>
#include <string>
#include <vector>

#include "opengrid/core/structured/grid3d.hpp"

namespace og {

enum class HyperbolicEdge { Free, ConstantX, ConstantY, ConstantZ, Periodic };

struct HyperbolicOptions {
  std::array<HyperbolicEdge, 4> edges = {HyperbolicEdge::Free, HyperbolicEdge::Free,
                                         HyperbolicEdge::Free, HyperbolicEdge::Free};
  double dissipation = 0.2;  // per layer, 0..1 (fraction of the move to the neighbour mean)
  int normal_smoothing = 4;  // Jacobi sweeps of the marching directions per layer
};

struct HyperbolicResult {
  Grid3D grid;
  bool mirrored = false;
  int inverted_layers = 0;  // layers with a non-positive cell Jacobian
  std::vector<std::string> warnings;
};

// fractions: s_0 = 0 < ... < s_{nk-1} = 1 of the total `distance`. Throws std::invalid_argument
// for invalid fractions, a zero distance, or periodic edges whose lines do not coincide.
HyperbolicResult extrude_hyperbolic(const Grid2D& base, double distance,
                                    const std::vector<double>& fractions,
                                    const HyperbolicOptions& options);

}  // namespace og
