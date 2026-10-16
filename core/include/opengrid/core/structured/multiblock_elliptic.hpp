// Elliptic smoothing of a planar multi-block structured grid as one grid:
// block interfaces, block corners and boundary nodes move.
//
// Method: the grid is a set of quadrilaterals over one global node array (shared nodes merged
// by the caller). Around every node the incident quads are chained counter-clockwise; a node
// with a closed chain of four quads is *regular* and gets a 3 x 3 stencil in a local frame
// (E, N, W, S and the diagonals) assembled from the quads, whichever domains they come from.
// Regular free nodes satisfy the Poisson equation of the single-block smoother
// (structured/elliptic.hpp, Laplace control plus wall sources) in their local frame, so grid
// lines continue smoothly across interfaces. Grid lines are traced through regular free nodes
// across domains and relaxed implicitly (tridiagonal solves, coefficients frozen per sweep);
// lines end at fixed, sliding or irregular nodes, and a closed line is opened at its first
// node. Lines are coloured greedily (lines of one colour never read each other's nodes) and
// lines of one colour are solved in parallel, so the result does not depend on the thread
// count; in a single block the colouring is the zebra ADI ordering.
//  * Irregular nodes (3 or 5+ quads around an interior node, e.g. the corners of an O-grid
//    core) have no (i, j) stencil: they move to the mean of the average of their edge
//    neighbours and the average of the parallelogram completions p + m - d of their quads
//    (both exact on a uniform rectangular grid; cf. the treatment of irregular points by
//    W. A. Spekreijse, J. Comput. Phys. 118 (1995) 38-61, which uses the local quads).
//  * Sliding boundary nodes move along their boundary curve chain to the foot of the grid line
//    leaving the boundary (orthogonal projection of the inward neighbour), under-relaxed by 1/2
//    (a Neumann/slip condition); each stays between its boundary neighbours and within half
//    its initial intervals of its initial position (so nodes cannot pile up into a corner).
//  * Wall control (fixed wall nodes): Hilgenstock feedback sources as in the single-block
//    smoother, p (tangential) and q (normal) per wall node, carried along the grid line that
//    leaves the wall (through regular free nodes, across interfaces) with weight
//    exp(-decay * lines from the wall), rotated into each node's frame.
// The iteration (line sweeps, point updates, wall feedback) is accelerated by Anderson mixing
// (util/anderson.hpp) over node positions and wall sources; sliding nodes are put back on
// their curve after each extrapolation.
// References: J. F. Thompson, Z. U. A. Warsi, C. W. Mastin, Numerical Grid Generation (1985),
// ch. VI (multiblock grids); P. R. Eiseman, "Multiblock grid generation with
// automatic zoning", NASA CP-3291 (1995); A. Hilgenstock, Numerical Grid Generation in CFD '88
// (1988) 137-146; W. A. Spekreijse (1995), above.
#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "opengrid/core/geom/curve_chain.hpp"
#include "opengrid/core/structured/elliptic.hpp"

namespace og {

struct MultiblockSlide {
  std::int64_t node = 0;
  std::size_t chain = 0;  // index into MultiblockInput::chains
};

struct MultiblockWall {
  std::int64_t node = 0;
  double spacing = 0.0;  // first-cell height; <= 0: the initial height at this node
  double decay = 0.1;
  double relaxation = 0.1;
  double max_source = 100.0;
};

struct MultiblockInput {
  std::vector<double> xy;  // 2 per node
  double z = 0.0;
  std::vector<std::array<std::int64_t, 4>> quads;  // counter-clockwise in the xy plane
  std::vector<std::uint8_t> fixed;                 // per node: 1 = never moves
  std::vector<std::shared_ptr<const CurveChain>> chains;
  std::vector<MultiblockSlide> slides;
  std::vector<MultiblockWall> walls;
};

struct MultiblockOptions {
  double omega = 1.0;
  int max_iterations = 2000;
  double tolerance = 1e-6;
  int acceleration = 20;  // Anderson depth (0: plain relaxation)
  // Sliding nodes stay within this fraction of their initial boundary intervals of their
  // initial position (3D; the 2D smoother uses 1/2 of the intervals towards each neighbour).
  double slide_limit = 0.5;
};

// Node classes reported per node.
enum class MultiblockNode : std::uint8_t { Fixed = 0, Regular = 1, Irregular = 2, Sliding = 3 };

struct MultiblockResult {
  std::vector<double> xy;
  std::vector<std::uint8_t> node_class;  // MultiblockNode
  std::vector<double> residuals;
  int iterations = 0;
  bool converged = false;
  bool cancelled = false;
  int lines = 0;
  int colours = 0;
  std::vector<std::string> warnings;
};

// Throws std::invalid_argument for malformed input (indices out of range, a node shared by
// quads that do not chain around it).
MultiblockResult smooth_multiblock(const MultiblockInput& input, const MultiblockOptions& options,
                                   SmoothingProgress* progress = nullptr);

}  // namespace og
