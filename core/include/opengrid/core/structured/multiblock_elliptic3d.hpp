// Elliptic smoothing of a multi-block structured 3D grid as one grid: block
// interfaces, interface edges and block corners move; boundary nodes slide on their surfaces
// and curves.
//
// Method: the 3D form of structured/multiblock_elliptic.hpp. The grid is a set of hexahedra
// over one global node array. Around a node the incident hexahedra are assembled into octants
// of a local frame by walking across shared faces (a hexahedron's three edges from the node map
// to signed frame axes; the hexahedron across a face keeps the two in-face edges and reverses
// the third), starting from the hexahedron in which the node is the first corner (a single
// block keeps its (i, j, k) frame). A node with eight hexahedra filling the eight octants is
// *regular* and gets a 3 x 3 x 3 stencil; it satisfies the 3D Laplace / Poisson system of
// the single-block smoother (structured/elliptic3d.hpp: cofactors of the covariant metric,
// central differences, wall sources P_m) in its frame; its metric coefficients are computed
// once per iteration (frozen for the three line sweeps through it), the cross-derivative terms
// with the current positions. Grid lines are traced through regular
// free nodes along the three frame axes across blocks and relaxed implicitly; lines are
// coloured greedily over their 3 x 3 x 3 neighbourhoods and lines of one colour are solved in
// parallel (thread-count independent).
//  * Irregular interior nodes (along edges shared by 3 or 5+ blocks, block corners of odd
//    valence): mean of the average of the edge neighbours and the average of the trilinear
//    (parallelepiped) completions a + b + c - ab - bc - ca + abc of the incident hexahedra.
//  * Sliding nodes on a boundary face (four hexahedra in a half space) move towards the
//    projection of their inward neighbour onto the face's surface; on a boundary edge (two
//    hexahedra) along the edge's curve chain towards the projection of the inward diagonal
//    neighbour; under-relaxed by 1/2. The step from the initial (u, v) or curve parameter is
//    limited in parameter space (with the metric at the initial point) to slide_limit times
//    the smallest initial distance to the boundary neighbours: the limited point lies on the
//    geometry and limiting it again changes nothing, which keeps the iteration convergent.
//  * Wall nodes (fixed, on a boundary face): Hilgenstock feedback sources (normal from the
//    first-cell height, two tangential from the lean of the first grid line), carried along
//    the grid line leaving the wall with weight exp(-decay * lines), rotated into each frame.
// Anderson acceleration over positions and wall sources; sliding nodes are put back on their
// surface or curve after each extrapolation.
// References: J. F. Thompson, Z. U. A. Warsi, C. W. Mastin, Numerical Grid Generation (1985),
// ch. VI; P. R. Eiseman, NASA CP-3291 (1995); A. Hilgenstock (1988); W. A. Spekreijse,
// J. Comput. Phys. 118 (1995) 38-61.
#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "opengrid/core/geom/curve_chain.hpp"
#include "opengrid/core/geom/surface.hpp"
#include "opengrid/core/structured/multiblock_elliptic.hpp"

namespace og {

enum class SlideTarget : std::uint8_t { Surface = 0, Chain = 1 };

struct Multiblock3DSlide {
  std::int64_t node = 0;
  SlideTarget target = SlideTarget::Surface;
  std::size_t index = 0;  // into surfaces or chains
};

struct Multiblock3DInput {
  std::vector<double> xyz;  // 3 per node
  // Corner c = i + 2 j + 4 k of each hexahedron, right-handed (positive volume).
  std::vector<std::array<std::int64_t, 8>> hexes;
  std::vector<std::uint8_t> fixed;
  std::vector<std::shared_ptr<const Surface>> surfaces;
  std::vector<std::shared_ptr<const CurveChain>> chains;
  std::vector<Multiblock3DSlide> slides;
  std::vector<MultiblockWall> walls;  // nodes on boundary faces
  // Per node (empty: none): nodes of symmetry planes. A face node sliding on a plane and
  // flagged here is updated like an interior node with its neighbours mirrored across the
  // plane, not limited; edges between two such planes likewise, and an edge between one and
  // another surface slides to the foot of its neighbour in the plane (a 2D boundary node).
  std::vector<std::uint8_t> symmetric;
};

struct Multiblock3DResult {
  std::vector<double> xyz;
  std::vector<std::uint8_t> node_class;  // MultiblockNode
  std::vector<double> residuals;
  int iterations = 0;
  bool converged = false;
  bool cancelled = false;
  int lines = 0;
  int colours = 0;
  std::vector<std::string> warnings;
};

// Throws std::invalid_argument for malformed input.
Multiblock3DResult smooth_multiblock_3d(const Multiblock3DInput& input,
                                        const MultiblockOptions& options,
                                        SmoothingProgress* progress = nullptr);

}  // namespace og
