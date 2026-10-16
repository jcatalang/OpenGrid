// Wall layers: the first N grid layers off a wall marched along the normals of the layer below
// (a boundary-layer block inside a multi-block grid), so that grid lines stay normal to the
// wall over N layers, not only at the first cell.
//
// Method: algebraic normal marching, as the boundary-layer stage of hybrid hyperbolic /
// elliptic generation. The grid lines leaving the wall are found in the cell complex (quads or
// hexahedra over one node array): from a wall node, every edge to a free node that is not a
// wall node starts a line, and a line continues through a node to the one edge neighbour that
// shares no cell with the edge it arrived by (the opposite node of a regular node, whatever the
// blocks' frames). Layer k is placed from layer k - 1: the direction of a line is the sum of
// the unit normals of the cell facets at its layer k - 1 node that lie in layer k - 1 (for each
// cell around the line's edge, the facet at that node across the edge), so it bisects the
// front's angles (e.g. a sharp trailing edge, or where a wall meets an internal wall); it is
// blended with its neighbours' directions (weight 1/10, the dissipation of hyperbolic marching)
// and kept in a symmetry plane where the node lies on one. Heights: spacing * growth^(k - 1).
// The nodes beyond the last layer are shifted along each line by the last layer's displacement,
// decaying linearly in initial arc length to zero at the line's end, as the starting grid for
// the elliptic smoothing (which then controls the first outer cell: height spacing * growth^N,
// orthogonal to the last layer). Layers that would overlap, lines that end before N layers and
// cells of the layers that fold are reported (exceptions), never clipped.
// References: J. L. Steger, D. S. Chaussee, "Generation of body-fitted coordinates using
// hyperbolic partial differential equations", SIAM J. Sci. Stat. Comput. 1 (1980) 431-437
// (marching off a wall along the front normals); W. M. Chan, "Hyperbolic methods for surface
// and field grid generation", ch. 5 of J. F. Thompson, B. K. Soni, N. P. Weatherill (eds.),
// Handbook of Grid Generation (1999) (smoothing of the marching directions).
#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "opengrid/core/geom/vec3.hpp"

namespace og {

struct WallLayerNode {
  std::int64_t node = 0;
  double spacing = 0.0;  // first height; <= 0: the initial distance to the first line node
  int layers = 0;
  double growth = 1.2;
  std::size_t source = 0;  // the caller's index (e.g. of its wall entry)
};

struct WallLayerPlane {
  Vec3 point;
  Vec3 normal;  // unit
};

struct WallLayerInput {
  int dim = 2;  // 2: quads (counter-clockwise, first four entries); 3: hexahedra (corner
                // c = i + 2 j + 4 k, right-handed)
  std::vector<std::array<std::int64_t, 8>> cells;
  std::vector<std::uint8_t> movable;  // per node: free nodes the layers may move
  std::vector<std::uint8_t> on_wall;  // per node: every wall node (with layers or without)
  std::vector<std::int32_t> plane;    // per node (empty: none): index into planes, or -1
  std::vector<WallLayerPlane> planes;
  std::vector<WallLayerNode> walls;  // wall nodes with layers > 0
};

struct WallLayerLine {
  std::int64_t base = 0;            // the wall node
  std::size_t source = 0;           // WallLayerNode::source
  std::vector<std::int64_t> nodes;  // layers 1..N
  double next_spacing = 0.0;        // height of the first cell beyond the layers
};

struct WallLayerResult {
  std::vector<WallLayerLine> lines;
};

// Moves the layer nodes in x (and shifts the free nodes beyond them along their lines). Throws
// std::invalid_argument (with the position) when the layers of two lines overlap, a line ends
// before its layers, or a cell of the layers folds.
WallLayerResult march_wall_layers(const WallLayerInput& input, std::vector<Vec3>& x);

}  // namespace og
