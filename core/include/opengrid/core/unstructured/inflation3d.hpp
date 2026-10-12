// 3D boundary-layer inflation: prism / hexahedron layers advanced from wall faces of a closed
// surface (the 3D counterpart of inflation.hpp).
//
// Method (advancing layers: Y. Kallinderis, S. Ward, AIAA J. 31 (1993) 1440-1447;
// S. Pirzadeh, AIAA J. 34 (1996) 43-49; with local termination of each column):
//  * The surface is closed and its faces (triangles / quads) are oriented out of the region;
//    the wall faces are inflated into the region. Each wall node carries a column. A node whose
//    faces are all walls is free; a node that also touches non-wall faces of one plane (group
//    `node_group`, e.g. a symmetry plane) slides in that plane; any other node is fixed (layers
//    taper to zero towards it, as at open wall ends in 2D). A sliding node whose group faces are
//    not coplanar is fixed (with a warning).
//  * Layer k has height h_k = first_height * growth_rate^k. Every column advances along its
//    marching normal: the angle-weighted mean of the inward normals of its front faces (projected
//    into the plane for sliding nodes), smoothed over the front (Jacobi sweeps accepted while the
//    normal stays visible from every adjacent face, n . n_f >= min_visibility), scaled by
//    1 / min(n . n_f) (at most max_scale) so the layer keeps the height h_k from every adjacent
//    face. Normals are recomputed from the current front for every layer.
//  * A column's isotropic size is the mean length of its front edges (at least the wall spacing
//    it started from), capped by max_size. A column stops for good when: max_layers (or
//    total_thickness) is reached; with transition, when h_k reaches transition_ratio times that
//    size; when it would stand more than max_slope times the lateral gap above a stopped
//    neighbour, and more than one layer; when a facing front or boundary face lies ahead closer
//    than proximity times its size; when its new column edge crosses a front or boundary face,
//    or a new front face crosses one (floating-point tests with a small tolerance); when its
//    marching normal is not visible from all adjacent faces; or when a cell of its layer would be
//    inverted (or, for full prisms / hexahedra, have a corner scaled Jacobian below
//    min_jacobian; then the face's column with the most tilted normal stops first). These checks
//    are repeated until no column has to stop, so every committed layer is valid.
//  * Cells: prisms on triangles, hexahedra on quads. Where some columns of a face stop, the cell
//    degenerates: a triangle gives a pyramid (one stopped) or a tetrahedron (two); a quad gives a
//    pyramid (three stopped) or else a transition cell split from its centroid into pyramids
//    (quad faces) and tetrahedra (triangle faces), and from then on the face continues as two
//    triangles (so the front has no steep quads at steps). All shared faces are kept exactly,
//    so the layers are conforming.
// The final front (one or two polygons per wall face, same orientation as the wall face) and the
// non-wall faces enclose the core region, which the caller fills with tetrahedra. Lateral faces
// of the layers on sliding planes are returned with the non-wall face they lie on, so the caller
// can re-mesh those faces around them.
#pragma once

#include <array>
#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace og {

struct Inflation3DOptions {
  double first_height = 1e-3;
  double growth_rate = 1.2;
  int max_layers = 30;
  std::optional<double> total_thickness;
  bool transition = true;
  double transition_ratio = 1.0;
  std::optional<double> max_size;
  double proximity = 1.0;
  double min_jacobian = 0.2;
  double max_slope = 0.5;
  int normal_smoothing = 5;
  double max_scale = 2.0;
  double min_visibility = 0.2;
};

struct Inflation3DInput {
  std::vector<std::array<double, 3>> points;
  std::vector<std::vector<std::int64_t>> faces;  // triangles / quads, oriented out of the region
  std::vector<bool> wall;                        // per face
  std::vector<std::int64_t> face_group;          // per face (e.g. its domain); >= 0
  std::vector<std::int64_t> node_group;          // per point: plane group it may slide in, or -1
};

struct Inflation3DResult {
  std::vector<std::array<double, 3>> points;  // input points first, then layer nodes
  std::vector<std::array<std::int64_t, 4>> tets;
  std::vector<std::array<std::int64_t, 5>> pyramids;  // base normal towards the apex
  std::vector<std::array<std::int64_t, 6>> prisms;    // first triangle's normal towards the second
  std::vector<std::array<std::int64_t, 8>> hexes;     // first quad's normal towards the second
  std::vector<std::vector<std::int64_t>> front;       // final front polygons
  std::vector<std::int64_t> front_face;               // the wall face of each front polygon
  std::vector<std::vector<std::int64_t>> sides;       // lateral faces on sliding planes
  std::vector<std::int64_t> side_face;                // the non-wall face each lies on
  std::vector<std::int64_t> top;                      // per input point: its final column node
  std::vector<double> height;                         // per input point: nominal layer height
  std::vector<int> layers;                            // per input point: number of layers
  std::vector<std::string> reason;  // per input point: why its column stopped ("": not a wall node)
  int max_layers = 0;
  std::map<std::string, std::int64_t> stops;  // columns stopped, by reason
  std::vector<std::string> warnings;
};

// Throws std::invalid_argument for invalid options, malformed faces (not 3 or 4 nodes, indices
// out of range, sizes not matching) or no wall face.
Inflation3DResult inflate_3d(const Inflation3DInput& input, const Inflation3DOptions& options);

}  // namespace og
