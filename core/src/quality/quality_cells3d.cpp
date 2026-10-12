#include "opengrid/core/quality/quality_cells3d.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numbers>
#include <stdexcept>
#include <utility>
#include <vector>

#include "opengrid/core/geom/vec3.hpp"

namespace og {

namespace {

// Corners (node, three neighbours with a positive triple product in a valid cell), the factor
// scaling the ideal cell's corner value to 1, outward faces, and the ideal smallest / largest
// dihedral angles, per cell size.
struct CellType {
  std::vector<std::array<int, 4>> corners;
  double factor = 1.0;
  std::vector<std::vector<int>> faces;
  double ideal_lo = 90.0;
  double ideal_hi = 90.0;
  std::vector<std::pair<int, int>> edge_faces;  // the two faces at every edge
};

CellType make(std::vector<std::array<int, 4>> corners, double factor,
              std::vector<std::vector<int>> faces, double lo, double hi) {
  CellType t{std::move(corners), factor, std::move(faces), lo, hi, {}};
  // Faces traversing an edge in opposite directions share it.
  std::map<std::pair<int, int>, int> owner;
  for (std::size_t f = 0; f < t.faces.size(); ++f) {
    const auto& face = t.faces[f];
    for (std::size_t i = 0; i < face.size(); ++i) {
      owner[{face[i], face[(i + 1) % face.size()]}] = static_cast<int>(f);
    }
  }
  for (const auto& [edge, f] : owner) {
    if (edge.first < edge.second) {
      t.edge_faces.emplace_back(f, owner.at({edge.second, edge.first}));
    }
  }
  return t;
}

const CellType& cell_type(int n) {
  const double tet_angle = std::acos(1.0 / 3.0) * 180.0 / std::numbers::pi;
  static const CellType tet =
      make({{0, 1, 2, 3}, {1, 2, 0, 3}, {2, 0, 1, 3}, {3, 0, 2, 1}}, std::numbers::sqrt2,
           {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {0, 3, 2}}, tet_angle, tet_angle);
  static const CellType pyramid =
      make({{0, 1, 3, 4}, {1, 2, 0, 4}, {2, 3, 1, 4}, {3, 0, 2, 4}}, std::numbers::sqrt2,
           {{0, 3, 2, 1}, {0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}},
           std::atan(std::numbers::sqrt2) * 180.0 / std::numbers::pi, 180.0 - tet_angle);
  static const CellType prism =
      make({{0, 1, 2, 3}, {1, 2, 0, 4}, {2, 0, 1, 5}, {3, 5, 4, 0}, {4, 3, 5, 1}, {5, 4, 3, 2}},
           2.0 / std::numbers::sqrt3,
           {{0, 2, 1}, {3, 4, 5}, {0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}}, 60.0, 90.0);
  static const CellType hex = make(
      {{0, 1, 3, 4},
       {1, 2, 0, 5},
       {2, 3, 1, 6},
       {3, 0, 2, 7},
       {4, 7, 5, 0},
       {5, 4, 6, 1},
       {6, 5, 7, 2},
       {7, 6, 4, 3}},
      1.0, {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}},
      90.0, 90.0);
  if (n == 4) {
    return tet;
  }
  if (n == 5) {
    return pyramid;
  }
  if (n == 6) {
    return prism;
  }
  if (n == 8) {
    return hex;
  }
  throw std::invalid_argument("cell_quality_3d: a cell has 4, 5, 6 or 8 nodes");
}

Vec3 at(const std::vector<std::array<double, 3>>& points, std::int64_t i) {
  const auto& p = points[static_cast<std::size_t>(i)];
  return {p[0], p[1], p[2]};
}

double scaled_jacobian(const std::vector<std::array<double, 3>>& points, const std::int64_t* nodes,
                       const CellType& type) {
  double worst = 1.0;
  for (const auto& c : type.corners) {
    const Vec3 o = at(points, nodes[c[0]]);
    const Vec3 e1 = at(points, nodes[c[1]]) - o;
    const Vec3 e2 = at(points, nodes[c[2]]) - o;
    const Vec3 e3 = at(points, nodes[c[3]]) - o;
    const double lengths = norm(e1) * norm(e2) * norm(e3);
    if (lengths <= 0.0) {
      return -1.0;
    }
    worst = std::min(worst, dot(cross(e1, e2), e3) / lengths * type.factor);
  }
  return worst;
}

}  // namespace

double cell_quality_3d(const std::vector<std::array<double, 3>>& points, const std::int64_t* nodes,
                       int n) {
  return scaled_jacobian(points, nodes, cell_type(n));
}

CellQuality3D cell_quality_3d_full(const std::vector<std::array<double, 3>>& points,
                                   const std::int64_t* nodes, int n) {
  const CellType& type = cell_type(n);
  CellQuality3D q;
  q.scaled_jacobian = scaled_jacobian(points, nodes, type);
  std::vector<Vec3> normal;
  normal.reserve(type.faces.size());
  for (const auto& face : type.faces) {
    const Vec3 p0 = at(points, nodes[face[0]]);
    const Vec3 p1 = at(points, nodes[face[1]]);
    const Vec3 p2 = at(points, nodes[face[2]]);
    if (face.size() == 3) {
      normal.push_back(cross(p1 - p0, p2 - p0) * 0.5);
    } else {
      normal.push_back(cross(p2 - p0, at(points, nodes[face[3]]) - p1) * 0.5);
    }
  }
  q.min_angle = 180.0;
  q.max_angle = 0.0;
  for (const auto& [a, b] : type.edge_faces) {
    const double scale = std::max(
        norm(normal[static_cast<std::size_t>(a)]) * norm(normal[static_cast<std::size_t>(b)]),
        1e-300);
    const double cosine = std::clamp(
        dot(normal[static_cast<std::size_t>(a)], normal[static_cast<std::size_t>(b)]) / scale, -1.0,
        1.0);
    const double dihedral = 180.0 - std::acos(cosine) * 180.0 / std::numbers::pi;
    q.min_angle = std::min(q.min_angle, dihedral);
    q.max_angle = std::max(q.max_angle, dihedral);
  }
  q.skewness = std::max({(q.max_angle - type.ideal_hi) / (180.0 - type.ideal_hi),
                         (type.ideal_lo - q.min_angle) / type.ideal_lo, 0.0});
  return q;
}

}  // namespace og
