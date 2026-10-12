#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <vector>

#include "opengrid/core/unstructured/inflation3d.hpp"

namespace {

using Index = std::int64_t;
using P3 = std::array<double, 3>;

std::size_t u(Index i) {
  return static_cast<std::size_t>(i);
}

// Box surfaces with n cells per edge (quads, or each quad split into two triangles), oriented
// with sign * outward normals; points merged by coordinates.
struct Builder {
  og::Inflation3DInput in;
  std::map<std::array<long long, 3>, Index> ids;

  Index point(const P3& p) {
    const std::array<long long, 3> key = {std::llround(p[0] * 1e9), std::llround(p[1] * 1e9),
                                          std::llround(p[2] * 1e9)};
    auto it = ids.find(key);
    if (it != ids.end()) {
      return it->second;
    }
    const auto id = static_cast<Index>(in.points.size());
    in.points.push_back(p);
    in.node_group.push_back(-1);
    ids[key] = id;
    return id;
  }

  // Returns the indices of the added faces, in order axis, side.
  std::vector<std::vector<std::size_t>> box(const P3& lo, const P3& hi, int n, double sign,
                                            bool triangles, bool wall) {
    std::vector<std::vector<std::size_t>> per_face;
    for (int axis = 0; axis < 3; ++axis) {
      for (int side = 0; side < 2; ++side) {
        const int a = (axis + 1) % 3;
        const int b = (axis + 2) % 3;
        std::vector<std::size_t> added;
        for (int i = 0; i < n; ++i) {
          for (int j = 0; j < n; ++j) {
            std::array<Index, 4> q{};
            const std::array<std::array<int, 2>, 4> corners = {
                {{i, j}, {i + 1, j}, {i + 1, j + 1}, {i, j + 1}}};
            for (std::size_t k = 0; k < 4; ++k) {
              P3 p{};
              p[u(axis)] = lo[u(axis)];
              if (side == 1) {
                p[u(axis)] = hi[u(axis)];
              }
              p[u(a)] = lo[u(a)] + (hi[u(a)] - lo[u(a)]) * corners[k][0] / n;
              p[u(b)] = lo[u(b)] + (hi[u(b)] - lo[u(b)]) * corners[k][1] / n;
              q[k] = point(p);
            }
            // (a, b, axis) is right-handed, so (a, b) order has normal +axis.
            const bool flip = (side == 0) == (sign > 0);
            std::vector<std::vector<Index>> cells;
            if (triangles) {
              cells = {{q[0], q[1], q[2]}, {q[0], q[2], q[3]}};
            } else {
              cells = {{q[0], q[1], q[2], q[3]}};
            }
            for (auto cell : cells) {
              if (flip) {
                std::reverse(cell.begin(), cell.end());
              }
              added.push_back(in.faces.size());
              in.faces.push_back(cell);
              in.wall.push_back(wall);
              in.face_group.push_back(axis * 2 + side + static_cast<int>(wall) * 10);
            }
          }
        }
        per_face.push_back(added);
      }
    }
    return per_face;
  }
};

double det(const P3& o, const P3& a, const P3& b, const P3& c) {
  const P3 e1 = {a[0] - o[0], a[1] - o[1], a[2] - o[2]};
  const P3 e2 = {b[0] - o[0], b[1] - o[1], b[2] - o[2]};
  const P3 e3 = {c[0] - o[0], c[1] - o[1], c[2] - o[2]};
  return e1[0] * (e2[1] * e3[2] - e2[2] * e3[1]) - e1[1] * (e2[0] * e3[2] - e2[2] * e3[0]) +
         e1[2] * (e2[0] * e3[1] - e2[1] * e3[0]);
}

void expect_positive(const og::Inflation3DResult& r) {
  const auto& p = r.points;
  for (const auto& t : r.tets) {
    EXPECT_GT(det(p[u(t[0])], p[u(t[1])], p[u(t[2])], p[u(t[3])]), 0.0);
  }
  for (const auto& y : r.pyramids) {
    for (int k = 0; k < 4; ++k) {
      EXPECT_GT(det(p[u(y[u(k)])], p[u(y[u((k + 1) % 4)])], p[u(y[u((k + 3) % 4)])], p[u(y[4])]),
                0.0);
    }
  }
  for (const auto& w : r.prisms) {
    EXPECT_GT(det(p[u(w[0])], p[u(w[1])], p[u(w[2])], p[u(w[3])]), 0.0);
    EXPECT_GT(det(p[u(w[3])], p[u(w[5])], p[u(w[4])], p[u(w[0])]), 0.0);
  }
  for (const auto& h : r.hexes) {
    EXPECT_GT(det(p[u(h[0])], p[u(h[1])], p[u(h[3])], p[u(h[4])]), 0.0);
    EXPECT_GT(det(p[u(h[6])], p[u(h[5])], p[u(h[7])], p[u(h[2])]), 0.0);
  }
}

// Every cell meets the quality floor (scaled Jacobian and skewness).
void expect_floor(const og::Inflation3DResult& r, const og::QualityFloor& floor = {}) {
  auto check = [&](const auto& cells, int n) {
    for (const auto& c : cells) {
      const og::CellQuality3D q = og::cell_quality_3d_full(r.points, c.data(), n);
      EXPECT_TRUE(floor.passes(q)) << "cell of " << n << " nodes: scaled Jacobian "
                                   << q.scaled_jacobian << ", skewness " << q.skewness;
    }
  };
  check(r.tets, 4);
  check(r.pyramids, 5);
  check(r.prisms, 6);
  check(r.hexes, 8);
}

// No floor beyond valid cells (to test the stopping mechanics on their own).
og::QualityFloor no_floor() {
  og::QualityFloor f;
  f.min_jacobian = -1.0;
  f.max_skewness = 1.0;
  return f;
}

double dist(const P3& a, const P3& b) {
  return std::hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

og::Inflation3DOptions options(int layers) {
  og::Inflation3DOptions o;
  o.first_height = 0.01;
  o.growth_rate = 1.2;
  o.max_layers = layers;
  o.transition = false;
  return o;
}

// A prism of diamond section (sharp edges at x = 0 and x = 2, wedge half angle atan(t / 1)) from
// z = 0 to z = 1, all of it wall, faces oriented into the body (out of the region).
void diamond(Builder& b, double t, int nx, int nz) {
  auto y_of = [&](int i, double sign) {
    return sign * t * (1.0 - std::abs(i - nx) / static_cast<double>(nx));
  };
  // Section loop counter-clockwise seen from +z: lower chain x = 0..2, then upper chain back.
  std::vector<std::array<double, 2>> loop;
  loop.reserve(4 * static_cast<std::size_t>(nx));
  for (int i = 0; i < 2 * nx; ++i) {
    loop.push_back({static_cast<double>(i) / nx, y_of(i, -1.0)});
  }
  for (int i = 2 * nx; i > 0; --i) {
    loop.push_back({static_cast<double>(i) / nx, y_of(i, 1.0)});
  }
  const std::size_t n = loop.size();
  auto add = [&](const std::vector<Index>& cell) {
    b.in.faces.push_back(cell);
    b.in.wall.push_back(true);
    b.in.face_group.push_back(100);
  };
  for (int k = 0; k < nz; ++k) {
    const double z0 = static_cast<double>(k) / nz;
    const double z1 = static_cast<double>(k + 1) / nz;
    for (std::size_t i = 0; i < n; ++i) {
      const auto& p = loop[i];
      const auto& q = loop[(i + 1) % n];
      // Outward (out of the body) is (q - p) x z for a counter-clockwise loop; into the body:
      // reverse.
      add({b.point({p[0], p[1], z0}), b.point({p[0], p[1], z1}), b.point({q[0], q[1], z1}),
           b.point({q[0], q[1], z0})});
    }
  }
  // End caps: strips between the lower and upper chains, as triangles.
  for (const double z : {0.0, 1.0}) {
    for (int i = 0; i < 2 * nx; ++i) {
      const Index l0 = b.point({static_cast<double>(i) / nx, y_of(i, -1.0), z});
      const Index l1 = b.point({static_cast<double>(i + 1) / nx, y_of(i + 1, -1.0), z});
      const Index u0 = b.point({static_cast<double>(i) / nx, y_of(i, 1.0), z});
      const Index u1 = b.point({static_cast<double>(i + 1) / nx, y_of(i + 1, 1.0), z});
      std::vector<std::vector<Index>> cells;
      if (i == 0) {
        cells = {{l0, l1, u1}};
      } else if (i == 2 * nx - 1) {
        cells = {{l0, l1, u0}};
      } else {
        cells = {{l0, l1, u1}, {l0, u1, u0}};
      }
      for (auto cell : cells) {
        // Counter-clockwise seen from +z points out at z = 1 and in at z = 0.
        if (z > 0.5) {
          std::reverse(cell.begin(), cell.end());
        }
        add(cell);
      }
    }
  }
}

// Signed volume (divergence theorem) of polygons, quads split from their centroid.
double signed_volume(const std::vector<P3>& p, const std::vector<std::vector<Index>>& polys) {
  double v = 0.0;
  for (const auto& poly : polys) {
    P3 c{};
    for (const Index i : poly) {
      for (std::size_t k = 0; k < 3; ++k) {
        c[k] += p[u(i)][k] / static_cast<double>(poly.size());
      }
    }
    for (std::size_t i = 0; i < poly.size(); ++i) {
      v += det({0, 0, 0}, c, p[u(poly[i])], p[u(poly[(i + 1) % poly.size()])]) / 6.0;
    }
  }
  return v;
}

// Total volume of the cells (each from its outward faces).
double cell_volume(const og::Inflation3DResult& r) {
  double v = 0.0;
  for (const auto& t : r.tets) {
    v += signed_volume(
        r.points, {{t[0], t[2], t[1]}, {t[0], t[1], t[3]}, {t[1], t[2], t[3]}, {t[0], t[3], t[2]}});
  }
  for (const auto& y : r.pyramids) {
    v += signed_volume(r.points, {{y[0], y[3], y[2], y[1]},
                                  {y[0], y[1], y[4]},
                                  {y[1], y[2], y[4]},
                                  {y[2], y[3], y[4]},
                                  {y[3], y[0], y[4]}});
  }
  for (const auto& w : r.prisms) {
    v += signed_volume(r.points, {{w[0], w[2], w[1]},
                                  {w[3], w[4], w[5]},
                                  {w[0], w[1], w[4], w[3]},
                                  {w[1], w[2], w[5], w[4]},
                                  {w[2], w[0], w[3], w[5]}});
  }
  for (const auto& h : r.hexes) {
    v += signed_volume(r.points, {{h[0], h[3], h[2], h[1]},
                                  {h[4], h[5], h[6], h[7]},
                                  {h[0], h[1], h[5], h[4]},
                                  {h[1], h[2], h[6], h[5]},
                                  {h[2], h[3], h[7], h[6]},
                                  {h[3], h[0], h[4], h[7]}});
  }
  return v;
}

// Every edge of the front is used once in each direction (a closed, consistently oriented
// surface), and the layer cells fill the space between the wall and the front exactly.
void expect_closed_layers(const og::Inflation3DInput& in, const og::Inflation3DResult& r) {
  std::map<std::pair<Index, Index>, int> uses;
  for (const auto& poly : r.front) {
    for (std::size_t i = 0; i < poly.size(); ++i) {
      uses[{poly[i], poly[(i + 1) % poly.size()]}] += 1;
    }
  }
  for (const auto& [e, k] : uses) {
    EXPECT_EQ(k, 1);
    EXPECT_EQ(uses.count({e.second, e.first}), 1U);
  }
  std::vector<std::vector<Index>> wall;
  for (std::size_t f = 0; f < in.faces.size(); ++f) {
    if (in.wall[f]) {
      wall.push_back(in.faces[f]);
    }
  }
  const double between = signed_volume(r.points, wall) - signed_volume(r.points, r.front);
  EXPECT_GT(between, 0.0);
  EXPECT_NEAR(cell_volume(r), between, 1e-12 * std::abs(between) + 1e-15);
}

TEST(Inflation3D, PrismLayersOnClosedBody) {
  Builder b;
  b.box({0, 0, 0}, {4, 4, 4}, 4, 1.0, true, false);
  b.box({1.5, 1.5, 1.5}, {2.5, 2.5, 2.5}, 12, -1.0, true, true);
  // Without fans first: one prism per wall face and layer.
  auto plain = options(4);
  plain.fan_threshold_deg = 180.0;
  const auto r = og::inflate_3d(b.in, plain);
  std::size_t walls = 0;
  for (const bool w : b.in.wall) {
    walls += static_cast<std::size_t>(w);
  }
  EXPECT_EQ(r.prisms.size(), 4 * walls);
  EXPECT_TRUE(r.tets.empty());
  EXPECT_TRUE(r.pyramids.empty());
  EXPECT_EQ(r.max_layers, 4);
  EXPECT_EQ(r.front.size(), walls);
  expect_positive(r);
  expect_floor(r);
  expect_closed_layers(b.in, r);
  // With fans at the 90-degree edges (two fan columns each) and caps at the corners: every wall
  // node gets all layers (also with the transition to the core size, which must not stop a
  // corner's centre column for want of neighbours on other nodes).
  auto with_transition = options(4);
  with_transition.transition = true;
  const auto fans = og::inflate_3d(b.in, with_transition);
  EXPECT_TRUE(fans.stops.empty());
  expect_positive(fans);
  expect_floor(fans);
  expect_closed_layers(b.in, fans);
  EXPECT_EQ(fans.fan_columns, 12 * 11 * 2 + 8 * (3 * 2 + 1));
  for (std::size_t i = 0; i < b.in.points.size(); ++i) {
    if (!fans.reason[i].empty()) {
      EXPECT_EQ(fans.layers[i], 4);
    }
  }
  // One layer: every face-interior node rises exactly first_height normal to its face (the
  // marching normal is scaled by 1 / (n . n_face)). Four layers: the central node (normal
  // untouched by the edges) rise first_height * (1.2^4 - 1) / 0.2 straight up.
  const auto one = og::inflate_3d(b.in, options(1));
  int checked = 0;
  int central = 0;
  for (std::size_t i = 0; i < b.in.points.size(); ++i) {
    const P3& p = b.in.points[i];
    if (!(p[0] > 1.5 && p[0] < 2.5 && p[1] > 1.5 && p[1] < 2.5 && p[2] == 1.5)) {
      continue;
    }
    EXPECT_EQ(r.layers[i], 4);
    EXPECT_NEAR(one.points[u(one.top[i])][2], 1.5 - 0.01, 1e-14);
    ++checked;
    if (std::abs(p[0] - 2.0) < 0.01 && std::abs(p[1] - 2.0) < 0.01) {
      const P3& top = r.points[u(r.top[i])];
      EXPECT_NEAR(top[2], 1.5 - 0.01 * (std::pow(1.2, 4) - 1.0) / 0.2, 1e-12);
      EXPECT_NEAR(dist({p[0], p[1], 0.0}, {top[0], top[1], 0.0}), 0.0, 1e-6);
      ++central;
    }
  }
  EXPECT_EQ(checked, 121);
  EXPECT_EQ(central, 1);
}

TEST(Inflation3D, HexLayersOnQuads) {
  Builder b;
  b.box({0, 0, 0}, {4, 4, 4}, 4, 1.0, false, false);
  b.box({1.5, 1.5, 1.5}, {2.5, 2.5, 2.5}, 5, -1.0, false, true);
  auto o = options(3);
  o.fan_threshold_deg = 180.0;  // the box's edges without fans: one hexahedron per face and layer
  const auto r = og::inflate_3d(b.in, o);
  EXPECT_EQ(r.hexes.size(), 3U * 6U * 25U);
  EXPECT_TRUE(r.prisms.empty());
  expect_positive(r);
  expect_floor(r);
  for (const auto& f : r.front) {
    EXPECT_EQ(f.size(), 4U);
  }
}

TEST(Inflation3D, FixedNodesTaperAndTransitionCells) {
  // Only the top face of the body is a wall: its boundary nodes touch non-wall faces of other
  // planes, so they are fixed and the layers taper towards them.
  Builder b;
  b.box({0, 0, 0}, {4, 4, 4}, 4, 1.0, false, false);
  const auto faces = b.box({1.5, 1.5, 1.5}, {2.5, 2.5, 2.5}, 8, -1.0, false, false);
  for (const std::size_t f : faces[5]) {  // axis z, side max
    b.in.wall[f] = true;
  }
  // The taper mechanics without the quality floor (its transition cells are thin: one layer of
  // height 0.01 over cells 0.125 wide).
  auto o = options(6);
  o.floor = no_floor();
  const auto r = og::inflate_3d(b.in, o);
  EXPECT_GT(r.stops.at("fixed"), 0);
  EXPECT_GT(r.stops.at("slope"), 0);
  EXPECT_FALSE(r.hexes.empty());
  EXPECT_FALSE(r.pyramids.empty() && r.tets.empty());
  expect_positive(r);
  expect_floor(r, no_floor());
  for (std::size_t i = 0; i < b.in.points.size(); ++i) {
    const P3& p = b.in.points[i];
    if (p[2] == 2.5 && (p[0] == 1.5 || p[0] == 2.5 || p[1] == 1.5 || p[1] == 2.5)) {
      EXPECT_EQ(r.layers[i], 0);
      EXPECT_EQ(r.reason[i], "fixed");
    }
  }
  // With the quality floor (the default) no cell misses it: a column whose cells would miss it
  // stops instead. Here every step towards the fixed rim does (one layer of 0.01 over 0.125),
  // so the columns stop and the face keeps no layers.
  const auto strict = og::inflate_3d(b.in, options(6));
  expect_positive(strict);
  expect_floor(strict);
  EXPECT_GT(strict.stops.at("quality"), 0);
}

TEST(Inflation3D, SlidingNodesStayInTheirPlane) {
  // The bottom face of the body is not a wall but a sliding plane (z = 1.5) for its nodes.
  Builder b;
  b.box({0, 0, 0}, {4, 4, 4}, 4, 1.0, true, false);
  const auto faces = b.box({1.5, 1.5, 1.5}, {2.5, 2.5, 2.5}, 4, -1.0, true, true);
  for (const std::size_t f : faces[4]) {  // axis z, side min
    b.in.wall[f] = false;
    b.in.face_group[f] = 99;
    for (const Index p : b.in.faces[f]) {
      b.in.node_group[u(p)] = 99;
    }
  }
  const auto r = og::inflate_3d(b.in, options(3));
  EXPECT_FALSE(r.sides.empty());
  for (std::size_t i = 0; i < r.sides.size(); ++i) {
    EXPECT_EQ(b.in.face_group[u(r.side_face[i])], 99);
    for (const Index p : r.sides[i]) {
      EXPECT_DOUBLE_EQ(r.points[u(p)][2], 1.5);
    }
  }
  for (std::size_t i = 0; i < b.in.points.size(); ++i) {
    if (b.in.node_group[i] == 99) {
      EXPECT_DOUBLE_EQ(r.points[u(r.top[i])][2], 1.5);
      if (b.in.points[i][0] > 1.5 && b.in.points[i][0] < 2.5 && b.in.points[i][1] > 1.5 &&
          b.in.points[i][1] < 2.5) {
        continue;
      }
      EXPECT_EQ(r.layers[i], 3);
    }
  }
  expect_positive(r);
  expect_floor(r);
}

TEST(Inflation3D, ProximityStopsFacingFronts) {
  // Two bodies 0.1 apart: the facing columns stop before the gap closes.
  Builder b;
  b.box({0, 0, 0}, {5, 3, 3}, 5, 1.0, true, false);
  b.box({1, 1, 1}, {2, 2, 2}, 5, -1.0, true, true);
  b.box({2.1, 1, 1}, {3.1, 2, 2}, 5, -1.0, true, true);
  auto o = options(10);
  const auto r = og::inflate_3d(b.in, o);
  EXPECT_GT(r.stops.at("proximity"), 0);
  expect_positive(r);
  expect_floor(r);
  for (std::size_t i = 0; i < b.in.points.size(); ++i) {
    const P3& p = b.in.points[i];
    if ((p[0] == 2.0 || p[0] == 2.1) && p[1] > 1 && p[1] < 2 && p[2] > 1 && p[2] < 2) {
      EXPECT_LT(dist(p, r.points[u(r.top[i])]), 0.05);
    }
  }
}

TEST(Inflation3D, InvalidInput) {
  Builder b;
  b.box({0, 0, 0}, {1, 1, 1}, 2, 1.0, true, false);
  EXPECT_THROW(og::inflate_3d(b.in, options(3)), std::invalid_argument);  // no wall
  b.in.wall[0] = true;
  auto o = options(3);
  o.first_height = -1.0;
  EXPECT_THROW(og::inflate_3d(b.in, o), std::invalid_argument);
  o = options(3);
  b.in.faces[0] = {0, 1};
  EXPECT_THROW(og::inflate_3d(b.in, o), std::invalid_argument);
}

TEST(Inflation3D, FansAtSharpRidges) {
  Builder b;
  b.box({-2, -2, -2}, {4, 2, 3}, 3, 1.0, true, false);
  diamond(b, 0.2, 8, 4);
  auto o = options(5);
  o.first_height = 0.005;
  // Default: fans at every edge turning more than 60 degrees: the sharp edges (turning 180 -
  // 2 atan(0.2) = 157.4 degrees: five fan columns in 30-degree steps) and the 90-degree rims of
  // the end caps (two each); where they meet, a cap around a centre column closes the corner.
  // Every wall node gets all layers.
  const auto r = og::inflate_3d(b.in, o);
  expect_positive(r);
  expect_floor(r);
  expect_closed_layers(b.in, r);
  EXPECT_GT(r.fan_columns, 5 * 2 * 5);
  EXPECT_TRUE(r.stops.empty());
  for (std::size_t i = 0; i < b.in.points.size(); ++i) {
    if (!r.reason[i].empty()) {
      EXPECT_EQ(r.layers[i], 5) << "node " << i << " stopped: " << r.reason[i];
    }
  }
  // The first layer of a fan is made of prisms lying on their sides; hexahedra above.
  EXPECT_GT(r.prisms.size(), 0U);
  // Fans only at the sharp edges: the four ridge ends (where they meet the end caps) do not
  // grow and the fans narrow to them; the interior ridge nodes still get all layers (without
  // the quality floor: the steps towards the ridge ends are thin).
  o.fan_threshold_deg = 120.0;
  o.floor = no_floor();
  const auto ends = og::inflate_3d(b.in, o);
  expect_positive(ends);
  expect_floor(ends, no_floor());
  expect_closed_layers(b.in, ends);
  EXPECT_EQ(ends.fan_columns, 5 * 2 * 3);
  EXPECT_EQ(ends.stops.at("ridge end"), 4);
  int ridge_nodes = 0;
  for (std::size_t i = 0; i < b.in.points.size(); ++i) {
    const P3& p = b.in.points[i];
    if ((p[0] == 0.0 || p[0] == 2.0) && p[2] > 0.0 && p[2] < 1.0) {
      ++ridge_nodes;
      EXPECT_EQ(ends.layers[i], 5) << "ridge node " << i << " stopped: " << ends.reason[i];
    }
  }
  EXPECT_EQ(ridge_nodes, 6);
  o.floor = {};
  const auto strict_ends = og::inflate_3d(b.in, o);
  expect_positive(strict_ends);
  expect_floor(strict_ends);
  // Without fans the sharp edges stop at once (their normals are not visible from both sides).
  o.fan_threshold_deg = 180.0;
  const auto none = og::inflate_3d(b.in, o);
  EXPECT_EQ(none.fan_columns, 0);
  for (std::size_t i = 0; i < b.in.points.size(); ++i) {
    const P3& p = b.in.points[i];
    if ((p[0] == 0.0 || p[0] == 2.0) && p[2] > 0.0 && p[2] < 1.0) {
      EXPECT_LT(none.layers[i], 5);
    }
  }
}

}  // namespace
