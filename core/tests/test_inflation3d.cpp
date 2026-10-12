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

TEST(Inflation3D, PrismLayersOnClosedBody) {
  Builder b;
  b.box({0, 0, 0}, {4, 4, 4}, 4, 1.0, true, false);
  b.box({1.5, 1.5, 1.5}, {2.5, 2.5, 2.5}, 12, -1.0, true, true);
  const auto r = og::inflate_3d(b.in, options(4));
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
  const auto r = og::inflate_3d(b.in, options(3));
  EXPECT_EQ(r.hexes.size(), 3U * 6U * 25U);
  EXPECT_TRUE(r.prisms.empty());
  expect_positive(r);
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
  const auto r = og::inflate_3d(b.in, options(6));
  EXPECT_GT(r.stops.at("fixed"), 0);
  EXPECT_GT(r.stops.at("slope"), 0);
  EXPECT_FALSE(r.hexes.empty());
  EXPECT_FALSE(r.pyramids.empty() && r.tets.empty());
  expect_positive(r);
  for (std::size_t i = 0; i < b.in.points.size(); ++i) {
    const P3& p = b.in.points[i];
    if (p[2] == 2.5 && (p[0] == 1.5 || p[0] == 2.5 || p[1] == 1.5 || p[1] == 2.5)) {
      EXPECT_EQ(r.layers[i], 0);
      EXPECT_EQ(r.reason[i], "fixed");
    }
  }
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

}  // namespace
