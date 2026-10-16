// Wall layers (structured/wall_layers.hpp): marching along the front normals.
#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <stdexcept>
#include <vector>

#include "opengrid/core/structured/wall_layers.hpp"

namespace {

constexpr std::int64_t kNone = -1;

// An ni x nj grid of an annulus (wall r = 1 at j = 0, r = 2 at the top; closed around, or a
// quarter of it) whose interior nodes are sheared along the arc, so its grid lines leave the
// wall obliquely.
struct Annulus {
  int ni = 0;
  int nj = 0;
  bool closed = true;
  std::vector<og::Vec3> x;
  og::WallLayerInput input;
  [[nodiscard]] std::int64_t id(int i, int j) const {
    return j * ni + (i % ni);
  }
};

Annulus annulus(int ni, int nj, int layers, double spacing, double growth, bool closed = true) {
  Annulus a;
  a.ni = ni;
  a.nj = nj;
  a.closed = closed;
  double sweep = 0.5 * std::numbers::pi;
  int cells_around = ni - 1;
  if (closed) {
    sweep = 2.0 * std::numbers::pi;
    cells_around = ni;
  }
  for (int j = 0; j < nj; ++j) {
    for (int i = 0; i < ni; ++i) {
      const double s = static_cast<double>(j) / (nj - 1);
      const double t = static_cast<double>(i) / cells_around;
      const double shear = 0.3 * s * (1.0 - s) * std::sin(2.0 * std::numbers::pi * t);
      const double theta = sweep * t + shear;
      const double r = 1.0 + s;
      a.x.emplace_back(r * std::cos(theta), r * std::sin(theta), 0.0);
    }
  }
  const auto n = static_cast<std::size_t>(ni) * static_cast<std::size_t>(nj);
  a.input.dim = 2;
  a.input.movable.assign(n, 0);
  a.input.on_wall.assign(n, 0);
  for (int j = 1; j < nj - 1; ++j) {
    for (int i = 0; i < ni; ++i) {
      if (closed || (i > 0 && i < ni - 1)) {
        a.input.movable[static_cast<std::size_t>(a.id(i, j))] = 1;
      }
    }
  }
  for (int i = 0; i < ni; ++i) {
    a.input.on_wall[static_cast<std::size_t>(a.id(i, 0))] = 1;
    a.input.walls.push_back({a.id(i, 0), spacing, layers, growth, static_cast<std::size_t>(i)});
  }
  for (int j = 0; j + 1 < nj; ++j) {
    for (int i = 0; i < cells_around; ++i) {
      // Counter-clockwise in the xy plane: out along r first (theta grows with i).
      a.input.cells.push_back({a.id(i, j), a.id(i, j + 1), a.id(i + 1, j + 1), a.id(i + 1, j),
                               kNone, kNone, kNone, kNone});
    }
  }
  return a;
}

TEST(WallLayers, LinesLeaveACurvedWallNormallyWithGeometricHeights) {
  Annulus a = annulus(40, 31, 12, 1e-4, 1.25);
  const og::WallLayerResult result = og::march_wall_layers(a.input, a.x);
  ASSERT_EQ(result.lines.size(), 40U);
  for (const og::WallLayerLine& line : result.lines) {
    ASSERT_EQ(line.nodes.size(), 12U);
    EXPECT_NEAR(line.next_spacing, 1e-4 * std::pow(1.25, 12), 1e-15);
    const og::Vec3 base = a.x[static_cast<std::size_t>(line.base)];
    const double theta = std::atan2(base.y, base.x);
    double height = 0.0;
    for (std::size_t k = 0; k < line.nodes.size(); ++k) {
      height += 1e-4 * std::pow(1.25, static_cast<double>(k));
      const og::Vec3 p = a.x[static_cast<std::size_t>(line.nodes[k])];
      // On the radial line through the wall node, at the summed height: normal to the wall.
      EXPECT_NEAR(std::atan2(p.y, p.x), theta, 1e-9);
      EXPECT_NEAR(std::hypot(p.x, p.y) - 1.0, height, 1e-12);
    }
  }
  // The grid beyond moved along its lines and stays valid: the next node is still outside.
  for (int i = 0; i < 40; ++i) {
    const og::Vec3 last = a.x[static_cast<std::size_t>(a.id(i, 12))];
    const og::Vec3 next = a.x[static_cast<std::size_t>(a.id(i, 13))];
    EXPECT_GT(std::hypot(next.x, next.y), std::hypot(last.x, last.y));
  }
}

TEST(WallLayers, FrontEndsUseOneSidedNormals) {
  // A quarter annulus: where the wall meets the side boundaries the last lines see one facet,
  // so they lean by at most about half a cell's angle (2.25 degrees here).
  Annulus a = annulus(21, 31, 12, 1e-4, 1.25, false);
  const og::WallLayerResult result = og::march_wall_layers(a.input, a.x);
  ASSERT_EQ(result.lines.size(), 19U);  // the wall's inner nodes
  for (const og::WallLayerLine& line : result.lines) {
    const og::Vec3 base = a.x[static_cast<std::size_t>(line.base)];
    const og::Vec3 top = a.x[static_cast<std::size_t>(line.nodes.back())];
    const double lean = std::fabs(std::atan2(top.y, top.x) - std::atan2(base.y, base.x));
    EXPECT_LT(lean, 2.5 * std::numbers::pi / 180.0);
  }
}

TEST(WallLayersReport, TooFewNodesAndOverlaps) {
  Annulus a = annulus(12, 9, 8, 1e-3, 1.2);  // only 7 free nodes on each line
  EXPECT_THROW(og::march_wall_layers(a.input, a.x), std::invalid_argument);
  Annulus b = annulus(12, 9, 4, 1e-3, 1.2);
  // The top is a wall too, with layers that meet the bottom's.
  for (int i = 0; i < 12; ++i) {
    b.input.on_wall[static_cast<std::size_t>(b.id(i, 8))] = 1;
    b.input.walls.push_back({b.id(i, 8), 1e-3, 4, 1.2, 100});
  }
  EXPECT_THROW(og::march_wall_layers(b.input, b.x), std::invalid_argument);
  Annulus c = annulus(12, 9, 3, 1e-3, 0.5);
  EXPECT_THROW(og::march_wall_layers(c.input, c.x), std::invalid_argument);  // growth < 1
}

TEST(WallLayers, HexahedraOffAPlaneWallWithAShearedGrid) {
  // n^3 nodes over the unit cube, interior nodes sheared in x; the wall is z = 0.
  const int m = 9;
  auto id = [&](int i, int j, int k) { return (static_cast<std::int64_t>(k) * m + j) * m + i; };
  std::vector<og::Vec3> x;
  for (int k = 0; k < m; ++k) {
    for (int j = 0; j < m; ++j) {
      for (int i = 0; i < m; ++i) {
        const double z = static_cast<double>(k) / (m - 1);
        double shear = 0.0;
        if (i > 0 && i < m - 1 && j > 0 && j < m - 1 && k > 0 && k < m - 1) {
          shear = 0.05 * z;
        }
        x.emplace_back(static_cast<double>(i) / (m - 1) + shear, static_cast<double>(j) / (m - 1),
                       z);
      }
    }
  }
  og::WallLayerInput input;
  input.dim = 3;
  const auto n = x.size();
  input.movable.assign(n, 0);
  input.on_wall.assign(n, 0);
  for (int k = 1; k < m - 1; ++k) {
    for (int j = 1; j < m - 1; ++j) {
      for (int i = 1; i < m - 1; ++i) {
        input.movable[static_cast<std::size_t>(id(i, j, k))] = 1;
      }
    }
  }
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < m; ++i) {
      input.on_wall[static_cast<std::size_t>(id(i, j, 0))] = 1;
      input.walls.push_back({id(i, j, 0), 0.01, 4, 1.3, 0});
    }
  }
  for (int k = 0; k + 1 < m; ++k) {
    for (int j = 0; j + 1 < m; ++j) {
      for (int i = 0; i + 1 < m; ++i) {
        std::array<std::int64_t, 8> h{};
        for (unsigned c = 0; c < 8; ++c) {
          h[c] = id(i + static_cast<int>(c & 1U), j + static_cast<int>((c >> 1U) & 1U),
                    k + static_cast<int>((c >> 2U) & 1U));
        }
        input.cells.push_back(h);
      }
    }
  }
  const std::vector<og::Vec3> initial = x;
  const og::WallLayerResult result = og::march_wall_layers(input, x);
  ASSERT_EQ(result.lines.size(), static_cast<std::size_t>((m - 2) * (m - 2)));
  for (const og::WallLayerLine& line : result.lines) {
    const og::Vec3 base = x[static_cast<std::size_t>(line.base)];
    double height = 0.0;
    for (std::size_t k = 0; k < line.nodes.size(); ++k) {
      height += 0.01 * std::pow(1.3, static_cast<double>(k));
      const og::Vec3 p = x[static_cast<std::size_t>(line.nodes[k])];
      EXPECT_NEAR(p.x, base.x, 1e-14);  // straight up from the wall, not sheared
      EXPECT_NEAR(p.y, base.y, 1e-14);
      EXPECT_NEAR(p.z, height, 1e-14);
    }
  }
  // Nodes beyond the layers moved along their lines, less towards the top.
  const og::Vec3 near = x[static_cast<std::size_t>(id(4, 4, 5))];
  const og::Vec3 far = x[static_cast<std::size_t>(id(4, 4, 7))];
  const double moved_near = std::fabs(near.x - initial[static_cast<std::size_t>(id(4, 4, 5))].x);
  const double moved_far = std::fabs(far.x - initial[static_cast<std::size_t>(id(4, 4, 7))].x);
  EXPECT_GT(moved_near, moved_far);
}

}  // namespace
