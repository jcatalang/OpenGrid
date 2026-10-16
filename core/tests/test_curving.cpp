// Extension of boundary displacements into the mesh (highorder/curving.hpp).
#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "opengrid/core/highorder/curving.hpp"

namespace {

std::size_t at(int i, int j, int ni) {
  return static_cast<std::size_t>(j) * static_cast<std::size_t>(ni) + static_cast<std::size_t>(i);
}

// (ni x nj) rectilinear grid with row heights `heights` (j = 0 at y = 0), i fastest; edges
// along both directions.
og::DisplacementExtensionInput grid(int ni, const std::vector<double>& heights) {
  og::DisplacementExtensionInput in;
  const int nj = static_cast<int>(heights.size()) + 1;
  double y = 0.0;
  for (int j = 0; j < nj; ++j) {
    for (int i = 0; i < ni; ++i) {
      in.points.emplace_back(static_cast<double>(i) / (ni - 1), y, 0.0);
    }
    if (j + 1 < nj) {
      y += heights[static_cast<std::size_t>(j)];
    }
  }
  auto id = [ni](int i, int j) { return static_cast<std::int64_t>(j) * ni + i; };
  for (int j = 0; j < nj; ++j) {
    for (int i = 0; i < ni; ++i) {
      if (i + 1 < ni) {
        in.edges.push_back({id(i, j), id(i + 1, j)});
      }
      if (j + 1 < nj) {
        in.edges.push_back({id(i, j), id(i, j + 1)});
      }
    }
  }
  const std::size_t n = in.points.size();
  in.fixed.assign(n, 0);
  in.displacement.assign(n, og::Vec3{});
  in.radius.assign(n, 0.0);
  return in;
}

}  // namespace

TEST(Curving, UniformChainIsLinear) {
  og::DisplacementExtensionInput in;
  for (int k = 0; k <= 10; ++k) {
    in.points.emplace_back(0.1 * k, 0.0, 0.0);
    if (k > 0) {
      in.edges.push_back({k - 1, k});
      in.edges.push_back({k, k - 1});  // duplicates are merged
    }
  }
  in.fixed.assign(11, 0);
  in.fixed[0] = 1;
  in.fixed[10] = 1;
  in.displacement.assign(11, og::Vec3{});
  in.displacement[0] = {0.0, 1.0, 0.0};
  in.radius.assign(11, 0.0);
  in.radius[0] = 10.0;
  const og::DisplacementExtension out = og::extend_displacement(in);
  EXPECT_EQ(out.seeds, 1);
  EXPECT_EQ(out.solved, 9);
  for (int k = 0; k <= 10; ++k) {
    EXPECT_NEAR(out.displacement[static_cast<std::size_t>(k)].y, 1.0 - 0.1 * k, 1e-12);
  }
}

TEST(Curving, ThinLayersMoveWithTheWall) {
  // 20 layers of height 1e-5 under a far row: the wall row is pushed up by a bump (zero at its
  // ends, at most 1e-3 = 100 layer heights); every layer must keep a positive height, and the
  // layer heights change by a small fraction only.
  const int ni = 21;
  std::vector<double> heights(20, 1e-5);
  for (int k = 0; k < 8; ++k) {
    heights.push_back(heights.back() * 2.5);
  }
  og::DisplacementExtensionInput in = grid(ni, heights);
  const int nj = static_cast<int>(heights.size()) + 1;
  for (int i = 0; i < ni; ++i) {
    const auto wall = static_cast<std::size_t>(i);
    const double x = in.points[wall].x;
    in.fixed[wall] = 1;
    in.displacement[wall] = {0.0, 1e-3 * 4.0 * x * (1.0 - x), 0.0};
    in.radius[wall] = 0.5;
    const std::size_t top = at(i, nj - 1, ni);
    in.fixed[top] = 1;
  }
  for (int j = 1; j < nj - 1; ++j) {  // the sides are boundaries too
    in.fixed[at(0, j, ni)] = 1;
    in.fixed[at(ni - 1, j, ni)] = 1;
  }
  const og::DisplacementExtension out = og::extend_displacement(in);
  double worst = 0.0;
  for (int j = 0; j < 20; ++j) {
    for (int i = 1; i < ni - 1; ++i) {
      const std::size_t a = at(i, j, ni);
      const std::size_t b = at(i, j + 1, ni);
      const double before = in.points[b].y - in.points[a].y;
      const double after = before + out.displacement[b].y - out.displacement[a].y;
      ASSERT_GT(after, 0.0) << "layer " << j << " column " << i;
      worst = std::max(worst, std::abs(after / before - 1.0));
    }
  }
  EXPECT_LT(worst, 0.05);
}

TEST(Curving, RegionAndValidation) {
  const int ni = 41;
  og::DisplacementExtensionInput in = grid(ni, std::vector<double>(40, 0.025));
  const auto seed = static_cast<std::size_t>(20);  // middle of the bottom row
  in.fixed[seed] = 1;
  in.displacement[seed] = {0.0, 0.01, 0.0};
  in.radius[seed] = 0.1;  // four edges
  const og::DisplacementExtension out = og::extend_displacement(in);
  EXPECT_EQ(out.seeds, 1);
  EXPECT_GT(out.solved, 0);
  EXPECT_LT(out.solved, 80);  // a diamond of graph radius 4 around the seed
  for (std::size_t k = 0; k < in.points.size(); ++k) {
    if (std::abs(in.points[k].x - 0.5) + in.points[k].y > 0.1 + 1e-12) {
      EXPECT_EQ(out.displacement[k].y, 0.0);
    }
  }
  og::DisplacementExtensionInput bad = in;
  bad.edges.push_back({3, 3});
  EXPECT_THROW(og::extend_displacement(bad), std::invalid_argument);
  bad = in;
  bad.radius.pop_back();
  EXPECT_THROW(og::extend_displacement(bad), std::invalid_argument);
}
