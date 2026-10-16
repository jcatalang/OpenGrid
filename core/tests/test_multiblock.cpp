// Phases 10a / 10b: multi-block elliptic smoothing (structured/multiblock_elliptic*.hpp).
#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <numbers>
#include <utility>
#include <vector>

#include "opengrid/core/geom/arc.hpp"
#include "opengrid/core/geom/curve_chain.hpp"
#include "opengrid/core/structured/elliptic.hpp"
#include "opengrid/core/structured/elliptic3d.hpp"
#include "opengrid/core/structured/multiblock_elliptic.hpp"
#include "opengrid/core/structured/multiblock_elliptic3d.hpp"
#include "opengrid/core/structured/tfi.hpp"
#include "opengrid/core/util/parallel.hpp"

namespace {

std::size_t ix(int value) {
  return static_cast<std::size_t>(value);
}

// A channel with a bump on the bottom, TFI.
og::Grid2D bump_grid(int ni, int nj) {
  std::vector<og::Vec3> bottom;
  std::vector<og::Vec3> top;
  std::vector<og::Vec3> left;
  std::vector<og::Vec3> right;
  for (int i = 0; i < ni; ++i) {
    const double s = static_cast<double>(i) / (ni - 1);
    bottom.emplace_back(3.0 * s, 0.25 * std::pow(std::sin(std::numbers::pi * s), 4), 0.0);
    top.emplace_back(3.0 * s, 1.0, 0.0);
  }
  for (int j = 0; j < nj; ++j) {
    const double t = std::pow(static_cast<double>(j) / (nj - 1), 1.5);
    left.emplace_back(0.0, t, 0.0);
    right.emplace_back(3.0, t, 0.0);
  }
  return og::tfi_linear(og::Boundary2D{bottom, right, top, left});
}

std::array<std::int64_t, 4> quad(const std::vector<std::int64_t>& ids, int ni, int i, int j) {
  auto at = [&](int a, int b) { return ids[(ix(b) * ix(ni)) + ix(a)]; };
  return {at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)};
}

// The grid split at column k into two blocks sharing that column (merged nodes).
struct Split {
  og::MultiblockInput input;
  std::vector<std::int64_t> ids_a;
  std::vector<std::int64_t> ids_b;
};

Split split_grid(const og::Grid2D& g, int k) {
  Split s;
  const int ni = g.ni;
  const int nj = g.nj;
  const int na = k + 1;
  const int nb = ni - k;
  std::int64_t next = 0;
  s.ids_a.resize(ix(na) * ix(nj));
  s.ids_b.resize(ix(nb) * ix(nj));
  for (int j = 0; j < nj; ++j) {
    for (int i = 0; i < na; ++i) {
      s.ids_a[(ix(j) * ix(na)) + ix(i)] = next++;
      s.input.xy.push_back(g.at(i, j).x);
      s.input.xy.push_back(g.at(i, j).y);
    }
  }
  for (int j = 0; j < nj; ++j) {
    s.ids_b[ix(j) * ix(nb)] = s.ids_a[(ix(j) * ix(na)) + ix(k)];
    for (int i = 1; i < nb; ++i) {
      s.ids_b[(ix(j) * ix(nb)) + ix(i)] = next++;
      s.input.xy.push_back(g.at(k + i, j).x);
      s.input.xy.push_back(g.at(k + i, j).y);
    }
  }
  s.input.fixed.assign(static_cast<std::size_t>(next), 0);
  for (const auto& [ids, n, first, last] :
       {std::tuple{&s.ids_a, na, true, false}, std::tuple{&s.ids_b, nb, false, true}}) {
    for (int j = 0; j < nj - 1; ++j) {
      for (int i = 0; i < n - 1; ++i) {
        s.input.quads.push_back(quad(*ids, n, i, j));
      }
    }
    for (int i = 0; i < n; ++i) {
      s.input.fixed[static_cast<std::size_t>((*ids)[ix(i)])] = 1;
      s.input.fixed[static_cast<std::size_t>((*ids)[(ix(nj - 1) * ix(n)) + ix(i)])] = 1;
    }
    for (int j = 0; j < nj; ++j) {
      if (first) {
        s.input.fixed[static_cast<std::size_t>((*ids)[ix(j) * ix(n)])] = 1;
      }
      if (last) {
        s.input.fixed[static_cast<std::size_t>((*ids)[(ix(j) * ix(n)) + ix(n - 1)])] = 1;
      }
    }
  }
  return s;
}

TEST(Multiblock, SplitGridMatchesSingleBlockWithWall) {
  const int ni = 41;
  const int nj = 21;
  const og::Grid2D g = bump_grid(ni, nj);
  og::EllipticOptions single;
  single.control = og::ControlType::Laplace;
  single.bottom.enabled = true;
  single.bottom.spacing = 0.01;
  single.max_iterations = 20000;
  single.tolerance = 1e-12;
  single.acceleration = 5;
  const og::EllipticResult ref = og::smooth_elliptic(g, single);
  ASSERT_TRUE(ref.converged);

  Split s = split_grid(g, 17);
  for (const auto* ids : {&s.ids_a, &s.ids_b}) {
    const int n = static_cast<int>(ids->size()) / nj;
    for (int i = 0; i < n; ++i) {
      s.input.walls.push_back({(*ids)[static_cast<std::size_t>(i)], 0.01, 0.1, 0.1, 100.0});
    }
  }
  og::MultiblockOptions options;
  options.tolerance = 1e-12;
  options.max_iterations = 20000;
  const og::MultiblockResult r = og::smooth_multiblock(s.input, options);
  ASSERT_TRUE(r.converged);
  EXPECT_EQ(r.colours, 4);  // the zebra ADI ordering of one block
  double worst = 0.0;
  for (int j = 0; j < nj; ++j) {
    for (int i = 0; i < ni; ++i) {
      std::int64_t id = 0;
      if (i <= 17) {
        id = s.ids_a[(ix(j) * 18U) + ix(i)];
      } else {
        id = s.ids_b[(ix(j) * ix(ni - 17)) + ix(i - 17)];
      }
      const auto u = static_cast<std::size_t>(id);
      worst = std::max(worst, std::hypot(r.xy[2 * u] - ref.grid.at(i, j).x,
                                         r.xy[2 * u + 1] - ref.grid.at(i, j).y));
    }
  }
  EXPECT_LT(worst, 1e-9);
}

// Unit disc as an O-grid: a core block and four shells; corners of the core are irregular
// (three quads around them). Outer nodes slide on the circle.
TEST(Multiblock, OGridIrregularCornersAndSliding) {
  const int n = 13;  // points along each side
  const int nr = 9;  // radial points of the shells
  const double c = 0.4;
  // Node ids: core n x n, then each shell (n x nr) shares its inner row with the core side.
  og::MultiblockInput input;
  std::vector<std::int64_t> core(ix(n) * ix(n));
  auto add = [&](double x, double y) {
    input.xy.push_back(x);
    input.xy.push_back(y);
    return static_cast<std::int64_t>(input.xy.size() / 2 - 1);
  };
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      const double s = -c + 2.0 * c * i / (n - 1);
      const double t = -c + 2.0 * c * j / (n - 1);
      core[(ix(j) * ix(n)) + ix(i)] = add(s, t);
    }
  }
  for (int j = 0; j < n - 1; ++j) {
    for (int i = 0; i < n - 1; ++i) {
      input.quads.push_back(quad(core, n, i, j));
    }
  }
  // Core boundary in counter-clockwise order (4 (n - 1) nodes), starting at (-c, -c).
  std::vector<std::int64_t> ring;
  ring.reserve(4 * ix(n - 1));
  for (int i = 0; i < n - 1; ++i) {
    ring.push_back(core[static_cast<std::size_t>(i)]);
  }
  for (int j = 0; j < n - 1; ++j) {
    ring.push_back(core[(ix(j) * ix(n)) + ix(n - 1)]);
  }
  for (int i = n - 1; i > 0; --i) {
    ring.push_back(core[(ix(n - 1) * ix(n)) + ix(i)]);
  }
  for (int j = n - 1; j > 0; --j) {
    ring.push_back(core[ix(j) * ix(n)]);
  }
  const auto m = ring.size();
  // Shell rows r = 1..nr-1 between the core ring and the circle (one ring of nodes per row).
  std::vector<std::vector<std::int64_t>> rows{ring};
  for (int r = 1; r < nr; ++r) {
    std::vector<std::int64_t> row;
    const double f = static_cast<double>(r) / (nr - 1);
    for (std::size_t k = 0; k < m; ++k) {
      const auto u = static_cast<std::size_t>(ring[k]);
      const double x0 = input.xy[2 * u];
      const double y0 = input.xy[2 * u + 1];
      const double a = std::atan2(y0, x0);
      row.push_back(add((1 - f) * x0 + f * std::cos(a), (1 - f) * y0 + f * std::sin(a)));
    }
    rows.push_back(row);
  }
  for (int r = 0; r + 1 < nr; ++r) {
    for (std::size_t k = 0; k < m; ++k) {
      const std::size_t k1 = (k + 1) % m;
      // Counter-clockwise: the ring runs counter-clockwise and row r + 1 lies outside row r.
      const auto& inner = rows[static_cast<std::size_t>(r)];
      const auto& outer = rows[ix(r) + 1];
      input.quads.push_back({inner[k], outer[k], outer[k1], inner[k1]});
    }
  }
  const std::size_t count = input.xy.size() / 2;
  input.fixed.assign(count, 0);
  // Outer circle: the four corner directions fixed, the rest sliding on four arcs.
  for (int side = 0; side < 4; ++side) {
    const double a0 = -0.75 * std::numbers::pi + side * 0.5 * std::numbers::pi;
    auto arc = std::make_shared<og::Arc>(og::Vec3{0, 0, 0}, 1.0, og::Vec3{0, 0, 1},
                                         og::Vec3{1, 0, 0}, a0, a0 + 0.5 * std::numbers::pi);
    input.chains.push_back(
        std::make_shared<og::CurveChain>(std::vector<og::ChainSegment>{{arc, 0.0, 1.0}}, 1e-12));
    for (int k = 0; k < n - 1; ++k) {
      const std::int64_t node = rows.back()[(ix(side) * ix(n - 1)) + ix(k)];
      if (k == 0) {
        input.fixed[static_cast<std::size_t>(node)] = 1;
      } else {
        input.slides.push_back({node, static_cast<std::size_t>(side)});
      }
    }
  }
  const og::MultiblockResult r = og::smooth_multiblock(input, og::MultiblockOptions{});
  ASSERT_TRUE(r.converged);
  int irregular = 0;
  int sliding = 0;
  for (const std::uint8_t kind : r.node_class) {
    irregular += static_cast<int>(kind == static_cast<std::uint8_t>(og::MultiblockNode::Irregular));
    sliding += static_cast<int>(kind == static_cast<std::uint8_t>(og::MultiblockNode::Sliding));
  }
  EXPECT_EQ(irregular, 4);
  EXPECT_EQ(sliding, 4 * (n - 2));
  for (const og::MultiblockSlide& s : input.slides) {
    const auto u = static_cast<std::size_t>(s.node);
    EXPECT_NEAR(std::hypot(r.xy[2 * u], r.xy[2 * u + 1]), 1.0, 1e-12);
  }
  // Symmetric layout: the irregular corners stay on the diagonals (to the solver tolerance;
  // the sweep order is not symmetric).
  for (std::size_t g = 0; g < count; ++g) {
    if (r.node_class[g] == static_cast<std::uint8_t>(og::MultiblockNode::Irregular)) {
      EXPECT_NEAR(std::fabs(r.xy[2 * g]), std::fabs(r.xy[2 * g + 1]), 1e-5);
    }
  }
  // Same result for any thread count.
  const int before = og::thread_count();
  og::set_thread_count(1);
  const og::MultiblockResult one = og::smooth_multiblock(input, og::MultiblockOptions{});
  og::set_thread_count(before);
  EXPECT_EQ(one.xy, r.xy);
  EXPECT_EQ(one.iterations, r.iterations);
}

TEST(Multiblock, RejectsNodesOutOfRange) {
  og::MultiblockInput input;
  input.xy = {0, 0, 1, 0, 1, 1};
  input.fixed = {1, 1, 1};
  input.quads.push_back({0, 1, 2, 7});
  EXPECT_THROW((void)og::smooth_multiblock(input, og::MultiblockOptions{}), std::invalid_argument);
}

}  // namespace

// 3D: a curved block split into two blocks sharing a plane of nodes matches the single-block
// smoother (Laplace, faces fixed), for any thread count.
TEST(Multiblock3D, SplitBlockMatchesSingleBlock) {
  const int ni = 13;
  const int nj = 11;
  const int nk = 9;
  og::Grid3D g(ni, nj, nk);
  for (int k = 0; k < nk; ++k) {
    for (int j = 0; j < nj; ++j) {
      for (int i = 0; i < ni; ++i) {
        const double s = static_cast<double>(i) / (ni - 1);
        const double t = static_cast<double>(j) / (nj - 1);
        const double r = std::pow(static_cast<double>(k) / (nk - 1), 1.3);
        double x = 2.0 * s;
        double y = t + 0.2 * std::sin(std::numbers::pi * s) * (1.0 - r);
        double z = r + 0.1 * std::sin(std::numbers::pi * t) * s;
        if (i > 0 && i < ni - 1 && j > 0 && j < nj - 1 && k > 0 && k < nk - 1) {
          x += 0.01 * std::sin(7.0 * i + 3.0 * j);
          y += 0.01 * std::cos(5.0 * j + 2.0 * k);
          z += 0.01 * std::sin(3.0 * k + i);
        }
        g.at(i, j, k) = {x, y, z};
      }
    }
  }
  og::Elliptic3DOptions single;
  single.control = og::ControlType::Laplace;
  single.max_iterations = 20000;
  single.tolerance = 1e-12;
  single.acceleration = 5;
  const og::Elliptic3DResult ref = og::smooth_elliptic_3d(g, single);
  ASSERT_TRUE(ref.converged);

  og::Multiblock3DInput input;
  auto id = [&](int i, int j, int k) {
    const std::size_t index = (((ix(k) * ix(nj)) + ix(j)) * ix(ni)) + ix(i);
    return static_cast<std::int64_t>(index);
  };
  for (const og::Vec3& p : g.x) {
    input.xyz.insert(input.xyz.end(), {p.x, p.y, p.z});
  }
  input.fixed.assign(g.x.size(), 0);
  for (int k = 0; k < nk; ++k) {
    for (int j = 0; j < nj; ++j) {
      for (int i = 0; i < ni; ++i) {
        if (i == 0 || j == 0 || k == 0 || i == ni - 1 || j == nj - 1 || k == nk - 1) {
          input.fixed[static_cast<std::size_t>(id(i, j, k))] = 1;
        }
      }
    }
  }
  // Two blocks: i in [0, 5] and [5, ni - 1] (the plane i = 5 shared).
  for (const auto& [i0, i1] : {std::pair{0, 5}, std::pair{5, ni - 1}}) {
    for (int k = 0; k < nk - 1; ++k) {
      for (int j = 0; j < nj - 1; ++j) {
        for (int i = i0; i < i1; ++i) {
          std::array<std::int64_t, 8> h{};
          for (unsigned c = 0; c < 8; ++c) {
            h[c] = id(i + static_cast<int>(c & 1U), j + static_cast<int>((c >> 1U) & 1U),
                      k + static_cast<int>((c >> 2U) & 1U));
          }
          input.hexes.push_back(h);
        }
      }
    }
  }
  og::MultiblockOptions options;
  options.tolerance = 1e-12;
  options.max_iterations = 20000;
  const og::Multiblock3DResult r = og::smooth_multiblock_3d(input, options);
  ASSERT_TRUE(r.converged);
  double worst = 0.0;
  for (std::size_t n = 0; n < g.x.size(); ++n) {
    worst = std::max(worst,
                     std::hypot(r.xyz[3 * n] - ref.grid.x[n].x, r.xyz[3 * n + 1] - ref.grid.x[n].y,
                                r.xyz[3 * n + 2] - ref.grid.x[n].z));
  }
  EXPECT_LT(worst, 1e-9);
  const int before = og::thread_count();
  og::set_thread_count(1);
  const og::Multiblock3DResult one = og::smooth_multiblock_3d(input, options);
  og::set_thread_count(before);
  EXPECT_EQ(one.xyz, r.xyz);
}

TEST(Multiblock3D, RejectsBadSlides) {
  og::Multiblock3DInput input;
  input.xyz.assign(24, 0.0);
  input.fixed.assign(8, 1);
  input.hexes.push_back({0, 1, 2, 3, 4, 5, 6, 7});
  input.slides.push_back({0, og::SlideTarget::Surface, 3});
  EXPECT_THROW((void)og::smooth_multiblock_3d(input, og::MultiblockOptions{}),
               std::invalid_argument);
}

// Phase 10d: internal walls (a fixed row inside the grid, controlled on both sides), the
// absolute tolerance of warm-started runs and the slide limit.
namespace {

// A unit square, (ni x nj) uniform, the middle row j = nj / 2 fixed (an internal wall).
og::MultiblockInput square_with_middle_row(int ni, int nj, double wall_spacing) {
  og::MultiblockInput input;
  std::vector<std::int64_t> ids(ix(ni) * ix(nj));
  for (int j = 0; j < nj; ++j) {
    for (int i = 0; i < ni; ++i) {
      const std::size_t index = (ix(j) * ix(ni)) + ix(i);
      ids[index] = static_cast<std::int64_t>(index);
      input.xy.push_back(static_cast<double>(i) / (ni - 1));
      input.xy.push_back((2.0 * j / (nj - 1)) - 1.0);
      const bool boundary = i == 0 || j == 0 || i == ni - 1 || j == nj - 1 || j == nj / 2;
      input.fixed.push_back(static_cast<std::uint8_t>(boundary));
    }
  }
  for (int j = 0; j + 1 < nj; ++j) {
    for (int i = 0; i + 1 < ni; ++i) {
      input.quads.push_back(quad(ids, ni, i, j));
    }
  }
  for (int i = 0; i < ni; ++i) {
    input.walls.push_back({ids[(ix(nj / 2) * ix(ni)) + ix(i)], wall_spacing, 0.1, 0.1, 100.0});
  }
  return input;
}

TEST(Multiblock, InternalWallControlsBothSides) {
  const int ni = 21;
  const int nj = 41;
  const double h = 0.01;  // the uniform grid has 0.05
  const og::MultiblockInput input = square_with_middle_row(ni, nj, h);
  og::MultiblockOptions options;
  options.tolerance = 1e-10;
  options.max_iterations = 20000;
  const og::MultiblockResult r = og::smooth_multiblock(input, options);
  ASSERT_TRUE(r.converged);
  auto y = [&](int i, int j) { return r.xy[(2 * ((ix(j) * ix(ni)) + ix(i))) + 1]; };
  for (int i = 4; i < ni - 4; ++i) {
    EXPECT_NEAR(y(i, (nj / 2) + 1) - y(i, nj / 2), h, 0.05 * h) << i;
    EXPECT_NEAR(y(i, nj / 2) - y(i, (nj / 2) - 1), h, 0.05 * h) << i;
    EXPECT_EQ(y(i, nj / 2), 0.0);
  }
}

TEST(Multiblock, AbsoluteToleranceStopsAWarmStart) {
  const int ni = 21;
  const int nj = 41;
  og::MultiblockInput input = square_with_middle_row(ni, nj, 0.01);
  input.walls.clear();  // wall sources restart from zero: positions only here
  og::MultiblockOptions options;
  options.tolerance = 1e-9;
  options.max_iterations = 20000;
  const og::MultiblockResult cold = og::smooth_multiblock(input, options);
  ASSERT_TRUE(cold.converged);
  // From the converged grid, relative to its own (tiny) first update the run would go on;
  // the absolute tolerance (the cold run's final accuracy) stops it at once.
  input.xy = cold.xy;
  options.absolute_tolerance = options.tolerance * cold.residuals.front();
  const og::MultiblockResult warm = og::smooth_multiblock(input, options);
  EXPECT_TRUE(warm.converged);
  EXPECT_LE(warm.iterations, 3);
}

TEST(Multiblock3D, InternalWallControlsBothSides) {
  const int n = 13;
  const int nk = 21;
  const double h = 0.02;  // the uniform grid has 0.1
  og::Multiblock3DInput input;
  auto id = [&](int i, int j, int k) {
    const std::size_t index = (((ix(k) * ix(n)) + ix(j)) * ix(n)) + ix(i);
    return static_cast<std::int64_t>(index);
  };
  for (int k = 0; k < nk; ++k) {
    for (int j = 0; j < n; ++j) {
      for (int i = 0; i < n; ++i) {
        input.xyz.insert(input.xyz.end(),
                         {static_cast<double>(i) / (n - 1), static_cast<double>(j) / (n - 1),
                          (2.0 * k / (nk - 1)) - 1.0});
        const bool fixed =
            i == 0 || j == 0 || k == 0 || i == n - 1 || j == n - 1 || k == nk - 1 || k == nk / 2;
        input.fixed.push_back(static_cast<std::uint8_t>(fixed));
        if (k == nk / 2) {
          input.walls.push_back({id(i, j, k), h, 0.1, 0.1, 100.0});
        }
      }
    }
  }
  for (int k = 0; k + 1 < nk; ++k) {
    for (int j = 0; j + 1 < n; ++j) {
      for (int i = 0; i + 1 < n; ++i) {
        std::array<std::int64_t, 8> hex{};
        for (unsigned c = 0; c < 8; ++c) {
          hex[c] = id(i + static_cast<int>(c & 1U), j + static_cast<int>((c >> 1U) & 1U),
                      k + static_cast<int>((c >> 2U) & 1U));
        }
        input.hexes.push_back(hex);
      }
    }
  }
  og::MultiblockOptions options;
  options.tolerance = 1e-10;
  options.max_iterations = 20000;
  const og::Multiblock3DResult r = og::smooth_multiblock_3d(input, options);
  ASSERT_TRUE(r.converged);
  auto z = [&](int i, int j, int k) {
    return r.xyz[(3 * static_cast<std::size_t>(id(i, j, k))) + 2];
  };
  for (int i = 3; i < n - 3; ++i) {
    EXPECT_NEAR(z(i, n / 2, (nk / 2) + 1) - z(i, n / 2, nk / 2), h, 0.05 * h) << i;
    EXPECT_NEAR(z(i, n / 2, nk / 2) - z(i, n / 2, (nk / 2) - 1), h, 0.05 * h) << i;
  }
}

}  // namespace
