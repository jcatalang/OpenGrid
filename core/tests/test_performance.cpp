// Phase 9c: parallel loops, Anderson acceleration, volume quality metrics.
#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "opengrid/core/quality/volume_quality.hpp"
#include "opengrid/core/util/anderson.hpp"
#include "opengrid/core/util/parallel.hpp"

namespace {

TEST(Parallel, ForCoversTheRangeOnce) {
  std::vector<int> hits(100000, 0);
  og::parallel_for(0, static_cast<std::int64_t>(hits.size()),
                   [&](std::int64_t i) { hits[static_cast<std::size_t>(i)] += 1; });
  for (const int h : hits) {
    ASSERT_EQ(h, 1);
  }
  EXPECT_GE(og::thread_count(), 1);
  EXPECT_THROW(og::set_thread_count(0), std::invalid_argument);
}

// Jacobi iteration for a 1D Poisson problem: slow plain, fast with Anderson acceleration, and
// identical results for any thread count.
std::vector<double> jacobi(int depth, int threads, int& iterations) {
  og::set_thread_count(threads);
  const std::size_t n = 200;
  std::vector<double> u(n, 0.0);
  u.back() = 1.0;  // Dirichlet ends 0 and 1
  og::AndersonAccelerator anderson(n, depth);
  std::vector<double> g(n);
  iterations = 0;
  for (int k = 0; k < 200000; ++k) {
    g = u;
    double change = 0.0;
    for (std::size_t i = 1; i + 1 < n; ++i) {
      g[i] = 0.5 * (u[i - 1] + u[i + 1]);
      change = std::max(change, std::fabs(g[i] - u[i]));
    }
    iterations = k + 1;
    if (change < 1e-13) {
      break;
    }
    if (depth > 0) {
      anderson.step(u, g);
    }
    u = g;
  }
  return u;
}

TEST(Anderson, AcceleratesAndIsThreadIndependent) {
  const int before = og::thread_count();
  int plain = 0;
  int fast = 0;
  int fast_one = 0;
  const std::vector<double> a = jacobi(0, before, plain);
  const std::vector<double> b = jacobi(10, before, fast);
  const std::vector<double> c = jacobi(10, 1, fast_one);
  og::set_thread_count(before);
  EXPECT_LT(fast * 5, plain) << fast << " vs " << plain;
  EXPECT_EQ(b, c);
  EXPECT_EQ(fast, fast_one);
  for (std::size_t i = 0; i < a.size(); ++i) {
    const double exact = static_cast<double>(i) / static_cast<double>(a.size() - 1);
    EXPECT_NEAR(b[i], exact, 1e-8);
  }
}

TEST(Anderson, GuardRestoresTheBestIterateAndSwitchesOff) {
  og::AndersonAccelerator anderson(2, 5);
  std::vector<double> u = {0.0, 0.0};
  std::vector<double> g = {1.0, 1.0};
  EXPECT_FALSE(anderson.guarded_step(u, g, 1.0));  // no history yet
  u = g;
  g = {1.5, 1.25};
  anderson.guarded_step(u, g, 0.5);  // the lowest residual: {1.5, 1.25} kept
  u = g;
  g = {50.0, -40.0};
  EXPECT_TRUE(anderson.guarded_step(u, g, 0.5 * og::AndersonAccelerator::kBlowup * 2.0));
  EXPECT_EQ(g, (std::vector<double>{1.5, 1.25}));
  EXPECT_EQ(anderson.restarts(), 1);
  EXPECT_TRUE(anderson.enabled());
  g = {70.0, 70.0};
  EXPECT_TRUE(anderson.guarded_step(u, g, 1e6));
  EXPECT_EQ(g, (std::vector<double>{1.5, 1.25}));
  EXPECT_FALSE(anderson.enabled());  // kMaxRestarts = 2: the plain iteration continues
  g = {2.0, 2.0};
  EXPECT_FALSE(anderson.guarded_step(u, g, 0.1));
  EXPECT_EQ(g, (std::vector<double>{2.0, 2.0}));
  EXPECT_NE(anderson.restart_note().find("without acceleration"), std::string::npos);
}

TEST(VolumeQuality, RegularCellsAndNeighbours) {
  // Two unit cubes side by side (hexahedra), and a regular tetrahedron apart.
  std::vector<double> p;
  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 3; ++i) {
        p.insert(p.end(), {static_cast<double>(i), static_cast<double>(j), static_cast<double>(k)});
      }
    }
  }
  auto id = [](int i, int j, int k) { return static_cast<std::int64_t>((k * 2 + j) * 3 + i); };
  std::vector<std::int64_t> hexes;
  for (int i = 0; i < 2; ++i) {
    for (const auto& [a, b, c] : std::array<std::array<int, 3>, 8>{{{i, 0, 0},
                                                                    {i + 1, 0, 0},
                                                                    {i + 1, 1, 0},
                                                                    {i, 1, 0},
                                                                    {i, 0, 1},
                                                                    {i + 1, 0, 1},
                                                                    {i + 1, 1, 1},
                                                                    {i, 1, 1}}}) {
      hexes.push_back(id(a, b, c));
    }
  }
  const auto base = static_cast<std::int64_t>(p.size() / 3);
  const double h = std::sqrt(2.0 / 3.0);
  p.insert(p.end(), {10.0, 0.0, 0.0, 11.0, 0.0, 0.0, 10.5, std::sqrt(3.0) / 2.0, 0.0, 10.5,
                     std::sqrt(3.0) / 6.0, h});
  const std::vector<std::int64_t> tets = {base, base + 1, base + 2, base + 3};
  const og::VolumeQuality q = og::volume_quality(
      p.data(), p.size() / 3,
      {{tets.data(), 1, 4}, {nullptr, 0, 5}, {nullptr, 0, 6}, {hexes.data(), 2, 8}});
  ASSERT_EQ(q.volume.size(), 3U);
  EXPECT_NEAR(q.scaled_jacobian[0], 1.0, 1e-12);
  EXPECT_NEAR(q.aspect_ratio[0], 1.0, 1e-12);
  EXPECT_NEAR(q.skewness[0], 0.0, 1e-9);
  EXPECT_NEAR(q.volume[0], 1.0 / (6.0 * std::sqrt(2.0)), 1e-12);  // regular, edge 1
  for (int c = 1; c < 3; ++c) {
    EXPECT_NEAR(q.scaled_jacobian[static_cast<std::size_t>(c)], 1.0, 1e-12);
    EXPECT_NEAR(q.volume[static_cast<std::size_t>(c)], 1.0, 1e-12);
    EXPECT_NEAR(q.min_angle[static_cast<std::size_t>(c)], 90.0, 1e-9);
    EXPECT_NEAR(q.non_orthogonality[static_cast<std::size_t>(c)], 0.0, 1e-6);
    EXPECT_NEAR(q.volume_ratio[static_cast<std::size_t>(c)], 1.0, 1e-12);
  }
}

}  // namespace
