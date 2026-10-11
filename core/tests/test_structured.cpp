#include <gtest/gtest.h>

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <thread>

#include "opengrid/core/dist/distribution.hpp"
#include "opengrid/core/quality/quality2d.hpp"
#include "opengrid/core/structured/elliptic.hpp"
#include "opengrid/core/structured/tfi.hpp"

using og::Vec3;

namespace {

std::vector<double> stretched(int n, double begin, double end) {
  og::DistributionSpec spec;
  spec.type = og::DistributionType::Tanh;
  spec.begin_spacing = begin;
  spec.end_spacing = end;
  return og::distribute(n, spec).s;
}

// Quarter annulus r in [r0, r1], theta in [0, pi/2], counter-clockwise: i radial (bottom is the
// ray theta = 0), j angular (left is the inner arc).
og::Boundary2D annulus(const std::vector<double>& rs, const std::vector<double>& ts, double r0,
                       double r1) {
  og::Boundary2D b;
  const double quarter = std::numbers::pi / 2.0;
  for (const double s : rs) {
    const double r = r0 + s * (r1 - r0);
    b.bottom.emplace_back(r, 0.0, 0.0);
    b.top.emplace_back(r * std::cos(quarter), r, 0.0);
  }
  for (const double t : ts) {
    const double th = t * quarter;
    b.left.emplace_back(r0 * std::cos(th), r0 * std::sin(th), 0.0);
    b.right.emplace_back(r1 * std::cos(th), r1 * std::sin(th), 0.0);
  }
  b.top.front() = b.left.back();
  b.top.back() = b.right.back();
  return b;
}

og::Boundary2D rectangle(const std::vector<double>& xs, const std::vector<double>& ys, double w,
                         double h) {
  og::Boundary2D b;
  for (const double s : xs) {
    b.bottom.emplace_back(s * w, 0.0, 0.0);
    b.top.emplace_back(s * w, h, 0.0);
  }
  for (const double t : ys) {
    b.left.emplace_back(0.0, t * h, 0.0);
    b.right.emplace_back(w, t * h, 0.0);
  }
  return b;
}

// Bump channel: bottom y = 0.15 sin(pi x), top y = 1, sides clustered at the bottom.
og::Boundary2D bump(int ni, int nj, double wall) {
  og::Boundary2D b;
  for (int i = 0; i < ni; ++i) {
    const double x = static_cast<double>(i) / (ni - 1);
    b.bottom.emplace_back(x, 0.15 * std::sin(std::numbers::pi * x), 0.0);
    b.top.emplace_back(x, 1.0, 0.0);
  }
  og::DistributionSpec spec;
  spec.type = og::DistributionType::Tanh;
  spec.begin_spacing = wall;
  const auto ys = og::distribute(nj, spec).s;
  for (const double y : ys) {
    b.left.emplace_back(0.0, y, 0.0);
    b.right.emplace_back(1.0, y, 0.0);
  }
  return b;
}

double min_jacobian(const og::Grid2D& g) {
  const auto q = og::quality_2d(g);
  return *std::min_element(q.scaled_jacobian.begin(), q.scaled_jacobian.end());
}

}  // namespace

TEST(Tfi, ExactOnQuarterAnnulus) {
  const auto rs = stretched(31, 1e-3, 0.05);
  const auto ts = stretched(41, 0.01, 0.04);
  const auto b = annulus(rs, ts, 1.0, 3.0);
  const og::Grid2D g = og::tfi_linear(b);
  double worst = 0.0;
  for (int j = 0; j < g.nj; ++j) {
    for (int i = 0; i < g.ni; ++i) {
      const double r = 1.0 + 2.0 * rs[static_cast<std::size_t>(i)];
      const double th = ts[static_cast<std::size_t>(j)] * std::numbers::pi / 2.0;
      worst = std::max(worst, og::distance(g.at(i, j), {r * std::cos(th), r * std::sin(th), 0.0}));
    }
  }
  EXPECT_LT(worst, 1e-12 * 3.0);
  EXPECT_GT(min_jacobian(g), 0.0);
}

TEST(Tfi, ReproducesBoundariesAndValidates) {
  const auto b = bump(21, 11, 0.01);
  for (const og::Grid2D& g : {og::tfi_linear(b), og::tfi_hermite(b, {})}) {
    for (int i = 0; i < g.ni; ++i) {
      EXPECT_EQ(g.at(i, 0).y, b.bottom[static_cast<std::size_t>(i)].y);
      EXPECT_EQ(g.at(i, g.nj - 1).y, 1.0);
    }
    EXPECT_GT(min_jacobian(g), 0.0);
  }
  og::Boundary2D bad = b;
  bad.top.pop_back();
  EXPECT_THROW(og::tfi_linear(bad), std::invalid_argument);
  bad = b;
  bad.right.back().x += 0.1;
  EXPECT_THROW(og::tfi_linear(bad), std::invalid_argument);
}

TEST(Tfi, HermiteWithoutWallsEqualsLinear) {
  const auto b = bump(21, 11, 0.02);
  og::HermiteOptions none;
  none.bottom = false;
  none.top = false;
  const og::Grid2D h = og::tfi_hermite(b, none);
  const og::Grid2D l = og::tfi_linear(b);
  for (std::size_t k = 0; k < h.x.size(); ++k) {
    EXPECT_LT(og::distance(h.x[k], l.x[k]), 1e-13);
  }
}

TEST(Tfi, HermiteIsOrthogonalAtTheWall) {
  const auto b = bump(41, 21, 0.01);
  const og::Grid2D linear = og::tfi_linear(b);
  const og::Grid2D hermite = og::tfi_hermite(b, {});
  auto worst_deviation = [](const og::Grid2D& g) {
    double worst = 0.0;
    for (int i = 5; i < g.ni - 5; ++i) {
      const Vec3 t = g.at(i + 1, 0) - g.at(i - 1, 0);
      const Vec3 n = g.at(i, 1) - g.at(i, 0);
      const double angle =
          std::acos(og::dot(t, n) / (og::norm(t) * og::norm(n))) * 180.0 / std::numbers::pi;
      worst = std::max(worst, std::fabs(angle - 90.0));
    }
    return worst;
  };
  EXPECT_GT(worst_deviation(linear), 10.0);
  EXPECT_LT(worst_deviation(hermite), 1.0);
}

TEST(Elliptic, UniformRectangleIsAFixedPoint) {
  const auto u = og::distribute(11, {}).s;
  const og::Grid2D g = og::tfi_linear(rectangle(u, og::distribute(7, {}).s, 2.0, 1.0));
  og::EllipticOptions options;
  options.max_iterations = 50;
  const auto result = og::smooth_elliptic(g, options);
  for (std::size_t k = 0; k < g.x.size(); ++k) {
    EXPECT_LT(og::distance(result.grid.x[k], g.x[k]), 1e-14);
  }
}

TEST(Elliptic, ThomasMiddlecoffPreservesSeparableClustering) {
  const og::Grid2D g =
      og::tfi_linear(rectangle(stretched(41, 1e-3, 0.05), stretched(31, 1e-3, 0.08), 3.0, 1.0));
  og::EllipticOptions options;
  options.control = og::ControlType::ThomasMiddlecoff;
  options.max_iterations = 200;
  const auto tm = og::smooth_elliptic(g, options);
  double moved = 0.0;
  for (std::size_t k = 0; k < g.x.size(); ++k) {
    moved = std::max(moved, og::distance(tm.grid.x[k], g.x[k]));
  }
  EXPECT_LT(moved, 1e-9);
  // Plain Laplace smoothing destroys the clustering.
  options.control = og::ControlType::Laplace;
  const auto laplace = og::smooth_elliptic(g, options);
  EXPECT_GT(og::distance(laplace.grid.at(1, 1), g.at(1, 1)), 1e-4);
}

TEST(Elliptic, LaplaceConvergesOnAnnulus) {
  for (const auto solver : {og::EllipticSolver::PointSor, og::EllipticSolver::LineSor}) {
    const auto b = annulus(stretched(31, 0.01, 0.05), og::distribute(41, {}).s, 1.0, 3.0);
    og::EllipticOptions options;
    options.solver = solver;
    options.max_iterations = 5000;
    options.tolerance = 1e-8;
    const auto result = og::smooth_elliptic(og::tfi_linear(b), options);
    EXPECT_TRUE(result.converged);
    EXPECT_LE(result.residuals.back(), 1e-8 * result.residuals.front());
    EXPECT_GT(min_jacobian(result.grid), 0.0);
  }
}

TEST(Elliptic, StegerSorensonWallSpacingAndOrthogonality) {
  const double wall = 2e-3;
  const auto b = bump(61, 41, wall);
  og::EllipticOptions options;
  options.control = og::ControlType::ThomasMiddlecoff;
  options.bottom.enabled = true;
  options.bottom.spacing = wall;
  options.max_iterations = 3000;
  options.tolerance = 1e-7;
  const auto result = og::smooth_elliptic(og::tfi_linear(b), options);
  EXPECT_TRUE(result.warnings.empty()) << result.warnings.front();
  const og::Grid2D& g = result.grid;
  EXPECT_GT(min_jacobian(g), 0.0);
  for (int i = 6; i < g.ni - 6; ++i) {
    const Vec3 t = g.at(i + 1, 0) - g.at(i - 1, 0);
    const Vec3 n = g.at(i, 1) - g.at(i, 0);
    const double angle =
        std::acos(og::dot(t, n) / (og::norm(t) * og::norm(n))) * 180.0 / std::numbers::pi;
    EXPECT_NEAR(angle, 90.0, 2.0) << "i=" << i;
    EXPECT_NEAR(og::norm(n), wall, 0.05 * wall) << "i=" << i;
  }
}

TEST(Elliptic, CancelAndInvalidInput) {
  const og::Grid2D g = og::tfi_linear(bump(11, 11, 0.05));
  og::SmoothingProgress progress;
  progress.cancel();
  const auto result = og::smooth_elliptic(g, {}, &progress);
  EXPECT_TRUE(result.cancelled);
  EXPECT_EQ(result.iterations, 0);
  og::Grid2D tilted = g;
  tilted.at(3, 3).z = 0.5;
  EXPECT_THROW(og::smooth_elliptic(tilted, {}), std::invalid_argument);
}

TEST(Elliptic, ProgressFromAnotherThread) {
  const og::Grid2D g = og::tfi_linear(bump(81, 61, 1e-3));
  og::EllipticOptions options;
  options.max_iterations = 100000;
  options.tolerance = 1e-30;
  og::SmoothingProgress progress;
  std::thread worker([&] { og::smooth_elliptic(g, options, &progress); });
  while (progress.iteration() < 5) {
    std::this_thread::yield();
  }
  progress.cancel();
  worker.join();
  EXPECT_GE(progress.residuals().size(), 5U);
  EXPECT_LT(progress.residuals().size(), 100000U);
}

TEST(Quality, UnitSquares) {
  const og::Grid2D g =
      og::tfi_linear(rectangle(og::distribute(5, {}).s, og::distribute(5, {}).s, 4.0, 4.0));
  const auto q = og::quality_2d(g);
  ASSERT_EQ(q.scaled_jacobian.size(), 16U);
  for (std::size_t c = 0; c < 16; ++c) {
    EXPECT_NEAR(q.scaled_jacobian[c], 1.0, 1e-15);
    EXPECT_NEAR(q.skewness[c], 0.0, 1e-13);
    EXPECT_NEAR(q.aspect_ratio[c], 1.0, 1e-15);
    EXPECT_NEAR(q.min_angle[c], 90.0, 1e-12);
    EXPECT_NEAR(q.max_angle[c], 90.0, 1e-12);
    EXPECT_NEAR(q.non_orthogonality[c], 0.0, 1e-6);
    EXPECT_NEAR(q.area_ratio[c], 1.0, 1e-15);
    EXPECT_NEAR(q.area[c], 1.0, 1e-15);
  }
}

TEST(Quality, ShearedStretchedAndInverted) {
  og::Grid2D g(2, 2);
  const double s = std::tan(std::numbers::pi / 6.0);  // 30 degree shear
  g.at(0, 0) = {0, 0, 0};
  g.at(1, 0) = {2, 0, 0};
  g.at(0, 1) = {s, 1, 0};
  g.at(1, 1) = {2 + s, 1, 0};
  const auto q = og::quality_2d(g);
  EXPECT_NEAR(q.min_angle[0], 60.0, 1e-12);
  EXPECT_NEAR(q.max_angle[0], 120.0, 1e-12);
  EXPECT_NEAR(q.skewness[0], 1.0 / 3.0, 1e-12);
  EXPECT_NEAR(q.scaled_jacobian[0], std::sin(std::numbers::pi / 3.0), 1e-12);
  EXPECT_NEAR(q.aspect_ratio[0], 2.0 / std::hypot(s, 1.0), 1e-12);
  std::swap(g.at(0, 0), g.at(1, 0));
  std::swap(g.at(0, 1), g.at(1, 1));
  EXPECT_LT(og::quality_2d(g).scaled_jacobian[0], 0.0);
}

TEST(Quality, AreaRatioAndNonOrthogonality) {
  // Geometric growth 1.5 in x: neighbour area ratio 1.5; skewed column gives non-orthogonality.
  const og::Grid2D g =
      og::tfi_linear(rectangle(og::distribute(6,
                                              [] {
                                                og::DistributionSpec spec;
                                                spec.type = og::DistributionType::Geometric;
                                                spec.growth_ratio = 1.5;
                                                return spec;
                                              }())
                                   .s,
                               og::distribute(3, {}).s, 1.0, 1.0));
  const auto q = og::quality_2d(g);
  EXPECT_NEAR(q.area_ratio[1], 1.5, 1e-12);
  EXPECT_NEAR(q.non_orthogonality[1], 0.0, 1e-6);
  og::Grid2D skew = g;
  skew.at(2, 1).x += 0.05;
  EXPECT_GT(og::quality_2d(skew).non_orthogonality[1], 1.0);
}
