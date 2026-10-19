#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <numbers>
#include <random>
#include <vector>

#include "opengrid/core/geom/arc.hpp"
#include "opengrid/core/geom/surface.hpp"
#include "opengrid/core/quality/quality3d.hpp"
#include "opengrid/core/structured/elliptic3d.hpp"
#include "opengrid/core/structured/hyperbolic.hpp"
#include "opengrid/core/structured/surface_elliptic.hpp"
#include "opengrid/core/structured/tfi3d.hpp"

namespace {

using og::Vec3;
constexpr double kPi = std::numbers::pi;

Vec3 cyl(double r, double t, double z) {
  return {r * std::cos(t), r * std::sin(t), z};
}

og::Grid3D sector(int ni, int nj, int nk) {
  og::Grid3D g(ni, nj, nk);
  for (int k = 0; k < nk; ++k) {
    for (int j = 0; j < nj; ++j) {
      for (int i = 0; i < ni; ++i) {
        g.at(i, j, k) = cyl(1.0 + static_cast<double>(i) / (ni - 1), (kPi / 3) * j / (nj - 1),
                            static_cast<double>(k) / (nk - 1));
      }
    }
  }
  return g;
}

void perturb(og::Grid3D& g, double amplitude, unsigned seed) {
  std::mt19937 rng(seed);  // NOLINT(bugprone-random-generator-seed): reproducible test
  std::uniform_real_distribution<double> u(-amplitude, amplitude);
  for (int k = 1; k < g.nk - 1; ++k) {
    for (int j = 1; j < g.nj - 1; ++j) {
      for (int i = 1; i < g.ni - 1; ++i) {
        g.at(i, j, k) += Vec3{u(rng), u(rng), u(rng)};
      }
    }
  }
}

double min_jacobian(const og::Grid3D& g) {
  const og::Quality3D q = og::quality_3d(g);
  double m = 1.0;
  for (const double v : q.scaled_jacobian) {
    m = std::min(m, v);
  }
  return m;
}

}  // namespace

TEST(Elliptic3D, LaplaceKeepsUniformCube) {
  og::Grid3D g(6, 5, 7);
  for (int k = 0; k < 7; ++k) {
    for (int j = 0; j < 5; ++j) {
      for (int i = 0; i < 6; ++i) {
        g.at(i, j, k) = {i * 0.2, j * 0.25, k / 6.0};
      }
    }
  }
  const og::Elliptic3DResult r = og::smooth_elliptic_3d(g, {});
  EXPECT_TRUE(r.converged);
  for (std::size_t n = 0; n < g.x.size(); ++n) {
    EXPECT_LT(og::norm(r.grid.x[n] - g.x[n]), 1e-12);
  }
}

TEST(Elliptic3D, ThomasMiddlecoffConvergesToTheSameGridFromAnyStart) {
  // The smoothed grid is the unique solution for the given faces: starting from the exact sector
  // or from a perturbed interior gives the same result.
  const og::Grid3D exact = sector(11, 15, 9);
  og::Grid3D noisy = exact;
  perturb(noisy, 0.03, 5);
  og::Elliptic3DOptions options;
  options.control = og::ControlType::ThomasMiddlecoff;
  options.max_iterations = 5000;
  options.tolerance = 1e-12;
  const og::Elliptic3DResult a = og::smooth_elliptic_3d(exact, options);
  const og::Elliptic3DResult b = og::smooth_elliptic_3d(noisy, options);
  EXPECT_TRUE(a.converged) << a.iterations;
  EXPECT_TRUE(b.converged) << b.iterations;
  double err = 0.0;
  for (std::size_t n = 0; n < exact.x.size(); ++n) {
    err = std::max(err, og::norm(a.grid.x[n] - b.grid.x[n]));
  }
  EXPECT_LT(err, 1e-8);
  EXPECT_GT(min_jacobian(b.grid), 0.95);
}

TEST(Elliptic3D, WallControlSpacingAndOrthogonality) {
  // Block over a bump: x in [0, 2], y in [0, 1], z from the bump to 1; uniform TFI start.
  const int ni = 21;
  const int nj = 9;
  const int nk = 17;
  auto bump = [](double x) { return 0.1 * std::sin(kPi * x / 2.0); };
  og::Grid3D shape(ni, nj, nk);
  for (int k = 0; k < nk; ++k) {
    for (int j = 0; j < nj; ++j) {
      for (int i = 0; i < ni; ++i) {
        const double x = 2.0 * i / (ni - 1);
        const double z0 = bump(x);
        shape.at(i, j, k) = {x, static_cast<double>(j) / (nj - 1), z0 + (1.0 - z0) * k / (nk - 1)};
      }
    }
  }
  const og::Grid3D initial = og::tfi_block(og::block_faces(shape));
  og::Elliptic3DOptions options;
  og::WallControl3D wall;
  wall.face = 4;
  wall.control.enabled = true;
  wall.control.spacing = 0.02;
  wall.control.decay = 0.4;
  wall.control.relaxation = 0.1;
  options.walls.push_back(wall);
  options.max_iterations = 4000;
  options.tolerance = 1e-9;
  const og::Elliptic3DResult r = og::smooth_elliptic_3d(initial, options);
  EXPECT_GT(min_jacobian(r.grid), 0.0);
  double worst_spacing = 0.0;
  double worst_angle = 0.0;
  for (int j = 1; j < nj - 1; ++j) {
    for (int i = 1; i < ni - 1; ++i) {
      const Vec3 line = r.grid.at(i, j, 1) - r.grid.at(i, j, 0);
      const Vec3 ti = r.grid.at(i + 1, j, 0) - r.grid.at(i - 1, j, 0);
      const Vec3 tj = r.grid.at(i, j + 1, 0) - r.grid.at(i, j - 1, 0);
      const Vec3 n = og::cross(ti, tj);
      worst_spacing = std::max(worst_spacing, std::fabs(og::norm(line) / 0.02 - 1.0));
      const double angle =
          std::acos(std::clamp(og::dot(line, n) / (og::norm(line) * og::norm(n)), -1.0, 1.0));
      worst_angle = std::max(worst_angle, angle * 180.0 / kPi);
    }
  }
  EXPECT_LT(worst_spacing, 0.02) << r.iterations;
  EXPECT_LT(worst_angle, 1.0);
  // Faces are untouched.
  EXPECT_EQ(og::norm(r.grid.at(5, 0, 7) - initial.at(5, 0, 7)), 0.0);
}

TEST(SurfaceElliptic, StaysOnCylinderAndSmooths) {
  auto arc =
      std::make_shared<og::Arc>(Vec3{0, 0, 0}, 2.0, Vec3{0, 0, 1}, Vec3{1, 0, 0}, 0.0, kPi / 2);
  const og::ExtrudedSurface cylinder(arc, {0, 0, 3});
  const int ni = 17;
  const int nj = 13;
  og::Grid2D exact(ni, nj);
  for (int j = 0; j < nj; ++j) {
    for (int i = 0; i < ni; ++i) {
      exact.at(i, j) = cyl(2.0, 0.5 * kPi * i / (ni - 1), 3.0 * j / (nj - 1));
    }
  }
  og::Grid2D noisy = exact;
  std::mt19937 rng(9);  // NOLINT(bugprone-random-generator-seed): reproducible test
  std::uniform_real_distribution<double> u(-0.05, 0.05);
  for (int j = 1; j < nj - 1; ++j) {
    for (int i = 1; i < ni - 1; ++i) {
      const double t = 0.5 * kPi * i / (ni - 1) + u(rng);
      noisy.at(i, j) = cyl(2.0, t, 3.0 * j / (nj - 1) + u(rng));
    }
  }
  og::EllipticOptions options;
  options.control = og::ControlType::ThomasMiddlecoff;
  options.max_iterations = 3000;
  options.tolerance = 1e-10;
  const og::EllipticResult r = og::smooth_on_surface(cylinder, noisy, options);
  EXPECT_TRUE(r.converged);
  double err = 0.0;
  for (std::size_t n = 0; n < exact.x.size(); ++n) {
    EXPECT_NEAR(std::hypot(r.grid.x[n].x, r.grid.x[n].y), 2.0, 1e-11);
    err = std::max(err, og::norm(r.grid.x[n] - exact.x[n]));
  }
  EXPECT_LT(err, 2e-3);
  og::EllipticOptions bad;
  bad.bottom.enabled = true;
  EXPECT_THROW(og::smooth_on_surface(cylinder, noisy, bad), std::invalid_argument);
}

TEST(Hyperbolic, FlatPlateHeightsAreExact) {
  og::Grid2D base(9, 7);
  for (int j = 0; j < 7; ++j) {
    for (int i = 0; i < 9; ++i) {
      base.at(i, j) = {i * 0.25, j * 0.2, 0.0};
    }
  }
  std::vector<double> s;
  double total = 0.0;
  double h = 1.0;
  std::vector<double> heights;
  for (int k = 0; k < 10; ++k) {
    heights.push_back(h);
    total += h;
    h *= 1.2;
  }
  s.push_back(0.0);
  double acc = 0.0;
  for (const double v : heights) {
    acc += v;
    s.push_back(acc / total);
  }
  s.back() = 1.0;
  const og::HyperbolicResult r = og::extrude_hyperbolic(base, 0.5, s, {});
  EXPECT_EQ(r.inverted_layers, 0);
  for (int k = 0; k < 11; ++k) {
    for (int j = 0; j < 7; ++j) {
      for (int i = 0; i < 9; ++i) {
        EXPECT_NEAR(r.grid.at(i, j, k).z, 0.5 * s[static_cast<std::size_t>(k)], 1e-12);
        EXPECT_NEAR(r.grid.at(i, j, k).x, base.at(i, j).x, 1e-12);
      }
    }
  }
}

TEST(Hyperbolic, CylinderOGridIsRadial) {
  // Cylinder wall, theta around (periodic i, counter-clockwise seen from +z so that x_i x x_j
  // points outward), z along j with constant-z ends.
  const int ni = 49;
  const int nj = 9;
  og::Grid2D wall(ni, nj);
  for (int j = 0; j < nj; ++j) {
    for (int i = 0; i < ni; ++i) {
      double t = 2.0 * kPi * i / (ni - 1);
      if (i == ni - 1) {
        t = 0.0;
      }
      wall.at(i, j) = cyl(0.5, t, 2.0 * j / (nj - 1));
    }
  }
  std::vector<double> s;
  for (int k = 0; k <= 30; ++k) {
    s.push_back((std::pow(1.15, k) - 1.0) / (std::pow(1.15, 30) - 1.0));
  }
  og::HyperbolicOptions options;
  options.edges = {og::HyperbolicEdge::Periodic, og::HyperbolicEdge::Periodic,
                   og::HyperbolicEdge::ConstantZ, og::HyperbolicEdge::ConstantZ};
  const og::HyperbolicResult r = og::extrude_hyperbolic(wall, 9.5, s, options);
  EXPECT_FALSE(r.mirrored);
  EXPECT_EQ(r.inverted_layers, 0);
  for (int k = 0; k <= 30; ++k) {
    const double radius = 0.5 + 9.5 * s[static_cast<std::size_t>(k)];
    for (int j = 0; j < nj; ++j) {
      for (int i = 0; i < ni; ++i) {
        const Vec3& p = r.grid.at(i, j, k);
        EXPECT_NEAR(std::hypot(p.x, p.y) / radius, 1.0, 1e-9);
        EXPECT_NEAR(p.z, wall.at(i, j).z, 1e-12);
      }
    }
  }
  EXPECT_GT(min_jacobian(r.grid), 0.0);
  // Marching inwards with a negative distance mirrors i.
  const og::HyperbolicResult in = og::extrude_hyperbolic(wall, -0.2, {0.0, 0.5, 1.0}, options);
  EXPECT_TRUE(in.mirrored);
  EXPECT_GT(min_jacobian(in.grid), 0.0);
  options.edges[0] = og::HyperbolicEdge::Periodic;
  og::Grid2D open = wall;
  open.at(ni - 1, 3) += Vec3{0.1, 0, 0};
  EXPECT_THROW(og::extrude_hyperbolic(open, 1.0, s, options), std::invalid_argument);
}
