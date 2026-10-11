#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <numbers>
#include <random>
#include <vector>

#include "opengrid/core/geom/arc.hpp"
#include "opengrid/core/geom/line.hpp"
#include "opengrid/core/geom/surface.hpp"
#include "opengrid/core/quality/quality3d.hpp"
#include "opengrid/core/structured/extrude.hpp"
#include "opengrid/core/structured/surface_tfi.hpp"
#include "opengrid/core/structured/tfi3d.hpp"

namespace {

using og::Vec3;
constexpr double kPi = std::numbers::pi;

Vec3 cyl(double r, double t, double z) {
  return {r * std::cos(t), r * std::sin(t), z};
}

// Annular sector r in [r0, r1], theta in [0, t1], z in [0, h] with uniform nodes: i = r, j = theta,
// k = z.
og::Grid3D sector(int ni, int nj, int nk, double r0, double r1, double t1, double h) {
  og::Grid3D g(ni, nj, nk);
  for (int k = 0; k < nk; ++k) {
    for (int j = 0; j < nj; ++j) {
      for (int i = 0; i < ni; ++i) {
        g.at(i, j, k) = cyl(r0 + (r1 - r0) * i / (ni - 1), t1 * j / (nj - 1), h * k / (nk - 1));
      }
    }
  }
  return g;
}

double total(const std::vector<double>& v) {
  double s = 0.0;
  for (const double x : v) {
    s += x;
  }
  return s;
}

void expect_projection(const og::Surface& s, int seed) {
  std::mt19937 rng(static_cast<unsigned>(seed));  // NOLINT(bugprone-random-generator-seed)
  std::uniform_real_distribution<double> uniform(0.02, 0.98);
  const double size = s.size();
  for (int trial = 0; trial < 40; ++trial) {
    const double u = uniform(rng);
    const double v = uniform(rng);
    const Vec3 p = s.eval(u, v);
    const og::SurfaceProjection q = s.project(p);
    EXPECT_LT(q.distance, 1e-10 * size) << "u " << u << " v " << v;
  }
  // Analytic derivatives against central differences.
  for (int trial = 0; trial < 10; ++trial) {
    const double u = uniform(rng);
    const double v = uniform(rng);
    const double h = 1e-6;
    const Vec3 fu = (s.eval(u + h, v) - s.eval(u - h, v)) / (2 * h);
    const Vec3 fv = (s.eval(u, v + h) - s.eval(u, v - h)) / (2 * h);
    EXPECT_LT(og::norm(fu - s.du(u, v)), 1e-6 * size);
    EXPECT_LT(og::norm(fv - s.dv(u, v)), 1e-6 * size);
  }
}

}  // namespace

TEST(Surfaces, ProjectionAndDerivatives) {
  const og::PlaneSurface plane({1, 2, 3}, {2, 0, 0.5}, {0, 1.5, 0});
  expect_projection(plane, 1);
  // Off-surface point: exact distance.
  const Vec3 n = plane.normal(0.3, 0.4);
  const auto q = plane.project(plane.eval(0.3, 0.4) + 0.25 * n);
  EXPECT_NEAR(q.distance, 0.25, 1e-12);

  auto arc =
      std::make_shared<og::Arc>(Vec3{0, 0, 0}, 2.0, Vec3{0, 0, 1}, Vec3{1, 0, 0}, 0.0, kPi / 2);
  const og::ExtrudedSurface cylinder(arc, {0, 0, 3});
  expect_projection(cylinder, 2);
  const auto qc = cylinder.project(cyl(2.5, 0.7, 1.2));
  EXPECT_NEAR(qc.distance, 0.5, 1e-12);

  auto inner =
      std::make_shared<og::Arc>(Vec3{0, 0, 0}, 1.0, Vec3{0, 0, 1}, Vec3{1, 0, 0}, 0.0, kPi / 2);
  auto outer =
      std::make_shared<og::Arc>(Vec3{0, 0, 1}, 2.0, Vec3{0, 0, 1}, Vec3{1, 0, 0}, 0.0, kPi / 2);
  const og::RuledSurface ruled(inner, outer);
  expect_projection(ruled, 3);

  // Revolving a line parallel to the axis gives a cylinder.
  auto generator = std::make_shared<og::Line>(Vec3{1.5, 0, 0}, Vec3{1.5, 0, 2});
  const og::RevolvedSurface revolved(generator, {0, 0, 0}, {0, 0, 1}, 0.0, kPi);
  expect_projection(revolved, 4);
  const auto qr = revolved.project(cyl(1.0, 1.0, 0.5));
  EXPECT_NEAR(qr.distance, 0.5, 1e-12);

  // Quarter cylinder as a rational quadratic NURBS (exact circle in u, linear in v).
  const double w = std::sqrt(0.5);
  std::vector<Vec3> control;
  std::vector<double> weights;
  for (const auto& [p, wi] : {std::pair{Vec3{1, 0, 0}, 1.0}, std::pair{Vec3{1, 1, 0}, w},
                              std::pair{Vec3{0, 1, 0}, 1.0}}) {
    for (const double z : {0.0, 2.0}) {
      control.push_back(p + Vec3{0, 0, z});
      weights.push_back(wi);
    }
  }
  const og::NurbsSurface nurbs(2, 1, {0, 0, 0, 1, 1, 1}, {0, 0, 1, 1}, 3, 2, control, weights);
  expect_projection(nurbs, 5);
  for (int k = 0; k <= 10; ++k) {
    const Vec3 p = nurbs.eval(k / 10.0, 0.3);
    EXPECT_NEAR(std::hypot(p.x, p.y), 1.0, 1e-14);
    EXPECT_NEAR(p.z, 0.6, 1e-14);
  }
  const auto qn = nurbs.project(cyl(1.3, 0.4, 1.0));
  EXPECT_NEAR(qn.distance, 0.3, 1e-12);
}

TEST(BlockTfi, AnnularSectorIsExactAndVolumeConverges) {
  const og::Grid3D exact = sector(9, 13, 7, 1.0, 2.0, kPi / 3, 1.5);
  const og::Grid3D tfi = og::tfi_block(og::block_faces(exact));
  double err = 0.0;
  for (std::size_t n = 0; n < exact.x.size(); ++n) {
    err = std::max(err, og::norm(tfi.x[n] - exact.x[n]));
  }
  EXPECT_LT(err, 1e-13);
  // Fine in theta: volume within 1e-6 of the analytic sector volume.
  const og::Grid3D fine =
      og::tfi_block(og::block_faces(sector(3, 1201, 3, 1.0, 2.0, kPi / 3, 1.5)));
  const og::Quality3D q = og::quality_3d(fine);
  const double analytic = 0.5 * (kPi / 3) * (4.0 - 1.0) * 1.5;
  EXPECT_NEAR(total(q.volume) / analytic, 1.0, 1e-6);
  for (const double j : q.scaled_jacobian) {
    EXPECT_GT(j, 0.0);
  }
}

TEST(BlockTfi, RejectsMismatchedFaces) {
  og::BlockFaces f = og::block_faces(sector(5, 5, 5, 1.0, 2.0, 1.0, 1.0));
  f.kmax.at(2, 0) = f.kmax.at(2, 0) + Vec3{0.1, 0, 0};
  EXPECT_THROW(og::tfi_block(f), std::invalid_argument);
  f = og::block_faces(sector(5, 5, 5, 1.0, 2.0, 1.0, 1.0));
  f.imin = og::Grid2D(4, 5);
  EXPECT_THROW(og::tfi_block(f), std::invalid_argument);
}

TEST(Quality3D, UnitCubeCells) {
  og::Grid3D g(5, 5, 5);
  for (int k = 0; k < 5; ++k) {
    for (int j = 0; j < 5; ++j) {
      for (int i = 0; i < 5; ++i) {
        g.at(i, j, k) = {i / 4.0, j / 4.0, k / 4.0};
      }
    }
  }
  const og::Quality3D q = og::quality_3d(g);
  ASSERT_EQ(q.volume.size(), 64U);
  for (std::size_t c = 0; c < 64; ++c) {
    EXPECT_NEAR(q.scaled_jacobian[c], 1.0, 1e-14);
    EXPECT_NEAR(q.skewness[c], 0.0, 1e-12);
    EXPECT_NEAR(q.aspect_ratio[c], 1.0, 1e-14);
    EXPECT_NEAR(q.volume[c], 1.0 / 64.0, 1e-16);
    EXPECT_NEAR(q.non_orthogonality[c], 0.0, 1e-6);
    EXPECT_NEAR(q.volume_ratio[c], 1.0, 1e-14);
  }
}

TEST(Extrude, TranslateRotateNormal) {
  // Base: unit square in the xy plane, counter-clockwise about +z.
  og::Grid2D base(6, 4);
  for (int j = 0; j < 4; ++j) {
    for (int i = 0; i < 6; ++i) {
      base.at(i, j) = {i / 5.0, j / 3.0, 0.0};
    }
  }
  const std::vector<double> s = {0.0, 0.1, 0.3, 0.6, 1.0};
  const og::Extruded up = og::extrude_translate(base, {0, 0, 2}, s);
  EXPECT_FALSE(up.mirrored);
  EXPECT_NEAR(total(og::quality_3d(up.grid).volume), 2.0, 1e-14);
  const og::Extruded down = og::extrude_translate(base, {0, 0, -2}, s);
  EXPECT_TRUE(down.mirrored);
  const og::Quality3D qd = og::quality_3d(down.grid);
  EXPECT_NEAR(total(qd.volume), 2.0, 1e-14);
  for (const double j : qd.scaled_jacobian) {
    EXPECT_GT(j, 0.0);
  }
  // Rotating the square x in [1, 2] (y = 0 plane) about z by 90 degrees: quarter annulus volume.
  og::Grid2D radial(5, 5);
  for (int j = 0; j < 5; ++j) {
    for (int i = 0; i < 5; ++i) {
      radial.at(i, j) = {1.0 + i / 4.0, 0.0, j / 4.0};
    }
  }
  std::vector<double> t;
  for (int k = 0; k <= 400; ++k) {
    t.push_back(k / 400.0);
  }
  const og::Extruded turned = og::extrude_rotate(radial, {0, 0, 0}, {0, 0, 1}, kPi / 2, t);
  const og::Quality3D qr = og::quality_3d(turned.grid);
  EXPECT_NEAR(total(qr.volume) / (0.25 * kPi * 3.0), 1.0, 1e-5);
  for (const double j : qr.scaled_jacobian) {
    EXPECT_GT(j, 0.0);
  }
  const og::Extruded thick = og::extrude_normal(base, 0.5, s);
  EXPECT_NEAR(total(og::quality_3d(thick.grid).volume), 0.5, 1e-14);
  EXPECT_THROW(og::extrude_translate(base, {0, 0, 0}, s), std::invalid_argument);
  EXPECT_THROW(og::extrude_translate(base, {0, 0, 1}, {0.0, 0.5}), std::invalid_argument);
}

TEST(SurfaceTfi, CylinderPatchStaysOnTheSurface) {
  auto arc =
      std::make_shared<og::Arc>(Vec3{0, 0, 0}, 2.0, Vec3{0, 0, 1}, Vec3{1, 0, 0}, 0.0, kPi / 2);
  const og::ExtrudedSurface cylinder(arc, {0, 0, 3});
  const int ni = 21;
  const int nj = 11;
  og::Boundary2D b;
  for (int i = 0; i < ni; ++i) {
    const double t = 0.5 * kPi * std::pow(static_cast<double>(i) / (ni - 1), 1.5);
    b.bottom.push_back(cyl(2.0, t, 0.0));
    b.top.push_back(cyl(2.0, t, 3.0));
  }
  for (int j = 0; j < nj; ++j) {
    const double z = 3.0 * j / (nj - 1);
    b.left.push_back(cyl(2.0, 0.0, z));
    b.right.push_back(cyl(2.0, 0.5 * kPi, z));
  }
  const og::SurfaceTfi r = og::tfi_on_surface(cylinder, b);
  EXPECT_LT(r.max_boundary_distance, 1e-12);
  for (const Vec3& p : r.grid.x) {
    EXPECT_NEAR(std::hypot(p.x, p.y), 2.0, 1e-12);
  }
  EXPECT_EQ(og::norm(r.grid.at(5, 0) - b.bottom[5]), 0.0);
  EXPECT_EQ(og::norm(r.grid.at(ni - 1, 4) - b.right[4]), 0.0);
}
