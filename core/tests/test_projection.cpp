#include <gtest/gtest.h>

#include <cmath>
#include <numbers>

#include "curve_fixtures.hpp"

using og::Vec3;

TEST(Projection, ArcAnalytic) {
  const og::Arc arc({0, 0, 0}, 2.0, {0, 0, 1}, {1, 0, 0}, 0.0, std::numbers::pi);
  // Off-curve point inside and outside the circle, at angle 1 rad.
  for (const double radius : {0.5, 1.9, 2.0, 3.7}) {
    const Vec3 p{radius * std::cos(1.0), radius * std::sin(1.0), 0.0};
    const og::CurveProjection proj = arc.project(p);
    EXPECT_NEAR(proj.t, 1.0 / std::numbers::pi, 1e-12) << radius;
    EXPECT_NEAR(proj.distance, std::fabs(radius - 2.0), 1e-12 * arc.length()) << radius;
  }
  // Beyond the end of the arc the end point is closest.
  const og::CurveProjection end = arc.project({-3.0, -1.0, 0.0});
  EXPECT_EQ(end.t, 1.0);
  EXPECT_NEAR(end.distance, std::hypot(1.0, 1.0), 1e-14);
}

TEST(Projection, PointsOnCurves) {
  for (const auto& [name, curve] : og::test::all_curves()) {
    const double scale = curve->length();
    for (int k = 0; k <= 37; ++k) {
      const double t = k / 37.0;
      const Vec3 p = curve->eval(t);
      const og::CurveProjection proj = curve->project(p);
      EXPECT_LT(proj.distance, 1e-12 * scale) << name << " t=" << t;
    }
  }
}

TEST(Projection, OffsetAlongNormal) {
  // Offset curve points along a direction orthogonal to the tangent by less than the local
  // radius of curvature: the projection must return the same foot point and the offset.
  for (const auto& [name, curve] : og::test::all_curves()) {
    const double scale = curve->length();
    for (int k = 1; k < 23; ++k) {
      const double t = k / 23.0;
      const Vec3 x1 = curve->d1(t);
      const Vec3 x2 = curve->d2(t);
      const double speed = og::norm(x1);
      const Vec3 tangent = x1 / speed;
      // Unit normal in the osculating plane, or any perpendicular for straight pieces.
      Vec3 normal = x2 - og::dot(x2, tangent) * tangent;
      double curvature = og::norm(normal) / (speed * speed);
      if (og::norm(normal) <= 1e-12 * og::norm(x2) + 1e-300) {
        normal = og::cross(tangent, {0, 0, 1});
        if (og::norm(normal) < 1e-6) {
          normal = og::cross(tangent, {1, 0, 0});
        }
        curvature = 0.0;
      }
      normal = normal / og::norm(normal);
      double offset = 1e-3 * scale;
      if (curvature > 0.0) {
        offset = std::min(offset, 0.1 / curvature);
      }
      // Offset to the convex side (away from the centre of curvature).
      const Vec3 foot = curve->eval(t);
      const Vec3 p = foot - offset * normal;
      const og::CurveProjection proj = curve->project(p);
      EXPECT_NEAR(proj.distance, offset, 1e-12 * scale) << name << " t=" << t;
      EXPECT_LT(og::distance(proj.point, foot), 1e-9 * scale) << name << " t=" << t;
    }
  }
}

TEST(Tessellation, ChordDeviationWithinTolerance) {
  for (const auto& [name, curve] : og::test::all_curves()) {
    const double tol = 1e-4;
    const std::vector<double> ts = curve->tessellate(tol);
    ASSERT_GE(ts.size(), 2U) << name;
    EXPECT_EQ(ts.front(), 0.0);
    EXPECT_EQ(ts.back(), 1.0);
    for (std::size_t i = 1; i < ts.size(); ++i) {
      ASSERT_LT(ts[i - 1], ts[i]) << name;
      const double mid = 0.5 * (ts[i - 1] + ts[i]);
      const Vec3 chord_mid = 0.5 * (curve->eval(ts[i - 1]) + curve->eval(ts[i]));
      EXPECT_LT(og::distance(curve->eval(mid), chord_mid), 2.0 * tol * curve->length()) << name;
    }
  }
}
