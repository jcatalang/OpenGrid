#include <gtest/gtest.h>

#include <cmath>
#include <numbers>
#include <stdexcept>

#include "curve_fixtures.hpp"

using og::Vec3;

namespace {

// Parameters well inside each smooth segment (finite differences must not straddle a knot).
std::vector<double> interior_samples(const og::Curve& curve) {
  std::vector<double> ts;
  const std::vector<double> bps = curve.breakpoints();
  for (std::size_t i = 0; i + 1 < bps.size(); ++i) {
    for (const double frac : {0.13, 0.5, 0.91}) {
      ts.push_back(bps[i] + frac * (bps[i + 1] - bps[i]));
    }
  }
  return ts;
}

}  // namespace

TEST(Curves, DerivativesMatchFiniteDifferences) {
  for (const auto& [name, curve] : og::test::all_curves()) {
    const std::vector<double> bps = curve->breakpoints();
    for (const double t : interior_samples(*curve)) {
      const std::size_t seg = og::find_segment(bps, t);
      const double h = 1e-4 * (bps[seg + 1] - bps[seg]);
      const Vec3 fd1 = (curve->eval(t + h) - curve->eval(t - h)) / (2.0 * h);
      const Vec3 fd2 = (curve->d1(t + h) - curve->d1(t - h)) / (2.0 * h);
      const double scale1 = og::norm(curve->d1(t)) + 1e-12;
      const double scale2 = og::norm(curve->d2(t)) + scale1;
      EXPECT_LT(og::norm(curve->d1(t) - fd1), 1e-6 * scale1) << name << " t=" << t;
      EXPECT_LT(og::norm(curve->d2(t) - fd2), 1e-5 * scale2) << name << " t=" << t;
    }
  }
}

TEST(Curves, EndpointsAndBreakpoints) {
  for (const auto& [name, curve] : og::test::all_curves()) {
    const std::vector<double> bps = curve->breakpoints();
    ASSERT_GE(bps.size(), 2U) << name;
    EXPECT_EQ(bps.front(), 0.0) << name;
    EXPECT_EQ(bps.back(), 1.0) << name;
    for (std::size_t i = 1; i < bps.size(); ++i) {
      EXPECT_LT(bps[i - 1], bps[i]) << name;
    }
  }
}

TEST(Line, EvaluationAndLength) {
  const og::Line line({1, 2, 3}, {4, 6, 3});
  EXPECT_DOUBLE_EQ(line.eval(0.5).x, 2.5);
  EXPECT_DOUBLE_EQ(line.length(), 5.0);
  EXPECT_THROW(og::Line({1, 1, 1}, {1, 1, 1}), std::invalid_argument);
}

TEST(Arc, PointsLieOnCircle) {
  const og::Arc arc({1, -2, 0.5}, 2.5, {0, 0, 2}, {3, 0, 1}, -0.4, 1.7);
  for (int k = 0; k <= 20; ++k) {
    const double t = k / 20.0;
    EXPECT_NEAR(og::distance(arc.eval(t), {1, -2, 0.5}), 2.5, 1e-14);
    EXPECT_NEAR(arc.eval(t).z, 0.5, 1e-14);
  }
  // x_axis is projected into the plane: start point at angle -0.4 from +x.
  EXPECT_NEAR(arc.eval(0.0).x, 1.0 + 2.5 * std::cos(-0.4), 1e-14);
  EXPECT_NEAR(arc.eval(0.0).y, -2.0 + 2.5 * std::sin(-0.4), 1e-14);
}

TEST(Arc, InvalidInput) {
  EXPECT_THROW(og::Arc({0, 0, 0}, 0.0, {0, 0, 1}, {1, 0, 0}, 0, 1), std::invalid_argument);
  EXPECT_THROW(og::Arc({0, 0, 0}, 1.0, {0, 0, 0}, {1, 0, 0}, 0, 1), std::invalid_argument);
  EXPECT_THROW(og::Arc({0, 0, 0}, 1.0, {0, 0, 1}, {0, 0, 3}, 0, 1), std::invalid_argument);
  EXPECT_THROW(og::Arc({0, 0, 0}, 1.0, {0, 0, 1}, {1, 0, 0}, 1, 1), std::invalid_argument);
  EXPECT_THROW(og::Arc({0, 0, 0}, 1.0, {0, 0, 1}, {1, 0, 0}, 0, 7), std::invalid_argument);
}

TEST(Polyline, VerticesAndInvalidInput) {
  const og::Polyline poly({{0, 0, 0}, {3, 0, 0}, {3, 4, 0}});
  EXPECT_DOUBLE_EQ(poly.length(), 7.0);
  EXPECT_NEAR(poly.eval(3.0 / 7.0).x, 3.0, 1e-15);
  EXPECT_NEAR(poly.eval(3.0 / 7.0).y, 0.0, 1e-15);
  EXPECT_THROW(og::Polyline({{0, 0, 0}}), std::invalid_argument);
  EXPECT_THROW(og::Polyline({{0, 0, 0}, {1, 0, 0}, {1, 0, 0}}), std::invalid_argument);
}

TEST(CubicSpline, InterpolatesPoints) {
  const std::vector<Vec3> pts = og::test::wavy_points();
  for (const auto method : {og::SplineMethod::Natural, og::SplineMethod::CatmullRom}) {
    const og::CubicSpline spline(pts, method);
    for (std::size_t i = 0; i < pts.size(); ++i) {
      EXPECT_LT(og::distance(spline.eval(spline.knots()[i]), pts[i]), 1e-14);
    }
  }
}

TEST(CubicSpline, NaturalIsC2WithZeroEndCurvature) {
  const og::CubicSpline spline(og::test::wavy_points(), og::SplineMethod::Natural);
  EXPECT_LT(og::norm(spline.d2(0.0)), 1e-12);
  EXPECT_LT(og::norm(spline.d2(1.0)), 1e-12);
  const auto& knots = spline.knots();
  for (std::size_t i = 1; i + 1 < knots.size(); ++i) {
    const double eps = 1e-10;
    const Vec3 left = spline.d2(knots[i] - eps);
    const Vec3 right = spline.d2(knots[i]);
    EXPECT_LT(og::norm(left - right), 1e-6 * (og::norm(right) + 1.0)) << "knot " << i;
  }
}

TEST(CubicSpline, ClampedTangentDirections) {
  const og::CubicSpline spline(og::test::wavy_points(), og::SplineMethod::Clamped, Vec3{0, 2, 0},
                               Vec3{1, 1, 0});
  const Vec3 start = spline.d1(0.0) / og::norm(spline.d1(0.0));
  const Vec3 end = spline.d1(1.0) / og::norm(spline.d1(1.0));
  EXPECT_NEAR(start.y, 1.0, 1e-14);
  EXPECT_NEAR(end.x, std::numbers::sqrt2 / 2.0, 1e-14);
  EXPECT_NEAR(end.y, std::numbers::sqrt2 / 2.0, 1e-14);
}

TEST(CubicSpline, InvalidInput) {
  const std::vector<Vec3> pts = og::test::wavy_points();
  EXPECT_THROW(og::CubicSpline(pts, og::SplineMethod::Clamped), std::invalid_argument);
  EXPECT_THROW(og::CubicSpline(pts, og::SplineMethod::Natural, Vec3{1, 0, 0}),
               std::invalid_argument);
  EXPECT_THROW(og::CubicSpline(pts, og::SplineMethod::Clamped, Vec3{0, 0, 0}),
               std::invalid_argument);
  EXPECT_THROW(og::CubicSpline({{0, 0, 0}, {0, 0, 0}}, og::SplineMethod::Natural),
               std::invalid_argument);
}

TEST(CubicSpline, TwoPointsIsStraight) {
  const og::CubicSpline spline({{0, 0, 0}, {2, 0, 0}}, og::SplineMethod::Natural);
  EXPECT_NEAR(spline.eval(0.25).x, 0.5, 1e-15);
  EXPECT_NEAR(spline.length(), 2.0, 1e-14);
}

TEST(Nurbs, CircleIsExact) {
  const auto circle = og::test::nurbs_circle({0.5, 0, 0}, 20.0);
  for (int k = 0; k <= 400; ++k) {
    const double t = k / 400.0;
    EXPECT_NEAR(og::distance(circle->eval(t), {0.5, 0, 0}), 20.0, 1e-13);
  }
  EXPECT_LT(og::distance(circle->eval(0.0), circle->eval(1.0)), 1e-14);
}

TEST(Nurbs, NonRationalReproducesControlPolygonEnds) {
  const auto curve = og::test::cubic_bspline();
  EXPECT_LT(og::distance(curve->eval(0.0), {0, 0, 0}), 1e-15);
  EXPECT_LT(og::distance(curve->eval(1.0), {5, 1, 0}), 1e-15);
  // Clamped end tangent is p (P1 - P0) / (u1 - u0) with the domain scaled to [0, 1].
  const Vec3 tangent = curve->d1(0.0);
  EXPECT_NEAR(tangent.x, 3.0 * 1.0 / 0.3, 1e-12);
  EXPECT_NEAR(tangent.y, 3.0 * 2.0 / 0.3, 1e-12);
}

TEST(Nurbs, LinearDegreeHasZeroSecondDerivative) {
  const og::Nurbs poly(1, {0, 0, 0.5, 1, 1}, {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}}, {1, 1, 1});
  EXPECT_NEAR(poly.eval(0.75).y, 0.5, 1e-15);
  EXPECT_EQ(og::norm(poly.d2(0.3)), 0.0);
  EXPECT_NEAR(poly.length(), 2.0, 1e-14);
}

TEST(Nurbs, InvalidInput) {
  const std::vector<Vec3> ctrl = {{0, 0, 0}, {1, 1, 0}, {2, 0, 0}};
  EXPECT_THROW(og::Nurbs(0, {0, 1, 2}, ctrl, {1, 1, 1}), std::invalid_argument);
  EXPECT_THROW(og::Nurbs(2, {0, 0, 0, 1, 1}, ctrl, {1, 1, 1}), std::invalid_argument);
  EXPECT_THROW(og::Nurbs(2, {0, 0, 0, 1, 1, 1}, ctrl, {1, 0, 1}), std::invalid_argument);
  EXPECT_THROW(og::Nurbs(2, {0, 0, 1, 0, 1, 1}, ctrl, {1, 1, 1}), std::invalid_argument);
  EXPECT_THROW(og::Nurbs(2, {0, 0, 0, 1, 1, 1}, ctrl, {1, 1}), std::invalid_argument);
}
