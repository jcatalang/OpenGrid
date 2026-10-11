// Shared curve factories for the geometry tests.
#pragma once

#include <cmath>
#include <memory>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

#include "opengrid/core/geom/arc.hpp"
#include "opengrid/core/geom/cubic_spline.hpp"
#include "opengrid/core/geom/line.hpp"
#include "opengrid/core/geom/nurbs.hpp"
#include "opengrid/core/geom/polyline.hpp"

namespace og::test {

// Upper surface of a NACA 0012 (closed trailing edge), cosine-spaced, LE -> TE.
inline std::vector<Vec3> naca0012_upper(int n) {
  std::vector<Vec3> pts;
  for (int i = 0; i < n; ++i) {
    const double beta = std::numbers::pi * i / (n - 1);
    const double x = 0.5 * (1.0 - std::cos(beta));
    const double y = 0.6 * (0.2969 * std::sqrt(x) - 0.1260 * x - 0.3516 * x * x +
                            0.2843 * x * x * x - 0.1036 * x * x * x * x);
    pts.emplace_back(x, y, 0.0);
  }
  return pts;
}

inline std::vector<Vec3> wavy_points() {
  return {{0.0, 0.0, 0.0}, {0.7, 0.4, 0.1}, {1.5, -0.2, 0.0}, {2.0, 0.9, -0.3}, {3.1, 1.0, 0.2}};
}

// Exact full circle as a rational quadratic NURBS with 9 control points (Piegl & Tiller §7.5).
inline std::unique_ptr<Nurbs> nurbs_circle(const Vec3& c, double r) {
  const double w = std::numbers::sqrt2 / 2.0;
  std::vector<Vec3> ctrl = {{1, 0, 0},   {1, 1, 0},  {0, 1, 0},  {-1, 1, 0}, {-1, 0, 0},
                            {-1, -1, 0}, {0, -1, 0}, {1, -1, 0}, {1, 0, 0}};
  for (Vec3& p : ctrl) {
    p = c + r * p;
  }
  std::vector<double> weights = {1, w, 1, w, 1, w, 1, w, 1};
  std::vector<double> knots = {0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1};
  return std::make_unique<Nurbs>(2, knots, ctrl, weights);
}

inline std::unique_ptr<Nurbs> cubic_bspline() {
  std::vector<Vec3> ctrl = {{0, 0, 0}, {1, 2, 0}, {2, -1, 0.5}, {3, 3, 0}, {4, 0, 1}, {5, 1, 0}};
  std::vector<double> weights = {1, 1, 1, 1, 1, 1};
  std::vector<double> knots = {0, 0, 0, 0, 0.3, 0.6, 1, 1, 1, 1};
  return std::make_unique<Nurbs>(3, knots, ctrl, weights);
}

inline std::unique_ptr<Nurbs> rational_cubic() {
  std::vector<Vec3> ctrl = {{0, 0, 0}, {1, 2, 0}, {2, -1, 0.5}, {3, 3, 0}, {4, 0, 1}, {5, 1, 0}};
  std::vector<double> weights = {1, 2, 0.5, 1.5, 1, 3};
  std::vector<double> knots = {0, 0, 0, 0, 0.3, 0.6, 1, 1, 1, 1};
  return std::make_unique<Nurbs>(3, knots, ctrl, weights);
}

// Named collection of one curve of every type for parametrised tests.
inline std::vector<std::pair<std::string, std::unique_ptr<Curve>>> all_curves() {
  std::vector<std::pair<std::string, std::unique_ptr<Curve>>> curves;
  curves.emplace_back("line", std::make_unique<Line>(Vec3{0, 0, 0}, Vec3{2, 1, 0.5}));
  curves.emplace_back(
      "arc", std::make_unique<Arc>(Vec3{1, 2, 0}, 3.0, Vec3{0, 0, 1}, Vec3{1, 0, 0}, 0.3, 2.5));
  curves.emplace_back(
      "arc_cw", std::make_unique<Arc>(Vec3{0, 0, 0}, 0.5, Vec3{1, 1, 1}, Vec3{1, 0, 0}, 1.0, -2.0));
  curves.emplace_back("polyline", std::make_unique<Polyline>(wavy_points()));
  curves.emplace_back("spline_natural",
                      std::make_unique<CubicSpline>(wavy_points(), SplineMethod::Natural));
  curves.emplace_back("spline_clamped",
                      std::make_unique<CubicSpline>(wavy_points(), SplineMethod::Clamped,
                                                    Vec3{0, 1, 0}, Vec3{1, 0, 0}));
  curves.emplace_back("spline_catmull_rom",
                      std::make_unique<CubicSpline>(wavy_points(), SplineMethod::CatmullRom));
  curves.emplace_back("naca_upper",
                      std::make_unique<CubicSpline>(naca0012_upper(81), SplineMethod::Clamped,
                                                    Vec3{0, 1, 0}, std::nullopt));
  curves.emplace_back("nurbs_circle", nurbs_circle({0.5, 0, 0}, 20.0));
  curves.emplace_back("nurbs_cubic", cubic_bspline());
  curves.emplace_back("nurbs_rational", rational_cubic());
  return curves;
}

}  // namespace og::test
