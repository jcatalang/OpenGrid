#include "opengrid/core/geom/cubic_spline.hpp"

#include <stdexcept>
#include <string>
#include <utility>

#include "opengrid/core/geom/polyline.hpp"

namespace og {
namespace {

Vec3 scaled_direction(const Vec3& tangent, double scale, const char* which) {
  const double len = norm(tangent);
  if (len == 0.0) {
    throw std::invalid_argument(std::string("spline: ") + which + " tangent must be non-zero");
  }
  return tangent * (scale / len);
}

}  // namespace

CubicSpline::CubicSpline(std::vector<Vec3> points, SplineMethod method,
                         std::optional<Vec3> start_tangent, std::optional<Vec3> end_tangent)
    : points_(std::move(points)), knots_(chord_length_knots(points_, "spline")) {
  double chord = 0.0;
  for (std::size_t i = 1; i < points_.size(); ++i) {
    chord += distance(points_[i], points_[i - 1]);
  }
  std::optional<Vec3> start;
  std::optional<Vec3> end;
  if (method == SplineMethod::Clamped) {
    if (!start_tangent && !end_tangent) {
      throw std::invalid_argument("spline: clamped method needs a start and/or end tangent");
    }
    if (start_tangent) {
      start = scaled_direction(*start_tangent, chord, "start");
    }
    if (end_tangent) {
      end = scaled_direction(*end_tangent, chord, "end");
    }
  } else if (start_tangent || end_tangent) {
    throw std::invalid_argument("spline: tangents are only allowed with the clamped method");
  }

  if (method == SplineMethod::CatmullRom) {
    catmull_rom();
  } else {
    solve_c2(start, end);
  }
}

void CubicSpline::solve_c2(const std::optional<Vec3>& start, const std::optional<Vec3>& end) {
  // Tridiagonal system a_i m_{i-1} + b_i m_i + c_i m_{i+1} = r_i for the tangents m_i.
  const std::size_t n = points_.size();
  std::vector<double> a(n, 0.0);
  std::vector<double> b(n, 0.0);
  std::vector<double> c(n, 0.0);
  std::vector<Vec3> r(n);
  std::vector<double> h(n - 1);
  std::vector<Vec3> d(n - 1);
  for (std::size_t i = 0; i + 1 < n; ++i) {
    h[i] = knots_[i + 1] - knots_[i];
    d[i] = (points_[i + 1] - points_[i]) / h[i];
  }
  if (start) {
    b[0] = 1.0;
    r[0] = *start;
  } else {
    b[0] = 2.0;
    c[0] = 1.0;
    r[0] = 3.0 * d[0];
  }
  for (std::size_t i = 1; i + 1 < n; ++i) {
    a[i] = h[i];
    b[i] = 2.0 * (h[i - 1] + h[i]);
    c[i] = h[i - 1];
    r[i] = 3.0 * (h[i] * d[i - 1] + h[i - 1] * d[i]);
  }
  if (end) {
    b[n - 1] = 1.0;
    r[n - 1] = *end;
  } else {
    a[n - 1] = 1.0;
    b[n - 1] = 2.0;
    r[n - 1] = 3.0 * d[n - 2];
  }
  // Thomas algorithm (diagonally dominant for both end conditions).
  for (std::size_t i = 1; i < n; ++i) {
    const double w = a[i] / b[i - 1];
    b[i] -= w * c[i - 1];
    r[i] -= w * r[i - 1];
  }
  tangents_.assign(n, Vec3{});
  tangents_[n - 1] = r[n - 1] / b[n - 1];
  for (std::size_t k = n - 1; k-- > 0;) {
    tangents_[k] = (r[k] - c[k] * tangents_[k + 1]) / b[k];
  }
}

void CubicSpline::catmull_rom() {
  const std::size_t n = points_.size();
  tangents_.assign(n, Vec3{});
  if (n == 2) {
    const Vec3 d = (points_[1] - points_[0]) / (knots_[1] - knots_[0]);
    tangents_[0] = d;
    tangents_[1] = d;
    return;
  }
  for (std::size_t i = 1; i + 1 < n; ++i) {
    tangents_[i] = (points_[i + 1] - points_[i - 1]) / (knots_[i + 1] - knots_[i - 1]);
  }
  // Zero second derivative at the ends: m_0 = (3 d_0 - m_1) / 2, and symmetrically.
  const Vec3 d_first = (points_[1] - points_[0]) / (knots_[1] - knots_[0]);
  const Vec3 d_last = (points_[n - 1] - points_[n - 2]) / (knots_[n - 1] - knots_[n - 2]);
  tangents_[0] = 0.5 * (3.0 * d_first - tangents_[1]);
  tangents_[n - 1] = 0.5 * (3.0 * d_last - tangents_[n - 2]);
}

Vec3 CubicSpline::eval(double t) const {
  const std::size_t i = find_segment(knots_, t);
  const double h = knots_[i + 1] - knots_[i];
  const double u = (t - knots_[i]) / h;
  const double u2 = u * u;
  const double u3 = u2 * u;
  const double h00 = 2.0 * u3 - 3.0 * u2 + 1.0;
  const double h10 = u3 - 2.0 * u2 + u;
  const double h01 = -2.0 * u3 + 3.0 * u2;
  const double h11 = u3 - u2;
  return h00 * points_[i] + (h10 * h) * tangents_[i] + h01 * points_[i + 1] +
         (h11 * h) * tangents_[i + 1];
}

Vec3 CubicSpline::d1(double t) const {
  const std::size_t i = find_segment(knots_, t);
  const double h = knots_[i + 1] - knots_[i];
  const double u = (t - knots_[i]) / h;
  const double u2 = u * u;
  const double h00 = 6.0 * u2 - 6.0 * u;
  const double h10 = 3.0 * u2 - 4.0 * u + 1.0;
  const double h01 = -6.0 * u2 + 6.0 * u;
  const double h11 = 3.0 * u2 - 2.0 * u;
  return (h00 / h) * points_[i] + h10 * tangents_[i] + (h01 / h) * points_[i + 1] +
         h11 * tangents_[i + 1];
}

Vec3 CubicSpline::d2(double t) const {
  const std::size_t i = find_segment(knots_, t);
  const double h = knots_[i + 1] - knots_[i];
  const double u = (t - knots_[i]) / h;
  const double h00 = 12.0 * u - 6.0;
  const double h10 = 6.0 * u - 4.0;
  const double h01 = -12.0 * u + 6.0;
  const double h11 = 6.0 * u - 2.0;
  return (h00 / (h * h)) * points_[i] + (h10 / h) * tangents_[i] +
         (h01 / (h * h)) * points_[i + 1] + (h11 / h) * tangents_[i + 1];
}

std::vector<double> CubicSpline::breakpoints() const {
  return knots_;
}

}  // namespace og
