// Abstract parametric curve x(t), t in [0, 1], and the algorithms shared by all curve types:
// arc length, arc-length reparametrization, closest-point projection and tessellation.
//
// Every curve is C2 inside each *segment* between consecutive breakpoints (knots of piecewise
// curves); algorithms integrate / search segment by segment so that derivative jumps at
// breakpoints never degrade accuracy.
#pragma once

#include <cstddef>
#include <mutex>
#include <vector>

#include "opengrid/core/geom/vec3.hpp"

namespace og {

struct CurveProjection {
  double t = 0.0;  // parameter of the closest point
  Vec3 point;      // closest point on the curve
  double distance = 0.0;
};

class Curve {
 public:
  Curve() = default;
  Curve(const Curve&) = delete;
  Curve& operator=(const Curve&) = delete;
  virtual ~Curve() = default;

  // Position and first / second derivatives with respect to t. At a breakpoint, the
  // derivatives of the segment to the right are returned (left one at t = 1).
  [[nodiscard]] virtual Vec3 eval(double t) const = 0;
  [[nodiscard]] virtual Vec3 d1(double t) const = 0;
  [[nodiscard]] virtual Vec3 d2(double t) const = 0;

  // Sorted parameters 0 = b_0 < b_1 < ... < b_k = 1 delimiting the smooth segments.
  [[nodiscard]] virtual std::vector<double> breakpoints() const;

  // Total arc length.
  [[nodiscard]] double length() const;
  // Signed arc length between t0 and t1 (negative if t1 < t0).
  [[nodiscard]] double length_between(double t0, double t1) const;
  // Parameter t with length_between(0, t) == s; s is clamped to [0, length()].
  // Accurate to 1e-13 * length() in arc length.
  [[nodiscard]] double t_at_length(double s) const;
  // Closest point on the curve to p (global search over all segments, then Newton polish).
  [[nodiscard]] CurveProjection project(const Vec3& p) const;
  // Adaptive polyline approximation: chord deviation <= rel_tol * length(). Returns the
  // parameters of the vertices (including 0 and 1).
  [[nodiscard]] std::vector<double> tessellate(double rel_tol) const;

 private:
  void build_length_table() const;
  [[nodiscard]] double speed(double t) const;
  [[nodiscard]] double polish_projection(const Vec3& p, double t, double lo, double hi) const;

  mutable std::once_flag length_once_;
  mutable std::vector<double> segment_params_;  // breakpoints
  mutable std::vector<double> cumulative_;      // arc length at each breakpoint
};

// Index i of the segment [knots[i], knots[i+1]] containing t, clamped to [0, knots.size()-2].
std::size_t find_segment(const std::vector<double>& knots, double t);

}  // namespace og
