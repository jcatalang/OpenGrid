#include "opengrid/core/geom/arc.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace og {

Arc::Arc(const Vec3& center, double radius, const Vec3& normal, const Vec3& x_axis,
         double start_angle, double end_angle)
    : center_(center), radius_(radius), a0_(start_angle), sweep_(end_angle - start_angle) {
  if (!(radius > 0.0)) {
    throw std::invalid_argument("arc: radius must be positive");
  }
  const double n_len = norm(normal);
  if (n_len == 0.0) {
    throw std::invalid_argument("arc: normal must be non-zero");
  }
  const Vec3 n = normal / n_len;
  const Vec3 u = x_axis - dot(x_axis, n) * n;
  const double u_len = norm(u);
  if (u_len <= 1e-12 * norm(x_axis) || u_len == 0.0) {
    throw std::invalid_argument("arc: x_axis must not be parallel to the normal");
  }
  u_ = u / u_len;
  v_ = cross(n, u_);
  if (sweep_ == 0.0) {
    throw std::invalid_argument("arc: start and end angles coincide");
  }
  if (std::fabs(sweep_) > 2.0 * std::numbers::pi * (1.0 + 1e-12)) {
    throw std::invalid_argument("arc: sweep exceeds a full circle");
  }
}

Vec3 Arc::eval(double t) const {
  const double theta = a0_ + t * sweep_;
  return center_ + radius_ * (std::cos(theta) * u_ + std::sin(theta) * v_);
}

Vec3 Arc::d1(double t) const {
  const double theta = a0_ + t * sweep_;
  return (radius_ * sweep_) * (-std::sin(theta) * u_ + std::cos(theta) * v_);
}

Vec3 Arc::d2(double t) const {
  const double theta = a0_ + t * sweep_;
  return (-radius_ * sweep_ * sweep_) * (std::cos(theta) * u_ + std::sin(theta) * v_);
}

}  // namespace og
