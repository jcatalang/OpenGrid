// Circular arc in the plane through `center` with unit normal n:
//   x(t) = c + r (cos(theta) u + sin(theta) v),  theta = a0 + t (a1 - a0),
// where u is `x_axis` made orthogonal to n and normalised, and v = n x u.
// Angles are in radians; a1 < a0 traverses clockwise about n. |a1 - a0| <= 2 pi.
#pragma once

#include "opengrid/core/geom/curve.hpp"

namespace og {

class Arc final : public Curve {
 public:
  // Throws std::invalid_argument on non-positive radius, zero normal, x_axis parallel to the
  // normal, a1 == a0, or |a1 - a0| > 2 pi.
  Arc(const Vec3& center, double radius, const Vec3& normal, const Vec3& x_axis, double start_angle,
      double end_angle);

  [[nodiscard]] Vec3 eval(double t) const override;
  [[nodiscard]] Vec3 d1(double t) const override;
  [[nodiscard]] Vec3 d2(double t) const override;

 private:
  Vec3 center_;
  double radius_;
  Vec3 u_;
  Vec3 v_;
  double a0_;
  double sweep_;
};

}  // namespace og
