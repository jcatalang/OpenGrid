// Straight line segment x(t) = p0 + t (p1 - p0).
#pragma once

#include "opengrid/core/geom/curve.hpp"

namespace og {

class Line final : public Curve {
 public:
  // Throws std::invalid_argument if p0 == p1.
  Line(const Vec3& p0, const Vec3& p1);

  [[nodiscard]] Vec3 eval(double t) const override;
  [[nodiscard]] Vec3 d1(double t) const override;
  [[nodiscard]] Vec3 d2(double t) const override;

  [[nodiscard]] const Vec3& p0() const {
    return p0_;
  }
  [[nodiscard]] const Vec3& p1() const {
    return p1_;
  }

 private:
  Vec3 p0_;
  Vec3 p1_;
};

}  // namespace og
