#include "opengrid/core/geom/line.hpp"

#include <stdexcept>

namespace og {

Line::Line(const Vec3& p0, const Vec3& p1) : p0_(p0), p1_(p1) {
  if (norm(p1 - p0) == 0.0) {
    throw std::invalid_argument("line: end points coincide");
  }
}

Vec3 Line::eval(double t) const {
  return p0_ + t * (p1_ - p0_);
}

Vec3 Line::d1(double /*t*/) const {
  return p1_ - p0_;
}

Vec3 Line::d2(double /*t*/) const {
  return {};
}

}  // namespace og
