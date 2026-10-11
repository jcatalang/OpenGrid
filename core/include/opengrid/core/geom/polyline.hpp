// Polyline through points P_0..P_n, parametrised by normalised cumulative chord length
// (constant speed equal to the total length). Every vertex is a breakpoint.
#pragma once

#include <vector>

#include "opengrid/core/geom/curve.hpp"

namespace og {

class Polyline final : public Curve {
 public:
  // Throws std::invalid_argument for fewer than 2 points or consecutive coincident points.
  explicit Polyline(std::vector<Vec3> points);

  [[nodiscard]] Vec3 eval(double t) const override;
  [[nodiscard]] Vec3 d1(double t) const override;
  [[nodiscard]] Vec3 d2(double t) const override;
  [[nodiscard]] std::vector<double> breakpoints() const override;

 private:
  std::vector<Vec3> points_;
  std::vector<double> knots_;
};

// Normalised cumulative chord-length parameters of a point sequence (shared with splines).
// Throws std::invalid_argument naming the first zero-length chord.
std::vector<double> chord_length_knots(const std::vector<Vec3>& points, const char* what);

}  // namespace og
