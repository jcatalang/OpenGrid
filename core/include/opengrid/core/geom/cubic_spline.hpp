// Interpolating C1/C2 cubic spline through points P_0..P_n in piecewise cubic Hermite form,
// parametrised by normalised cumulative chord length t_i.
//
// Methods (tangents m_i = dx/dt at the points):
//  * Natural:    C2 everywhere, x'' = 0 at both ends (tridiagonal system, Thomas algorithm).
//  * Clamped:    C2 everywhere; prescribed end tangent direction(s). An end without a
//                prescribed tangent uses the natural condition. The direction is scaled to
//                the chord length (|dx/dt| ~ total chord length for chord parametrisation).
//  * CatmullRom: C1; m_i = (P_{i+1} - P_{i-1}) / (t_{i+1} - t_{i-1}) (non-uniform Catmull-Rom),
//                end tangents chosen so that x'' = 0 at the ends.
// Reference: C. de Boor, A Practical Guide to Splines (1978), ch. IV; E. Catmull & R. Rom (1974).
#pragma once

#include <optional>
#include <vector>

#include "opengrid/core/geom/curve.hpp"

namespace og {

enum class SplineMethod { Natural, Clamped, CatmullRom };

class CubicSpline final : public Curve {
 public:
  // Throws std::invalid_argument for < 2 points, coincident consecutive points, zero tangents,
  // tangents given for a non-clamped method, or Clamped without any tangent.
  CubicSpline(std::vector<Vec3> points, SplineMethod method,
              std::optional<Vec3> start_tangent = std::nullopt,
              std::optional<Vec3> end_tangent = std::nullopt);

  [[nodiscard]] Vec3 eval(double t) const override;
  [[nodiscard]] Vec3 d1(double t) const override;
  [[nodiscard]] Vec3 d2(double t) const override;
  [[nodiscard]] std::vector<double> breakpoints() const override;

  [[nodiscard]] const std::vector<double>& knots() const {
    return knots_;
  }
  [[nodiscard]] const std::vector<Vec3>& tangents() const {
    return tangents_;
  }

 private:
  void solve_c2(const std::optional<Vec3>& start, const std::optional<Vec3>& end);
  void catmull_rom();

  std::vector<Vec3> points_;
  std::vector<double> knots_;
  std::vector<Vec3> tangents_;
};

}  // namespace og
