// Non-uniform rational B-spline curve of degree p with control points P_i, weights w_i > 0
// and knot vector U (size n_ctrl + p + 1). The valid parameter domain [U_p, U_{n_ctrl}] is
// mapped linearly onto t in [0, 1]; distinct interior knots are the breakpoints.
//
// Method: basis functions and their derivatives by Piegl & Tiller algorithm A2.3; rational
// derivatives from the homogeneous ones (A4.2). Reference: L. Piegl & W. Tiller, The NURBS
// Book, 2nd ed. (1997), ch. 2 and 4.
#pragma once

#include <vector>

#include "opengrid/core/geom/curve.hpp"

namespace og {

class Nurbs final : public Curve {
 public:
  // Throws std::invalid_argument on inconsistent sizes, decreasing knots, interior knot
  // multiplicity > degree, non-positive weights, or an empty parameter domain.
  Nurbs(int degree, std::vector<double> knots, std::vector<Vec3> control_points,
        std::vector<double> weights);

  [[nodiscard]] Vec3 eval(double t) const override;
  [[nodiscard]] Vec3 d1(double t) const override;
  [[nodiscard]] Vec3 d2(double t) const override;
  [[nodiscard]] std::vector<double> breakpoints() const override;

 private:
  struct Derivatives {
    Vec3 c0;
    Vec3 c1;
    Vec3 c2;
  };
  [[nodiscard]] Derivatives evaluate(double t) const;
  // Position and derivatives up to order `wanted` (0: position only; the others are zero).
  [[nodiscard]] Derivatives evaluate_up_to(double t, int wanted) const;
  [[nodiscard]] int find_span(double u) const;

  int degree_;
  std::vector<double> knots_;
  std::vector<Vec3> control_points_;
  std::vector<double> weights_;
  double u_lo_ = 0.0;
  double u_hi_ = 0.0;
};

}  // namespace og
