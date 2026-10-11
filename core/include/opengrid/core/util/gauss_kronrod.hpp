// Adaptive Gauss–Kronrod (G7/K15) quadrature.
//
// Method: 15-point Kronrod rule with embedded 7-point Gauss rule; the difference |K15 - G7| is
// used as error estimate and intervals are bisected recursively until it falls below the
// tolerance. Reference: R. Piessens et al., QUADPACK (1983), routine QK15.
#pragma once

#include <functional>

namespace og {

struct QuadratureOptions {
  double rel_tol = 1e-14;  // relative tolerance per subinterval
  double abs_tol = 0.0;    // absolute tolerance per subinterval
  int max_depth = 30;      // maximum bisection depth
};

// Integral of f over [a, b] (a > b gives the negated integral).
double integrate_gk15(const std::function<double(double)>& f, double a, double b,
                      const QuadratureOptions& options = {});

}  // namespace og
