// Globally adaptive Gauss–Kronrod (G7/K15) quadrature.
//
// Method: 15-point Kronrod rule with embedded 7-point Gauss rule; |K15 - G7| is the error
// estimate. The subinterval with the largest error is bisected until the total error is
// below max(rel_tol |I|, abs_tol, 50 eps * integral of |f|) or max_intervals is reached, so
// the cost stays bounded even for noisy integrands. Ties are broken by position, making the
// result deterministic. Reference: R. Piessens et al., QUADPACK (1983), routines QK15/QAG.
#pragma once

#include <functional>

namespace og {

struct QuadratureOptions {
  double rel_tol = 1e-14;   // relative tolerance per subinterval
  double abs_tol = 0.0;     // absolute tolerance per subinterval
  int max_intervals = 400;  // maximum number of subintervals
};

// Integral of f over [a, b] (a > b gives the negated integral).
double integrate_gk15(const std::function<double(double)>& f, double a, double b,
                      const QuadratureOptions& options = {});

}  // namespace og
