// Bracketed scalar root finding.
//
// Method: Brent's method (inverse quadratic interpolation / secant with bisection safeguard).
// Reference: R. P. Brent, Algorithms for Minimization without Derivatives (1973), ch. 4;
// W. H. Press et al., Numerical Recipes, 3rd ed., §9.3.
#pragma once

#include <functional>
#include <optional>

namespace og {

// Root of f in [a, b]; requires f(a) and f(b) of opposite sign (or zero). Returns nullopt if
// the interval does not bracket a root. x_tol is an absolute tolerance on the root (0 = to
// machine precision).
std::optional<double> brent_root(const std::function<double(double)>& f, double a, double b,
                                 double x_tol = 0.0, int max_iterations = 300);

}  // namespace og
