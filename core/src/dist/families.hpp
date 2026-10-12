// Internal: raw stretching functions and numerically stable helpers shared by the
// distribution implementations. Not part of the public API.
#pragma once

#include <cmath>
#include <limits>
#include <numbers>
#include <vector>

namespace og::dist_detail {

constexpr double kNegInf = -std::numeric_limits<double>::infinity();

// log(sinh(x)) for x >= 0, without overflow.
inline double log_sinh(double x) {
  if (x <= 0.0) {
    return kNegInf;
  }
  if (x > 20.0) {
    return x - std::numbers::ln2 + std::log1p(-std::exp(-2.0 * x));
  }
  return std::log(std::sinh(x));
}

// log(cosh(x)), without overflow.
inline double log_cosh(double x) {
  const double ax = std::fabs(x);
  return ax + std::log1p(std::exp(-2.0 * ax)) - std::numbers::ln2;
}

// Raw families evaluated at xi_i = i / n, i = 0..n, written into s (size n + 1).
// "Start derivative" S0 = ds/dxi at xi = 0, S1 = ds/dxi at xi = 1.
// They return false when the family cannot represent the request.
bool tanh_one_sided(double s0_derivative, int n, std::vector<double>& s);
bool sine_one_sided(double s0_derivative, int n, std::vector<double>& s);
// Two-sided families also return the last interval s_n - s_{n-1} computed without
// cancellation (from the complement 1 - s), which the matching iteration needs.
bool tanh_two_sided(double s0_derivative, double s1_derivative, int n, std::vector<double>& s,
                    double& last_spacing);
bool sine_two_sided(double s0_derivative, double s1_derivative, int n, std::vector<double>& s,
                    double& last_spacing);

// Geometric progression with first interval d and ratio r = exp(log_ratio); returns false if
// the request is outside floating-point range.
bool geometric_from_begin(double d, int n, std::vector<double>& s);
bool geometric_from_ratio(double ratio, int n, std::vector<double>& s);

}  // namespace og::dist_detail
