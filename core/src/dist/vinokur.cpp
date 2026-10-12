// Tanh and sine stretching functions.
//
// Method: M. Vinokur, "On one-dimensional stretching functions for finite-difference
// calculations", J. Comput. Phys. 50 (1983) 215-234. With xi in [0, 1] and prescribed end
// derivatives S0 = s'(0), S1 = s'(1):
//   two-sided:  A = sqrt(S1 / S0), B = 1 / sqrt(S0 S1), s = u / (A + (1 - A) u), where u is a
//               symmetric stretching function with u'(0) = u'(1) = 1 / B;
//   one-sided:  B = 1 / S0 and s is a function with s'(0) = 1 / B.
// Tanh family (Vinokur's functions):
//   two-sided  B > 1: u = 1/2 (1 + tanh(D (xi - 1/2)) / tanh(D / 2)),   sinh(D) / D = B
//              B < 1: u = 1/2 (1 + tan(D (xi - 1/2)) / tan(D / 2)),     sin(D) / D = B
//   one-sided  B > 1: s = 1 + tanh(D (xi - 1) / 2) / tanh(D / 2),        sinh(D) / D = B
//              B < 1: s = 1 + tan(D (xi - 1) / 2) / tan(D / 2),          sin(D) / D = B
// Sine family:
//   one-sided  B > 1: s = sinh(D xi) / sinh(D),  sinh(D) / D = B
//              B < 1: s = sin(D xi) / sin(D),    sin(D) / D = B, monotone only for D <= pi/2
//                     (B >= 2 / pi); otherwise the caller falls back to the tanh family
//   two-sided  B > 1: u = 1/2 (1 + sin(D (xi - 1/2)) / sin(D / 2)),   (D/2) cot(D/2) = 1 / B
//              B < 1: u = 1/2 (1 + sinh(D (xi - 1/2)) / sinh(D / 2)), (D/2) coth(D/2) = 1 / B
// All formulas are evaluated in cancellation-free product form (sum-to-product identities)
// and hyperbolic ones in log space, so tiny spacings keep full relative accuracy and nothing
// overflows. |B - 1| < 1e-13 is treated as uniform.
#include <cmath>
#include <numbers>

#include "families.hpp"
#include "opengrid/core/util/roots.hpp"

namespace og::dist_detail {
namespace {

constexpr double kUniformBand = 1e-13;
constexpr double kPi = std::numbers::pi;

// log(sinh(D) / D) for D > 0.
double log_sinhc(double d) {
  if (d < 1e-4) {
    const double d2 = d * d;
    return std::log1p(d2 / 6.0 + d2 * d2 / 120.0);
  }
  return log_sinh(d) - std::log(d);
}

// D > 0 with sinh(D) / D = b, b > 1.
double solve_sinhc(double b) {
  const double target = std::log(b);
  double hi = 2.0 * (target + 2.0);
  auto f = [target](double d) { return log_sinhc(d) - target; };
  while (f(hi) < 0.0) {
    hi *= 2.0;
  }
  const auto root = brent_root(f, 1e-12, hi);
  if (root) {
    return *root;
  }
  return hi;
}

// D in (0, pi) with sin(D) / D = b, 0 < b < 1.
double solve_sinc(double b) {
  auto f = [b](double d) { return std::sin(d) / d - b; };
  const auto root = brent_root(f, 1e-12, kPi);
  if (root) {
    return *root;
  }
  return kPi;
}

// D in (0, pi) with (D/2) cot(D/2) = target, 0 < target < 1.
double solve_half_cot(double target) {
  auto f = [target](double d) { return 0.5 * d / std::tan(0.5 * d) - target; };
  const auto root = brent_root(f, 1e-12, kPi * (1.0 - 1e-16));
  if (root) {
    return *root;
  }
  return kPi;
}

// D > 0 with (D/2) coth(D/2) = target, target > 1.
double solve_half_coth(double target) {
  auto f = [target](double d) { return 0.5 * d / std::tanh(0.5 * d) - target; };
  double hi = 2.0 * target + 2.0;
  while (f(hi) < 0.0) {
    hi *= 2.0;
  }
  const auto root = brent_root(f, 1e-12, hi);
  if (root) {
    return *root;
  }
  return hi;
}

double xi_at(int i, int n) {
  return static_cast<double>(i) / static_cast<double>(n);
}

std::size_t index(int i) {
  return static_cast<std::size_t>(i);
}

// Vinokur's asymmetry transform s = u / (A + (1 - A) u) applied to a symmetric u
// (u(1 - xi) = 1 - u(xi)). Points past the middle are computed from the complement
// 1 - s = A u(1 - xi) / (A + (1 - A) u(xi)), so spacings near xi = 1 keep full relative
// accuracy; the last interval is returned separately for the same reason.
void asymmetric(double a, std::vector<double>& u_and_s, double& last_spacing) {
  const std::vector<double> u = u_and_s;
  const std::size_t n = u.size() - 1;
  for (std::size_t i = 0; i <= n; ++i) {
    const double den = a + (1.0 - a) * u[i];
    if (2 * i <= n) {
      u_and_s[i] = u[i] / den;
    } else {
      u_and_s[i] = 1.0 - a * u[n - i] / den;
    }
  }
  last_spacing = a * u[1] / (a + (1.0 - a) * u[n - 1]);
}

void finish(std::vector<double>& s) {
  s.front() = 0.0;
  s.back() = 1.0;
}

}  // namespace

bool tanh_one_sided(double s0_derivative, int n, std::vector<double>& s) {
  s.assign(index(n) + 1, 0.0);
  const double b = 1.0 / s0_derivative;
  for (int i = 0; i <= n; ++i) {
    s[index(i)] = xi_at(i, n);
  }
  if (b > 1.0 + kUniformBand) {
    const double d = solve_sinhc(b);
    const double log_den = log_sinh(0.5 * d);
    for (int i = 1; i < n; ++i) {
      const double xi = xi_at(i, n);
      // 1 + tanh(D(xi-1)/2)/tanh(D/2) = sinh(D xi/2) / (sinh(D/2) cosh(D(1-xi)/2))
      s[index(i)] = std::exp(log_sinh(0.5 * d * xi) - log_den - log_cosh(0.5 * d * (1.0 - xi)));
    }
  } else if (b < 1.0 - kUniformBand) {
    const double d = solve_sinc(b);
    const double den = std::sin(0.5 * d);
    for (int i = 1; i < n; ++i) {
      const double xi = xi_at(i, n);
      s[index(i)] = std::sin(0.5 * d * xi) / (den * std::cos(0.5 * d * (1.0 - xi)));
    }
  }
  finish(s);
  return true;
}

bool sine_one_sided(double s0_derivative, int n, std::vector<double>& s) {
  s.assign(index(n) + 1, 0.0);
  const double b = 1.0 / s0_derivative;
  for (int i = 0; i <= n; ++i) {
    s[index(i)] = xi_at(i, n);
  }
  if (b > 1.0 + kUniformBand) {
    const double d = solve_sinhc(b);
    const double log_den = log_sinh(d);
    for (int i = 1; i < n; ++i) {
      s[index(i)] = std::exp(log_sinh(d * xi_at(i, n)) - log_den);
    }
  } else if (b < 1.0 - kUniformBand) {
    if (b < 2.0 / kPi) {
      return false;  // sin(D xi) / sin(D) is not monotone for D > pi / 2
    }
    const double d = solve_sinc(b);
    const double den = std::sin(d);
    for (int i = 1; i < n; ++i) {
      s[index(i)] = std::sin(d * xi_at(i, n)) / den;
    }
  }
  finish(s);
  return true;
}

bool tanh_two_sided(double s0_derivative, double s1_derivative, int n, std::vector<double>& s,
                    double& last_spacing) {
  s.assign(index(n) + 1, 0.0);
  const double a = std::sqrt(s1_derivative / s0_derivative);
  const double b = 1.0 / std::sqrt(s0_derivative * s1_derivative);
  for (int i = 0; i <= n; ++i) {
    s[index(i)] = xi_at(i, n);
  }
  if (b > 1.0 + kUniformBand) {
    const double d = solve_sinhc(b);
    const double log_den = std::numbers::ln2 + log_sinh(0.5 * d);
    for (int i = 1; i < n; ++i) {
      const double xi = xi_at(i, n);
      // 1/2 (1 + tanh(D(xi-1/2))/tanh(D/2)) = sinh(D xi) / (2 sinh(D/2) cosh(D(xi-1/2)))
      s[index(i)] = std::exp(log_sinh(d * xi) - log_den - log_cosh(d * (xi - 0.5)));
    }
  } else if (b < 1.0 - kUniformBand) {
    const double d = solve_sinc(b);
    const double den = 2.0 * std::sin(0.5 * d);
    for (int i = 1; i < n; ++i) {
      const double xi = xi_at(i, n);
      s[index(i)] = std::sin(d * xi) / (den * std::cos(d * (xi - 0.5)));
    }
  }
  finish(s);
  asymmetric(a, s, last_spacing);
  finish(s);
  return true;
}

bool sine_two_sided(double s0_derivative, double s1_derivative, int n, std::vector<double>& s,
                    double& last_spacing) {
  s.assign(index(n) + 1, 0.0);
  const double a = std::sqrt(s1_derivative / s0_derivative);
  const double b = 1.0 / std::sqrt(s0_derivative * s1_derivative);
  for (int i = 0; i <= n; ++i) {
    s[index(i)] = xi_at(i, n);
  }
  if (b > 1.0 + kUniformBand) {
    const double d = solve_half_cot(1.0 / b);
    const double den = std::sin(0.5 * d);
    for (int i = 1; i < n; ++i) {
      const double xi = xi_at(i, n);
      // 1/2 (1 + sin(D(xi-1/2))/sin(D/2)) = sin(D xi/2) cos(D(1-xi)/2) / sin(D/2)
      s[index(i)] = std::sin(0.5 * d * xi) * std::cos(0.5 * d * (1.0 - xi)) / den;
    }
  } else if (b < 1.0 - kUniformBand) {
    const double d = solve_half_coth(1.0 / b);
    const double log_den = log_sinh(0.5 * d);
    for (int i = 1; i < n; ++i) {
      const double xi = xi_at(i, n);
      s[index(i)] = std::exp(log_sinh(0.5 * d * xi) + log_cosh(0.5 * d * (1.0 - xi)) - log_den);
    }
  }
  finish(s);
  asymmetric(a, s, last_spacing);
  finish(s);
  return true;
}

}  // namespace og::dist_detail
