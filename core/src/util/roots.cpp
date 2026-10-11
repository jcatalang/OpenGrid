#include "opengrid/core/util/roots.hpp"

#include <cmath>
#include <limits>

namespace og {

std::optional<double> brent_root(const std::function<double(double)>& f, double a, double b,
                                 double x_tol, int max_iterations) {
  double fa = f(a);
  double fb = f(b);
  if (fa == 0.0) {
    return a;
  }
  if (fb == 0.0) {
    return b;
  }
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0) || std::isnan(fa) || std::isnan(fb)) {
    return std::nullopt;
  }
  constexpr double kEps = std::numeric_limits<double>::epsilon();
  double c = a;
  double fc = fa;
  double d = b - a;
  double e = d;
  for (int iteration = 0; iteration < max_iterations; ++iteration) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c = a;
      fc = fa;
      d = b - a;
      e = d;
    }
    if (std::fabs(fc) < std::fabs(fb)) {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    const double tol1 = 2.0 * kEps * std::fabs(b) + 0.5 * x_tol;
    const double xm = 0.5 * (c - b);
    if (std::fabs(xm) <= tol1 || fb == 0.0) {
      return b;
    }
    if (std::fabs(e) >= tol1 && std::fabs(fa) > std::fabs(fb)) {
      const double s = fb / fa;
      double p = 0.0;
      double q = 0.0;
      if (a == c) {
        p = 2.0 * xm * s;
        q = 1.0 - s;
      } else {
        const double qq = fa / fc;
        const double r = fb / fc;
        p = s * (2.0 * xm * qq * (qq - r) - (b - a) * (r - 1.0));
        q = (qq - 1.0) * (r - 1.0) * (s - 1.0);
      }
      if (p > 0.0) {
        q = -q;
      }
      p = std::fabs(p);
      const double min1 = 3.0 * xm * q - std::fabs(tol1 * q);
      const double min2 = std::fabs(e * q);
      double limit = min2;
      if (min1 < min2) {
        limit = min1;
      }
      if (2.0 * p < limit) {
        e = d;
        d = p / q;
      } else {
        d = xm;
        e = d;
      }
    } else {
      d = xm;
      e = d;
    }
    a = b;
    fa = fb;
    if (std::fabs(d) > tol1) {
      b += d;
    } else {
      b += std::copysign(tol1, xm);
    }
    fb = f(b);
    if (std::isnan(fb)) {
      return std::nullopt;
    }
  }
  return b;
}

}  // namespace og
