// Geometric distribution: interval i has length d r^i, so s_i = d (r^i - 1) / (r - 1).
//
// Method: given the first interval d, the ratio r solves d (r^n - 1) / (r - 1) = 1; this is
// solved for log r with Brent's method on log((r^n - 1) / (r - 1)) evaluated with expm1/log1p,
// and positions are computed in log space, so no power overflows. Reference: e.g. J. F.
// Thompson, B. K. Soni, N. P. Weatherill, Handbook of Grid Generation (1999), ch. 32.
#include <cmath>

#include "families.hpp"
#include "opengrid/core/util/roots.hpp"

namespace og::dist_detail {
namespace {

constexpr double kMaxLogRatio = 60.0;

// log((r^i - 1) / (r - 1)) with r = exp(x), i >= 1.
double log_geometric_sum(double x, int i) {
  const auto di = static_cast<double>(i);
  if (std::fabs(x) < 1e-9) {
    return std::log(di) + 0.5 * (di - 1.0) * x;
  }
  if (x > 0.0) {
    return di * x + std::log(-std::expm1(-di * x)) - std::log(std::expm1(x));
  }
  return std::log(-std::expm1(di * x)) - std::log(-std::expm1(x));
}

bool fill(double log_d, double x, int n, std::vector<double>& s) {
  s.assign(static_cast<std::size_t>(n) + 1, 0.0);
  for (int i = 1; i < n; ++i) {
    s[static_cast<std::size_t>(i)] = std::exp(log_d + log_geometric_sum(x, i));
  }
  s.back() = 1.0;
  for (std::size_t i = 1; i < s.size(); ++i) {
    if (!(s[i] > s[i - 1]) || !std::isfinite(s[i])) {
      return false;
    }
  }
  return true;
}

}  // namespace

bool geometric_from_begin(double d, int n, std::vector<double>& s) {
  const double target = -std::log(d);
  auto f = [target, n](double x) { return log_geometric_sum(x, n) - target; };
  if (f(-kMaxLogRatio) > 0.0 || f(kMaxLogRatio) < 0.0) {
    return false;
  }
  const auto x = brent_root(f, -kMaxLogRatio, kMaxLogRatio);
  if (!x) {
    return false;
  }
  return fill(std::log(d), *x, n, s);
}

bool geometric_from_ratio(double ratio, int n, std::vector<double>& s) {
  const double x = std::log(ratio);
  if (std::fabs(x) > kMaxLogRatio) {
    return false;
  }
  const double log_d = -log_geometric_sum(x, n);
  return fill(log_d, x, n, s);
}

}  // namespace og::dist_detail
