// Distribution dispatcher: input validation, discrete end-spacing matching for the Vinokur
// families, mirroring for end-only requests, and fallbacks with warnings.
//
// Matching: Vinokur's functions prescribe derivatives s'(0), s'(1); the discrete spacings
// s_1 - s_0 and s_n - s_{n-1} differ by O(1/n) for strong stretching. One-sided requests solve
// for log s'(0) with Brent's method; two-sided requests solve the 2x2 system for
// (log s'(0), log s'(1)) with a damped Newton iteration (finite-difference Jacobian).
#include "opengrid/core/dist/distribution.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#include "families.hpp"
#include "opengrid/core/util/roots.hpp"

namespace og {
namespace {

using dist_detail::geometric_from_begin;
using dist_detail::geometric_from_ratio;

constexpr double kMatchTolerance = 1e-13;  // relative, on log spacing
constexpr double kAccept = 1e-10;          // larger residuals are reported as warnings

std::string fmt(double value) {
  std::array<char, 32> buffer{};
  std::snprintf(buffer.data(), buffer.size(), "%.6g", value);
  return buffer.data();
}

std::vector<double> uniform(int n) {
  std::vector<double> s(static_cast<std::size_t>(n) + 1);
  for (int i = 0; i <= n; ++i) {
    s[static_cast<std::size_t>(i)] = static_cast<double>(i) / static_cast<double>(n);
  }
  s.back() = 1.0;
  return s;
}

std::vector<double> mirror(const std::vector<double>& s) {
  std::vector<double> out(s.size());
  for (std::size_t i = 0; i < s.size(); ++i) {
    out[i] = 1.0 - s[s.size() - 1 - i];
  }
  out.front() = 0.0;
  out.back() = 1.0;
  return out;
}

bool valid(const std::vector<double>& s) {
  if (s.size() < 2 || s.front() != 0.0 || s.back() != 1.0) {
    return false;
  }
  for (std::size_t i = 1; i < s.size(); ++i) {
    if (!std::isfinite(s[i]) || !(s[i] > s[i - 1])) {
      return false;
    }
  }
  return true;
}

void check_positive(const std::optional<double>& value, const char* what) {
  if (value && !(std::isfinite(*value) && *value > 0.0)) {
    throw std::invalid_argument(std::string("distribution: ") + what +
                                " must be positive and finite");
  }
}

using OneSided = bool (*)(double, int, std::vector<double>&);
using TwoSided = bool (*)(double, double, int, std::vector<double>&, double&);

// One-sided Vinokur family with discrete first spacing d. Returns false if the family cannot
// represent the request (the caller then falls back).
bool match_one_sided(OneSided family, double d, int n, std::vector<double>& s) {
  std::vector<double> work;
  auto residual = [&](double x) {
    if (!family(std::exp(x), n, work)) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return std::log(work[1]) - std::log(d);
  };
  const double x0 = std::log(d * n);
  double lo = x0 - 1.0;
  double hi = x0 + 1.0;
  for (int k = 0; k < 80 && !(residual(lo) <= 0.0); ++k) {
    lo -= std::fabs(lo - x0) + 1.0;
  }
  for (int k = 0; k < 80 && !(residual(hi) >= 0.0); ++k) {
    hi += std::fabs(hi - x0) + 1.0;
    if (std::isnan(residual(hi))) {
      return false;
    }
  }
  const auto x = brent_root(residual, lo, hi);
  if (!x || !family(std::exp(*x), n, s)) {
    return false;
  }
  return valid(s);
}

// Two-sided Vinokur family; returns the achieved maximum |log(achieved / requested)|.
double match_two_sided(TwoSided family, double d0, double d1, int n, std::vector<double>& s) {
  std::vector<double> work;
  double last_spacing = 0.0;
  auto residual = [&](const std::array<double, 2>& x, std::array<double, 2>& r) {
    family(std::exp(x[0]), std::exp(x[1]), n, work, last_spacing);
    r[0] = std::log(work[1] - work[0]) - std::log(d0);
    r[1] = std::log(last_spacing) - std::log(d1);
    return std::isfinite(r[0]) && std::isfinite(r[1]);
  };
  auto size = [](const std::array<double, 2>& r) {
    return std::max(std::fabs(r[0]), std::fabs(r[1]));
  };
  std::array<double, 2> x = {std::log(d0 * n), std::log(d1 * n)};
  std::array<double, 2> r{};
  if (!residual(x, r)) {
    s = uniform(n);
    return std::numeric_limits<double>::infinity();
  }
  for (int iteration = 0; iteration < 100 && size(r) > kMatchTolerance; ++iteration) {
    constexpr double kH = 1e-7;
    std::array<std::array<double, 2>, 2> jac{};
    for (std::size_t j = 0; j < 2; ++j) {
      std::array<double, 2> xp = x;
      xp[j] += kH;
      std::array<double, 2> rp{};
      residual(xp, rp);
      jac[0][j] = (rp[0] - r[0]) / kH;
      jac[1][j] = (rp[1] - r[1]) / kH;
    }
    const double det = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
    if (!(std::fabs(det) > 1e-300)) {
      break;
    }
    const std::array<double, 2> step = {-(jac[1][1] * r[0] - jac[0][1] * r[1]) / det,
                                        -(-jac[1][0] * r[0] + jac[0][0] * r[1]) / det};
    bool accepted = false;
    for (int halving = 0; halving < 20; ++halving) {
      const double t = std::ldexp(1.0, -halving);
      const std::array<double, 2> xn = {x[0] + t * step[0], x[1] + t * step[1]};
      std::array<double, 2> rn{};
      if (residual(xn, rn) && size(rn) < size(r)) {
        x = xn;
        r = rn;
        accepted = true;
        break;
      }
    }
    if (!accepted) {
      break;
    }
  }
  family(std::exp(x[0]), std::exp(x[1]), n, s, last_spacing);
  return size(r);
}

Distribution vinokur(DistributionType type, const DistributionSpec& spec, int n) {
  Distribution out;
  const char* name = to_string(type);
  OneSided one = dist_detail::tanh_one_sided;
  TwoSided two = dist_detail::tanh_two_sided;
  if (type == DistributionType::Sine) {
    one = dist_detail::sine_one_sided;
    two = dist_detail::sine_two_sided;
  }
  const bool has_begin = spec.begin_spacing.has_value();
  const bool has_end = spec.end_spacing.has_value();
  if (!has_begin && !has_end) {
    out.s = uniform(n);
    return out;
  }
  if (has_begin != has_end) {
    double d = 0.0;
    const char* side = "begin";
    if (has_begin) {
      d = *spec.begin_spacing;
    } else {
      d = *spec.end_spacing;
      side = "end";
    }
    if (n == 1 || d >= 1.0) {
      out.s = uniform(n);
      if (std::fabs(d * n - 1.0) > 1e-12) {
        out.warnings.push_back(std::string(name) + ": " + side + " spacing " + fmt(d) +
                               " cannot be achieved with " + std::to_string(n) +
                               " interval(s); using a uniform distribution");
      }
      return out;
    }
    std::vector<double> s;
    if (!match_one_sided(one, d, n, s)) {
      if (type == DistributionType::Sine && match_one_sided(dist_detail::tanh_one_sided, d, n, s)) {
        out.warnings.push_back(std::string("sine: ") + side + " spacing " + fmt(d) +
                               " is too large for the sine family; used tanh instead");
      } else {
        s = uniform(n);
        out.warnings.push_back(std::string(name) + ": could not achieve " + side + " spacing " +
                               fmt(d) + "; using a uniform distribution");
      }
    }
    if (has_begin) {
      out.s = std::move(s);
    } else {
      out.s = mirror(s);
    }
    return out;
  }

  const double d0 = *spec.begin_spacing;
  const double d1 = *spec.end_spacing;
  if (n == 1) {
    out.s = uniform(1);
    out.warnings.push_back(std::string(name) + ": a single interval cannot have spacings " +
                           fmt(d0) + " and " + fmt(d1) + "; using one interval");
    return out;
  }
  if (n == 2) {
    out.s = {0.0, d0 / (d0 + d1), 1.0};
    if (std::fabs(d0 + d1 - 1.0) > 1e-12) {
      out.warnings.push_back(std::string(name) + ": two intervals cannot have spacings " + fmt(d0) +
                             " and " + fmt(d1) + "; split proportionally");
    }
    return out;
  }
  if (d0 + d1 >= 1.0) {
    out.s = uniform(n);
    out.warnings.push_back(std::string(name) + ": begin + end spacing (" + fmt(d0 + d1) +
                           ") must be below 1; using a uniform distribution");
    return out;
  }
  std::vector<double> s;
  double error = match_two_sided(two, d0, d1, n, s);
  if (type == DistributionType::Sine && (error > kAccept || !valid(s))) {
    // The sine family cannot cluster beyond about (pi / 2n)^2 at an end.
    std::vector<double> t;
    const double tanh_error = match_two_sided(dist_detail::tanh_two_sided, d0, d1, n, t);
    if (tanh_error <= kAccept && valid(t)) {
      s = std::move(t);
      error = tanh_error;
      out.warnings.push_back("sine: spacings " + fmt(d0) + " / " + fmt(d1) +
                             " are beyond the sine family with " + std::to_string(n + 1) +
                             " points; used tanh instead");
    }
  }
  if (!valid(s)) {
    out.s = uniform(n);
    out.warnings.push_back(std::string(name) + ": could not achieve spacings " + fmt(d0) + " and " +
                           fmt(d1) + "; using a uniform distribution");
    return out;
  }
  out.s = s;
  if (error > kAccept) {
    const double a0 = s[1] - s[0];
    const double a1 = s[s.size() - 1] - s[s.size() - 2];
    out.warnings.push_back(std::string(name) + ": requested spacings " + fmt(d0) + " / " + fmt(d1) +
                           " not achievable with " + std::to_string(n + 1) + " points; achieved " +
                           fmt(a0) + " / " + fmt(a1));
  }
  return out;
}

Distribution geometric(const DistributionSpec& spec, int n) {
  Distribution out;
  std::vector<double> s;
  bool ok = true;
  std::string what;
  if (spec.growth_ratio) {
    ok = geometric_from_ratio(*spec.growth_ratio, n, s);
    what = "growth ratio " + fmt(*spec.growth_ratio);
    out.s = std::move(s);
  } else {
    double d = 0.0;
    if (spec.begin_spacing) {
      d = *spec.begin_spacing;
      what = "begin spacing " + fmt(d);
    } else if (spec.end_spacing) {
      d = *spec.end_spacing;
      what = "end spacing " + fmt(d);
    }
    if (n == 1 || d >= 1.0) {
      ok = n == 1 && std::fabs(d - 1.0) <= 1e-12;
      s = uniform(n);
    } else if (std::fabs(d * n - 1.0) <= 1e-15) {
      s = uniform(n);
    } else {
      ok = geometric_from_begin(d, n, s);
    }
    if (spec.end_spacing && ok) {
      out.s = mirror(s);
    } else {
      out.s = std::move(s);
    }
  }
  if (!ok || !valid(out.s)) {
    out.s = uniform(n);
    out.warnings.push_back("geometric: " + what + " cannot be achieved with " + std::to_string(n) +
                           " interval(s); using a uniform distribution");
  }
  return out;
}

}  // namespace

const char* to_string(DistributionType type) {
  switch (type) {
    case DistributionType::Uniform:
      return "uniform";
    case DistributionType::Geometric:
      return "geometric";
    case DistributionType::Tanh:
      return "tanh";
    case DistributionType::Sine:
      return "sine";
    case DistributionType::Explicit:
      return "explicit";
  }
  return "unknown";
}

Distribution distribute(int n_points, const DistributionSpec& spec) {
  if (n_points < 2) {
    throw std::invalid_argument("distribution: at least 2 points are required");
  }
  check_positive(spec.begin_spacing, "begin spacing");
  check_positive(spec.end_spacing, "end spacing");
  check_positive(spec.growth_ratio, "growth ratio");
  const int n = n_points - 1;
  Distribution out;
  switch (spec.type) {
    case DistributionType::Uniform:
      if (spec.begin_spacing || spec.end_spacing || spec.growth_ratio) {
        throw std::invalid_argument("distribution: uniform takes no spacing or growth ratio");
      }
      out.s = uniform(n);
      break;
    case DistributionType::Geometric: {
      const int given = static_cast<int>(spec.begin_spacing.has_value()) +
                        static_cast<int>(spec.end_spacing.has_value()) +
                        static_cast<int>(spec.growth_ratio.has_value());
      if (given != 1) {
        throw std::invalid_argument(
            "distribution: geometric needs exactly one of begin spacing, end spacing and "
            "growth ratio");
      }
      out = geometric(spec, n);
      break;
    }
    case DistributionType::Tanh:
    case DistributionType::Sine:
      if (spec.growth_ratio) {
        throw std::invalid_argument("distribution: growth ratio applies to geometric only");
      }
      out = vinokur(spec.type, spec, n);
      break;
    case DistributionType::Explicit:
      if (spec.begin_spacing || spec.end_spacing || spec.growth_ratio) {
        throw std::invalid_argument("distribution: explicit takes no spacing or growth ratio");
      }
      if (spec.positions.size() != static_cast<std::size_t>(n_points) ||
          spec.positions.front() != 0.0 || spec.positions.back() != 1.0 || !valid(spec.positions)) {
        throw std::invalid_argument(
            "distribution: explicit positions must be n_points strictly increasing values from 0 "
            "to 1");
      }
      out.s = spec.positions;
      break;
  }
  if (!valid(out.s)) {
    out.s = uniform(n);
    out.warnings.emplace_back("internal fallback to a uniform distribution");
  }
  return out;
}

double max_stretching_ratio(const std::vector<double>& positions) {
  double worst = 1.0;
  for (std::size_t i = 2; i < positions.size(); ++i) {
    const double a = positions[i - 1] - positions[i - 2];
    const double b = positions[i] - positions[i - 1];
    if (a > 0.0 && b > 0.0) {
      worst = std::max({worst, a / b, b / a});
    }
  }
  return worst;
}

}  // namespace og
