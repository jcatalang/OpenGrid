#include "opengrid/core/geom/curve.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "opengrid/core/util/gauss_kronrod.hpp"

namespace og {
namespace {

constexpr int kProjectionSamplesPerSegment = 32;
constexpr std::size_t kMaxProjectionCandidates = 8;
constexpr double kGolden = 0.3819660112501051;  // (3 - sqrt(5)) / 2

}  // namespace

std::size_t find_segment(const std::vector<double>& knots, double t) {
  const std::size_t last = knots.size() - 2;
  if (t <= knots.front()) {
    return 0;
  }
  if (t >= knots.back()) {
    return last;
  }
  const auto it = std::upper_bound(knots.begin(), knots.end(), t);
  const auto index = static_cast<std::size_t>(it - knots.begin());
  if (index == 0) {
    return 0;
  }
  return std::min(index - 1, last);
}

std::vector<double> Curve::breakpoints() const {
  return {0.0, 1.0};
}

double Curve::speed(double t) const {
  return norm(d1(t));
}

void Curve::build_length_table() const {
  std::call_once(length_once_, [this]() {
    segment_params_ = breakpoints();
    cumulative_.assign(segment_params_.size(), 0.0);
    for (std::size_t i = 0; i + 1 < segment_params_.size(); ++i) {
      const double piece = integrate_gk15([this](double t) { return speed(t); }, segment_params_[i],
                                          segment_params_[i + 1]);
      cumulative_[i + 1] = cumulative_[i] + piece;
    }
  });
}

double Curve::length() const {
  build_length_table();
  return cumulative_.back();
}

double Curve::length_between(double t0, double t1) const {
  if (t0 == t1) {
    return 0.0;
  }
  if (t1 < t0) {
    return -length_between(t1, t0);
  }
  build_length_table();
  // Integrate segment by segment so derivative jumps at breakpoints are never straddled.
  double total = 0.0;
  const std::size_t first = find_segment(segment_params_, t0);
  const std::size_t last = find_segment(segment_params_, t1);
  for (std::size_t i = first; i <= last; ++i) {
    double a = segment_params_[i];
    double b = segment_params_[i + 1];
    if (i == first) {
      a = t0;
    }
    if (i == last) {
      b = t1;
    }
    if (a == segment_params_[i] && b == segment_params_[i + 1]) {
      total += cumulative_[i + 1] - cumulative_[i];
    } else {
      total += integrate_gk15([this](double t) { return speed(t); }, a, b);
    }
  }
  return total;
}

double Curve::t_at_length(double s) const {
  build_length_table();
  const double total = cumulative_.back();
  if (s <= 0.0) {
    return 0.0;
  }
  if (s >= total) {
    return 1.0;
  }
  const std::size_t seg = find_segment(cumulative_, s);
  double lo = segment_params_[seg];
  double hi = segment_params_[seg + 1];
  const double seg_length = cumulative_[seg + 1] - cumulative_[seg];
  const double target = s - cumulative_[seg];
  if (seg_length <= 0.0) {
    return lo;
  }

  // Safeguarded Newton on g(t) = L(lo_seg, t) - target with g'(t) = |x'(t)|; g is tracked
  // incrementally by integrating between successive iterates.
  const double seg_start = lo;
  double t = lo + (hi - lo) * (target / seg_length);
  double g = integrate_gk15([this](double u) { return speed(u); }, seg_start, t) - target;
  const double tolerance = 1e-14 * total;
  for (int iteration = 0; iteration < 100; ++iteration) {
    if (std::fabs(g) <= tolerance) {
      break;
    }
    if (g < 0.0) {
      lo = t;
    } else {
      hi = t;
    }
    const double v = speed(t);
    double t_new = 0.5 * (lo + hi);
    if (v > 0.0) {
      const double newton = t - g / v;
      if (newton > lo && newton < hi) {
        t_new = newton;
      }
    }
    if (t_new == t) {
      break;
    }
    g += integrate_gk15([this](double u) { return speed(u); }, t, t_new);
    t = t_new;
    if (hi - lo <= 4.0 * std::numeric_limits<double>::epsilon()) {
      break;
    }
  }
  return t;
}

double Curve::polish_projection(const Vec3& p, double t, double lo, double hi) const {
  // Newton iterations on the stationarity condition f(t) = (x(t) - p) . x'(t) = 0. Steps are
  // accepted while |f| decreases (the distance itself is flat to round-off at the minimum) and
  // only where d^2 is locally convex (f' > 0), so a maximum is never approached.
  auto gradient = [this, &p](double u) { return dot(eval(u) - p, d1(u)); };
  double f = gradient(t);
  for (int iteration = 0; iteration < 30; ++iteration) {
    const Vec3 r = eval(t) - p;
    const Vec3 x1 = d1(t);
    const double fp = dot(x1, x1) + dot(r, d2(t));
    if (fp <= 0.0 || f == 0.0) {
      break;
    }
    const double t_new = std::clamp(t - f / fp, lo, hi);
    if (t_new == t) {
      break;
    }
    const double f_new = gradient(t_new);
    if (std::fabs(f_new) >= std::fabs(f)) {
      break;
    }
    t = t_new;
    f = f_new;
  }
  return t;
}

CurveProjection Curve::project(const Vec3& p) const {
  // 1. Sample every segment uniformly in t.
  const std::vector<double> bps = breakpoints();
  std::vector<double> ts;
  for (std::size_t i = 0; i + 1 < bps.size(); ++i) {
    for (int k = 0; k < kProjectionSamplesPerSegment; ++k) {
      const double frac = static_cast<double>(k) / kProjectionSamplesPerSegment;
      ts.push_back(bps[i] + frac * (bps[i + 1] - bps[i]));
    }
  }
  ts.push_back(1.0);
  std::vector<double> d2s(ts.size());
  for (std::size_t i = 0; i < ts.size(); ++i) {
    const Vec3 r = eval(ts[i]) - p;
    d2s[i] = dot(r, r);
  }

  // 2. Candidates: local minima of the sampled squared distance, best first.
  std::vector<std::size_t> candidates;
  for (std::size_t i = 0; i < ts.size(); ++i) {
    bool left_ok = true;
    bool right_ok = true;
    if (i > 0) {
      left_ok = d2s[i] <= d2s[i - 1];
    }
    if (i + 1 < ts.size()) {
      right_ok = d2s[i] <= d2s[i + 1];
    }
    if (left_ok && right_ok) {
      candidates.push_back(i);
    }
  }
  std::stable_sort(candidates.begin(), candidates.end(),
                   [&d2s](std::size_t a, std::size_t b) { return d2s[a] < d2s[b]; });
  if (candidates.size() > kMaxProjectionCandidates) {
    candidates.resize(kMaxProjectionCandidates);
  }

  // 3. Refine each candidate on its bracket: golden-section search, then Newton polish.
  CurveProjection best;
  best.distance = std::numeric_limits<double>::infinity();
  for (const std::size_t i : candidates) {
    double lo = ts[i];
    double hi = ts[i];
    if (i > 0) {
      lo = ts[i - 1];
    }
    if (i + 1 < ts.size()) {
      hi = ts[i + 1];
    }
    auto dist2 = [this, &p](double t) {
      const Vec3 r = eval(t) - p;
      return dot(r, r);
    };
    double a = lo;
    double b = hi;
    double x1 = a + kGolden * (b - a);
    double x2 = b - kGolden * (b - a);
    double f1 = dist2(x1);
    double f2 = dist2(x2);
    for (int iteration = 0; iteration < 200 && (b - a) > 1e-12 * (hi - lo + 1e-300); ++iteration) {
      if (f1 <= f2) {
        b = x2;
        x2 = x1;
        f2 = f1;
        x1 = a + kGolden * (b - a);
        f1 = dist2(x1);
      } else {
        a = x1;
        x1 = x2;
        f1 = f2;
        x2 = b - kGolden * (b - a);
        f2 = dist2(x2);
      }
    }
    double t = 0.5 * (a + b);
    // The sampled endpoints of the bracket may be better (minimum at a bracket end).
    if (dist2(lo) < dist2(t)) {
      t = lo;
    }
    if (dist2(hi) < dist2(t)) {
      t = hi;
    }
    t = polish_projection(p, t, lo, hi);
    const Vec3 point = eval(t);
    const double d = distance(point, p);
    if (d < best.distance) {
      best.t = t;
      best.point = point;
      best.distance = d;
    }
  }
  return best;
}

std::vector<double> Curve::tessellate(double rel_tol) const {
  const double tolerance = rel_tol * length();
  constexpr int kMinDepth = 2;
  constexpr int kMaxDepth = 14;
  std::vector<double> params;
  // Recursive bisection of [a, b] until the curve midpoint is within tolerance of the chord.
  auto refine = [&](auto&& self, double a, const Vec3& xa, double b, const Vec3& xb,
                    int depth) -> void {
    const double m = 0.5 * (a + b);
    const Vec3 xm = eval(m);
    const double deviation = distance(xm, 0.5 * (xa + xb));
    if (depth >= kMaxDepth || (depth >= kMinDepth && deviation <= tolerance)) {
      params.push_back(b);
      return;
    }
    self(self, a, xa, m, xm, depth + 1);
    self(self, m, xm, b, xb, depth + 1);
  };
  const std::vector<double> bps = breakpoints();
  params.push_back(bps.front());
  for (std::size_t i = 0; i + 1 < bps.size(); ++i) {
    refine(refine, bps[i], eval(bps[i]), bps[i + 1], eval(bps[i + 1]), 0);
  }
  return params;
}

}  // namespace og
