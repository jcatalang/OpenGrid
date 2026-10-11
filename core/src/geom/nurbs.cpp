#include "opengrid/core/geom/nurbs.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>

namespace og {
namespace {

std::size_t idx(int i) {
  return static_cast<std::size_t>(i);
}

}  // namespace

Nurbs::Nurbs(int degree, std::vector<double> knots, std::vector<Vec3> control_points,
             std::vector<double> weights)
    : degree_(degree),
      knots_(std::move(knots)),
      control_points_(std::move(control_points)),
      weights_(std::move(weights)) {
  if (degree_ < 1) {
    throw std::invalid_argument("nurbs: degree must be >= 1");
  }
  const std::size_t n_ctrl = control_points_.size();
  if (n_ctrl < idx(degree_) + 1) {
    throw std::invalid_argument("nurbs: need at least degree + 1 control points");
  }
  if (weights_.size() != n_ctrl) {
    throw std::invalid_argument("nurbs: " + std::to_string(weights_.size()) + " weights for " +
                                std::to_string(n_ctrl) + " control points");
  }
  if (knots_.size() != n_ctrl + idx(degree_) + 1) {
    throw std::invalid_argument("nurbs: expected " + std::to_string(n_ctrl + idx(degree_) + 1) +
                                " knots, got " + std::to_string(knots_.size()));
  }
  for (std::size_t i = 0; i < n_ctrl; ++i) {
    if (!(weights_[i] > 0.0)) {
      throw std::invalid_argument("nurbs: weights[" + std::to_string(i) + "] must be positive");
    }
  }
  for (std::size_t i = 1; i < knots_.size(); ++i) {
    if (knots_[i] < knots_[i - 1]) {
      throw std::invalid_argument("nurbs: knots[" + std::to_string(i) + "] decreases");
    }
  }
  u_lo_ = knots_[idx(degree_)];
  u_hi_ = knots_[n_ctrl];
  if (!(u_hi_ > u_lo_)) {
    throw std::invalid_argument("nurbs: empty parameter domain");
  }
  // Interior multiplicity must not exceed the degree (curve would be discontinuous).
  std::size_t run = 1;
  for (std::size_t i = 1; i < knots_.size(); ++i) {
    if (knots_[i] == knots_[i - 1]) {
      run += 1;
    } else {
      run = 1;
    }
    const bool interior = knots_[i] > u_lo_ && knots_[i] < u_hi_;
    if (interior && run > idx(degree_)) {
      throw std::invalid_argument("nurbs: interior knot " + std::to_string(knots_[i]) +
                                  " has multiplicity > degree");
    }
  }
}

int Nurbs::find_span(double u) const {
  const int n = static_cast<int>(control_points_.size()) - 1;
  if (u >= knots_[idx(n + 1)]) {
    // Last non-empty span.
    int span = n;
    while (span > degree_ && knots_[idx(span)] >= knots_[idx(span + 1)]) {
      span -= 1;
    }
    return span;
  }
  if (u <= knots_[idx(degree_)]) {
    int span = degree_;
    while (span < n && knots_[idx(span + 1)] <= u) {
      span += 1;
    }
    return span;
  }
  int low = degree_;
  int high = n + 1;
  int mid = (low + high) / 2;
  while (u < knots_[idx(mid)] || u >= knots_[idx(mid + 1)]) {
    if (u < knots_[idx(mid)]) {
      high = mid;
    } else {
      low = mid;
    }
    mid = (low + high) / 2;
  }
  return mid;
}

Nurbs::Derivatives Nurbs::evaluate(double t) const {
  const double scale = u_hi_ - u_lo_;
  const double u = u_lo_ + std::clamp(t, 0.0, 1.0) * scale;
  const int p = degree_;
  const int span = find_span(u);
  int n_der = 2;
  if (p < 2) {
    n_der = p;
  }

  // Algorithm A2.3: ders[k][j] = k-th derivative of N_{span-p+j, p}(u).
  const std::size_t pp = idx(p) + 1;
  std::vector<std::vector<double>> ndu(pp, std::vector<double>(pp, 0.0));
  std::vector<double> left(pp, 0.0);
  std::vector<double> right(pp, 0.0);
  ndu[0][0] = 1.0;
  for (int j = 1; j <= p; ++j) {
    left[idx(j)] = u - knots_[idx(span + 1 - j)];
    right[idx(j)] = knots_[idx(span + j)] - u;
    double saved = 0.0;
    for (int r = 0; r < j; ++r) {
      ndu[idx(j)][idx(r)] = right[idx(r + 1)] + left[idx(j - r)];
      const double temp = ndu[idx(r)][idx(j - 1)] / ndu[idx(j)][idx(r)];
      ndu[idx(r)][idx(j)] = saved + right[idx(r + 1)] * temp;
      saved = left[idx(j - r)] * temp;
    }
    ndu[idx(j)][idx(j)] = saved;
  }
  std::vector<std::vector<double>> ders(3, std::vector<double>(pp, 0.0));
  for (int j = 0; j <= p; ++j) {
    ders[0][idx(j)] = ndu[idx(j)][idx(p)];
  }
  std::vector<std::vector<double>> a(2, std::vector<double>(pp, 0.0));
  for (int r = 0; r <= p; ++r) {
    std::size_t s1 = 0;
    std::size_t s2 = 1;
    a[0][0] = 1.0;
    for (int k = 1; k <= n_der; ++k) {
      double d = 0.0;
      const int rk = r - k;
      const int pk = p - k;
      if (r >= k) {
        a[s2][0] = a[s1][0] / ndu[idx(pk + 1)][idx(rk)];
        d = a[s2][0] * ndu[idx(rk)][idx(pk)];
      }
      int j1 = 1;
      if (rk < -1) {
        j1 = -rk;
      }
      int j2 = p - r;
      if (r - 1 <= pk) {
        j2 = k - 1;
      }
      for (int j = j1; j <= j2; ++j) {
        a[s2][idx(j)] = (a[s1][idx(j)] - a[s1][idx(j - 1)]) / ndu[idx(pk + 1)][idx(rk + j)];
        d += a[s2][idx(j)] * ndu[idx(rk + j)][idx(pk)];
      }
      if (r <= pk) {
        a[s2][idx(k)] = -a[s1][idx(k - 1)] / ndu[idx(pk + 1)][idx(r)];
        d += a[s2][idx(k)] * ndu[idx(r)][idx(pk)];
      }
      ders[idx(k)][idx(r)] = d;
      std::swap(s1, s2);
    }
  }
  double factor = p;
  for (int k = 1; k <= n_der; ++k) {
    for (int j = 0; j <= p; ++j) {
      ders[idx(k)][idx(j)] *= factor;
    }
    factor *= p - k;
  }

  // Homogeneous derivatives A^(k) = sum N^(k) w P, w^(k) = sum N^(k) w.
  std::array<Vec3, 3> aw{};
  std::array<double, 3> w{};
  for (std::size_t k = 0; k < 3; ++k) {
    for (int j = 0; j <= p; ++j) {
      const std::size_t i = idx(span - p + j);
      const double nw = ders[k][idx(j)] * weights_[i];
      aw[k] += nw * control_points_[i];
      w[k] += nw;
    }
  }
  // Rational derivatives (A4.2), then chain rule du/dt = scale.
  Derivatives out;
  out.c0 = aw[0] / w[0];
  out.c1 = (aw[1] - w[1] * out.c0) / w[0];
  out.c2 = (aw[2] - 2.0 * w[1] * out.c1 - w[2] * out.c0) / w[0];
  out.c1 *= scale;
  out.c2 *= scale * scale;
  return out;
}

Vec3 Nurbs::eval(double t) const {
  return evaluate(t).c0;
}

Vec3 Nurbs::d1(double t) const {
  return evaluate(t).c1;
}

Vec3 Nurbs::d2(double t) const {
  return evaluate(t).c2;
}

std::vector<double> Nurbs::breakpoints() const {
  std::vector<double> bps = {0.0};
  for (const double k : knots_) {
    if (k > u_lo_ && k < u_hi_) {
      const double t = (k - u_lo_) / (u_hi_ - u_lo_);
      if (t > bps.back()) {
        bps.push_back(t);
      }
    }
  }
  bps.push_back(1.0);
  return bps;
}

}  // namespace og
