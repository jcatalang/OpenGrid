#include "opengrid/core/util/gauss_kronrod.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

namespace og {
namespace {

// Abscissae (positive half, descending) and weights of the 15-point Kronrod rule; the Gauss
// 7-point rule uses the odd-indexed abscissae xgk[1], xgk[3], xgk[5] and the centre xgk[7].
constexpr std::array<double, 8> kXgk = {
    0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
    0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
    0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
    0.207784955007898467600689403773245, 0.000000000000000000000000000000000};
constexpr std::array<double, 8> kWgk = {
    0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
    0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
    0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
    0.204432940075298892414161999234649, 0.209482141084727828012999174891714};
constexpr std::array<double, 4> kWg = {
    0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
    0.381830050505118944950369775488975, 0.417959183673469387755102040816327};

struct Estimate {
  double kronrod;
  double error;
  double abs_integral;  // integral of |f|, the scale for the round-off floor
};

Estimate qk15(const std::function<double(double)>& f, double a, double b) {
  const double center = 0.5 * (a + b);
  const double half = 0.5 * (b - a);
  const double f_center = f(center);
  double kronrod = f_center * kWgk[7];
  double gauss = f_center * kWg[3];
  double abs_sum = std::fabs(f_center) * kWgk[7];
  for (std::size_t j = 0; j < 7; ++j) {
    const double dx = half * kXgk[j];
    const double f_lo = f(center - dx);
    const double f_hi = f(center + dx);
    const double f_sum = f_lo + f_hi;
    kronrod += kWgk[j] * f_sum;
    abs_sum += kWgk[j] * (std::fabs(f_lo) + std::fabs(f_hi));
    if (j % 2 == 1) {
      gauss += kWg[j / 2] * f_sum;
    }
  }
  return {kronrod * half, std::fabs((kronrod - gauss) * half), abs_sum * std::fabs(half)};
}

}  // namespace

double integrate_gk15(const std::function<double(double)>& f, double a, double b,
                      const QuadratureOptions& options) {
  if (a == b) {
    return 0.0;
  }
  if (a > b) {
    return -integrate_gk15(f, b, a, options);
  }
  struct Interval {
    double a;
    double b;
    Estimate estimate;
  };
  // Max-heap on the error estimate; ties broken by the left end for determinism.
  auto worse = [](const Interval& x, const Interval& y) {
    if (x.estimate.error != y.estimate.error) {
      return x.estimate.error < y.estimate.error;
    }
    return x.a > y.a;
  };
  std::vector<Interval> heap;
  heap.push_back({a, b, qk15(f, a, b)});
  double total = heap.front().estimate.kronrod;
  double error = heap.front().estimate.error;
  double abs_total = heap.front().estimate.abs_integral;
  while (std::cmp_less(heap.size(), options.max_intervals)) {
    double tolerance = options.rel_tol * std::fabs(total);
    if (options.abs_tol > tolerance) {
      tolerance = options.abs_tol;
    }
    const double roundoff = 50.0 * std::numeric_limits<double>::epsilon() * abs_total;
    if (roundoff > tolerance) {
      tolerance = roundoff;
    }
    if (error <= tolerance) {
      break;
    }
    std::pop_heap(heap.begin(), heap.end(), worse);
    const Interval worst = heap.back();
    heap.pop_back();
    const double mid = 0.5 * (worst.a + worst.b);
    if (mid <= worst.a || mid >= worst.b) {
      heap.push_back(worst);
      std::push_heap(heap.begin(), heap.end(), worse);
      break;  // interval no longer representable
    }
    const Interval left{worst.a, mid, qk15(f, worst.a, mid)};
    const Interval right{mid, worst.b, qk15(f, mid, worst.b)};
    total += left.estimate.kronrod + right.estimate.kronrod - worst.estimate.kronrod;
    error += left.estimate.error + right.estimate.error - worst.estimate.error;
    abs_total +=
        left.estimate.abs_integral + right.estimate.abs_integral - worst.estimate.abs_integral;
    heap.push_back(left);
    std::push_heap(heap.begin(), heap.end(), worse);
    heap.push_back(right);
    std::push_heap(heap.begin(), heap.end(), worse);
  }
  // Re-sum in position order so the result does not depend on the update history.
  std::sort(heap.begin(), heap.end(),
            [](const Interval& x, const Interval& y) { return x.a < y.a; });
  double sum = 0.0;
  for (const Interval& interval : heap) {
    sum += interval.estimate.kronrod;
  }
  return sum;
}

}  // namespace og
