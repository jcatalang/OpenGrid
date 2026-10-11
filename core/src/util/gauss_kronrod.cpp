#include "opengrid/core/util/gauss_kronrod.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>

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

double adapt(const std::function<double(double)>& f, double a, double b, const Estimate& whole,
             const QuadratureOptions& options, int depth) {
  double tolerance = options.rel_tol * std::fabs(whole.kronrod);
  if (options.abs_tol > tolerance) {
    tolerance = options.abs_tol;
  }
  // Round-off floor (QUADPACK): an error estimate below ~50 eps |f| cannot be improved.
  const double roundoff = 50.0 * std::numeric_limits<double>::epsilon() * whole.abs_integral;
  if (roundoff > tolerance) {
    tolerance = roundoff;
  }
  if (whole.error <= tolerance || depth >= options.max_depth) {
    return whole.kronrod;
  }
  const double mid = 0.5 * (a + b);
  if (mid <= a || mid >= b) {
    return whole.kronrod;  // interval no longer representable
  }
  const Estimate left = qk15(f, a, mid);
  const Estimate right = qk15(f, mid, b);
  // Converged when the refined estimate agrees with the coarse one to the tolerance.
  if (std::fabs(left.kronrod + right.kronrod - whole.kronrod) <= tolerance) {
    return left.kronrod + right.kronrod;
  }
  return adapt(f, a, mid, left, options, depth + 1) + adapt(f, mid, b, right, options, depth + 1);
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
  return adapt(f, a, b, qk15(f, a, b), options, 0);
}

}  // namespace og
