#include "opengrid/core/util/predicates.hpp"

#include <cmath>
#include <ranges>
#include <vector>

namespace og {
namespace {

constexpr double kEpsilon = 1.1102230246251565e-16;  // 2^-53
constexpr double kSplitter = 134217729.0;            // 2^27 + 1
constexpr double kCcwErrBoundA = (3.0 + 16.0 * kEpsilon) * kEpsilon;
constexpr double kIccErrBoundA = (10.0 + 96.0 * kEpsilon) * kEpsilon;

using Expansion = std::vector<double>;  // nonoverlapping, increasing magnitude

void two_sum(double a, double b, double& x, double& y) {
  x = a + b;
  const double bv = x - a;
  const double av = x - bv;
  y = (a - av) + (b - bv);
}

void two_diff(double a, double b, double& x, double& y) {
  x = a - b;
  const double bv = a - x;
  const double av = x + bv;
  y = (a - av) + (bv - b);
}

void split(double a, double& hi, double& lo) {
  const double c = kSplitter * a;
  const double big = c - a;
  hi = c - big;
  lo = a - hi;
}

void two_product(double a, double b, double& x, double& y) {
  x = a * b;
  double ahi = 0.0;
  double alo = 0.0;
  double bhi = 0.0;
  double blo = 0.0;
  split(a, ahi, alo);
  split(b, bhi, blo);
  const double err1 = x - ahi * bhi;
  const double err2 = err1 - alo * bhi;
  const double err3 = err2 - ahi * blo;
  y = alo * blo - err3;
}

// Exact difference a - b as a two-component expansion.
Expansion difference(double a, double b) {
  double x = 0.0;
  double y = 0.0;
  two_diff(a, b, x, y);
  return {y, x};
}

// e * b (Shewchuk's scale_expansion_zeroelim).
Expansion scale(const Expansion& e, double b) {
  Expansion h;
  if (e.empty()) {
    return h;
  }
  double bhi = 0.0;
  double blo = 0.0;
  split(b, bhi, blo);
  double q = 0.0;
  double hh = 0.0;
  two_product(e[0], b, q, hh);
  if (hh != 0.0) {
    h.push_back(hh);
  }
  for (std::size_t i = 1; i < e.size(); ++i) {
    double product1 = 0.0;
    double product0 = 0.0;
    two_product(e[i], b, product1, product0);
    double sum = 0.0;
    two_sum(q, product0, sum, hh);
    if (hh != 0.0) {
      h.push_back(hh);
    }
    double q_new = 0.0;
    two_sum(product1, sum, q_new, hh);
    q = q_new;
    if (hh != 0.0) {
      h.push_back(hh);
    }
  }
  if (q != 0.0 || h.empty()) {
    h.push_back(q);
  }
  return h;
}

// e + f (linear expansion sum with zero elimination; Shewchuk's expansion_sum_zeroelim).
Expansion add(const Expansion& e, const Expansion& f) {
  Expansion h = e;
  for (const double b : f) {
    // grow_expansion(h, b)
    Expansion g;
    double q = b;
    for (const double hv : h) {
      double sum = 0.0;
      double err = 0.0;
      two_sum(q, hv, sum, err);
      q = sum;
      if (err != 0.0) {
        g.push_back(err);
      }
    }
    g.push_back(q);
    h = g;
  }
  return h;
}

Expansion negate(Expansion e) {
  for (double& v : e) {
    v = -v;
  }
  return e;
}

Expansion multiply(const Expansion& e, const Expansion& f) {
  Expansion result;
  for (const double v : f) {
    result = add(result, scale(e, v));
  }
  return result;
}

int sign(const Expansion& e) {
  // The most significant nonzero component (expansions are ordered by magnitude).
  for (const double v : std::views::reverse(e)) {
    if (v > 0.0) {
      return 1;
    }
    if (v < 0.0) {
      return -1;
    }
  }
  return 0;
}

int sign(double v) {
  if (v > 0.0) {
    return 1;
  }
  if (v < 0.0) {
    return -1;
  }
  return 0;
}

int orient2d_exact(const double* a, const double* b, const double* c) {
  const Expansion acx = difference(a[0], c[0]);
  const Expansion bcy = difference(b[1], c[1]);
  const Expansion acy = difference(a[1], c[1]);
  const Expansion bcx = difference(b[0], c[0]);
  return sign(add(multiply(acx, bcy), negate(multiply(acy, bcx))));
}

int incircle_exact(const double* a, const double* b, const double* c, const double* d) {
  const Expansion adx = difference(a[0], d[0]);
  const Expansion ady = difference(a[1], d[1]);
  const Expansion bdx = difference(b[0], d[0]);
  const Expansion bdy = difference(b[1], d[1]);
  const Expansion cdx = difference(c[0], d[0]);
  const Expansion cdy = difference(c[1], d[1]);
  const Expansion alift = add(multiply(adx, adx), multiply(ady, ady));
  const Expansion blift = add(multiply(bdx, bdx), multiply(bdy, bdy));
  const Expansion clift = add(multiply(cdx, cdx), multiply(cdy, cdy));
  const Expansion bc = add(multiply(bdx, cdy), negate(multiply(cdx, bdy)));
  const Expansion ca = add(multiply(cdx, ady), negate(multiply(adx, cdy)));
  const Expansion ab = add(multiply(adx, bdy), negate(multiply(bdx, ady)));
  return sign(add(add(multiply(alift, bc), multiply(blift, ca)), multiply(clift, ab)));
}

}  // namespace

int orient2d(const double* a, const double* b, const double* c) {
  const double detleft = (a[0] - c[0]) * (b[1] - c[1]);
  const double detright = (a[1] - c[1]) * (b[0] - c[0]);
  const double det = detleft - detright;
  const double errbound = kCcwErrBoundA * (std::fabs(detleft) + std::fabs(detright));
  if (det > errbound || -det > errbound) {
    return sign(det);
  }
  return orient2d_exact(a, b, c);
}

int incircle(const double* a, const double* b, const double* c, const double* d) {
  const double adx = a[0] - d[0];
  const double bdx = b[0] - d[0];
  const double cdx = c[0] - d[0];
  const double ady = a[1] - d[1];
  const double bdy = b[1] - d[1];
  const double cdy = c[1] - d[1];
  const double bdxcdy = bdx * cdy;
  const double cdxbdy = cdx * bdy;
  const double alift = adx * adx + ady * ady;
  const double cdxady = cdx * ady;
  const double adxcdy = adx * cdy;
  const double blift = bdx * bdx + bdy * bdy;
  const double adxbdy = adx * bdy;
  const double bdxady = bdx * ady;
  const double clift = cdx * cdx + cdy * cdy;
  const double det =
      alift * (bdxcdy - cdxbdy) + blift * (cdxady - adxcdy) + clift * (adxbdy - bdxady);
  const double permanent = (std::fabs(bdxcdy) + std::fabs(cdxbdy)) * alift +
                           (std::fabs(cdxady) + std::fabs(adxcdy)) * blift +
                           (std::fabs(adxbdy) + std::fabs(bdxady)) * clift;
  const double errbound = kIccErrBoundA * permanent;
  if (det > errbound || -det > errbound) {
    return sign(det);
  }
  return incircle_exact(a, b, c, d);
}

}  // namespace og
