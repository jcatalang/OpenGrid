#include "opengrid/core/geom/surface.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace og {
namespace {

constexpr double kFd = 1e-6;

std::size_t idx(int i) {
  return static_cast<std::size_t>(i);
}

// B-spline basis functions and first derivatives (Piegl & Tiller A2.3, n = 1) on knot span
// `span`: values[j] = N_{span-p+j}(u), derivs[j] = N'_{span-p+j}(u).
void basis(const std::vector<double>& knots, int p, int span, double u, std::vector<double>& values,
           std::vector<double>& derivs) {
  std::vector<std::vector<double>> ndu(idx(p + 1), std::vector<double>(idx(p + 1), 0.0));
  std::vector<double> left(idx(p + 1), 0.0);
  std::vector<double> right(idx(p + 1), 0.0);
  ndu[0][0] = 1.0;
  for (int j = 1; j <= p; ++j) {
    left[idx(j)] = u - knots[idx(span + 1 - j)];
    right[idx(j)] = knots[idx(span + j)] - u;
    double saved = 0.0;
    for (int r = 0; r < j; ++r) {
      ndu[idx(j)][idx(r)] = right[idx(r + 1)] + left[idx(j - r)];
      double temp = 0.0;
      if (ndu[idx(j)][idx(r)] != 0.0) {
        temp = ndu[idx(r)][idx(j - 1)] / ndu[idx(j)][idx(r)];
      }
      ndu[idx(r)][idx(j)] = saved + right[idx(r + 1)] * temp;
      saved = left[idx(j - r)] * temp;
    }
    ndu[idx(j)][idx(j)] = saved;
  }
  values.assign(idx(p + 1), 0.0);
  derivs.assign(idx(p + 1), 0.0);
  for (int j = 0; j <= p; ++j) {
    values[idx(j)] = ndu[idx(j)][idx(p)];
  }
  if (p == 0) {
    return;
  }
  for (int r = 0; r <= p; ++r) {
    // First derivative: p * (N_{r,p-1} / (u_{i+p} - u_i) - N_{r+1,p-1} / (u_{i+p+1} - u_{i+1})).
    double d = 0.0;
    if (r >= 1) {
      d += ndu[idx(r - 1)][idx(p - 1)] / ndu[idx(p)][idx(r - 1)];
    }
    if (r <= p - 1) {
      d -= ndu[idx(r)][idx(p - 1)] / ndu[idx(p)][idx(r)];
    }
    derivs[idx(r)] = p * d;
  }
}

int find_span(const std::vector<double>& knots, int p, int n_ctrl, double u) {
  const int n = n_ctrl - 1;
  if (u >= knots[idx(n + 1)]) {
    int span = n;
    while (span > p && knots[idx(span)] >= knots[idx(span + 1)]) {
      span -= 1;
    }
    return span;
  }
  int span = p;
  while (span < n && knots[idx(span + 1)] <= u) {
    span += 1;
  }
  return span;
}

void check_knots(const std::vector<double>& knots, int p, int n_ctrl, const char* which) {
  if (p < 1 || n_ctrl < p + 1 ||
      knots.size() != static_cast<std::size_t>(n_ctrl) + static_cast<std::size_t>(p) + 1) {
    throw std::invalid_argument(std::string("nurbs surface: inconsistent ") + which +
                                " degree / knots / control points");
  }
  for (std::size_t k = 1; k < knots.size(); ++k) {
    if (knots[k] < knots[k - 1]) {
      throw std::invalid_argument(std::string("nurbs surface: decreasing ") + which + " knots");
    }
  }
  if (!(knots[idx(n_ctrl)] > knots[idx(p)])) {
    throw std::invalid_argument(std::string("nurbs surface: empty ") + which + " domain");
  }
}

}  // namespace

// ------------------------------------------------------------------ Surface

Vec3 Surface::du(double u, double v) const {
  const double a = std::max(0.0, u - kFd);
  const double b = std::min(1.0, u + kFd);
  return (eval(b, v) - eval(a, v)) / (b - a);
}

Vec3 Surface::dv(double u, double v) const {
  const double a = std::max(0.0, v - kFd);
  const double b = std::min(1.0, v + kFd);
  return (eval(u, b) - eval(u, a)) / (b - a);
}

Vec3 Surface::normal(double u, double v) const {
  const Vec3 n = cross(du(u, v), dv(u, v));
  const double length = norm(n);
  if (length == 0.0) {
    return n;
  }
  return n / length;
}

double Surface::size() const {
  Vec3 lo = eval(0.0, 0.0);
  Vec3 hi = lo;
  for (int i = 0; i <= 8; ++i) {
    for (int j = 0; j <= 8; ++j) {
      const Vec3 p = eval(i / 8.0, j / 8.0);
      lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
      hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
    }
  }
  return norm(hi - lo);
}

SurfaceProjection Surface::project(const Vec3& p, int samples) const {
  samples = std::max(samples, 2);
  double best = std::numeric_limits<double>::infinity();
  double bu = 0.0;
  double bv = 0.0;
  for (int i = 0; i <= samples; ++i) {
    for (int j = 0; j <= samples; ++j) {
      const double u = static_cast<double>(i) / samples;
      const double v = static_cast<double>(j) / samples;
      const double d = distance(eval(u, v), p);
      if (d < best) {
        best = d;
        bu = u;
        bv = v;
      }
    }
  }
  return project_from(p, bu, bv);
}

SurfaceProjection Surface::project_from(const Vec3& p, double u, double v) const {
  auto gradient = [&](double uu, double vv, double& gu, double& gv) {
    const Vec3 r = eval(uu, vv) - p;
    gu = dot(r, du(uu, vv));
    gv = dot(r, dv(uu, vv));
  };
  for (int iteration = 0; iteration < 100; ++iteration) {
    double gu = 0.0;
    double gv = 0.0;
    gradient(u, v, gu, gv);
    // Hessian of 0.5 |x - p|^2 by central differences of the analytic gradient.
    const double hu = std::max(1e-7, 1e-6 * std::max(u, 1.0 - u));
    const double hv = std::max(1e-7, 1e-6 * std::max(v, 1.0 - v));
    const double u0 = std::max(0.0, u - hu);
    const double u1 = std::min(1.0, u + hu);
    const double v0 = std::max(0.0, v - hv);
    const double v1 = std::min(1.0, v + hv);
    double a0 = 0.0;
    double b0 = 0.0;
    double a1 = 0.0;
    double b1 = 0.0;
    gradient(u0, v, a0, b0);
    gradient(u1, v, a1, b1);
    const double huu = (a1 - a0) / (u1 - u0);
    const double hvu = (b1 - b0) / (u1 - u0);
    gradient(u, v0, a0, b0);
    gradient(u, v1, a1, b1);
    const double huv = (a1 - a0) / (v1 - v0);
    const double hvv = (b1 - b0) / (v1 - v0);
    const double h12 = 0.5 * (hvu + huv);
    const double det = huu * hvv - h12 * h12;
    double su = 0.0;
    double sv = 0.0;
    if (det > 0.0 && huu > 0.0) {
      su = -(hvv * gu - h12 * gv) / det;
      sv = -(huu * gv - h12 * gu) / det;
    } else {
      // Not convex here: Gauss-Newton step (J^T J positive semi-definite).
      const Vec3 xu = du(u, v);
      const Vec3 xv = dv(u, v);
      const double g11 = dot(xu, xu);
      const double g12 = dot(xu, xv);
      const double g22 = dot(xv, xv);
      const double gd = g11 * g22 - g12 * g12;
      if (gd <= 0.0) {
        break;
      }
      su = -(g22 * gu - g12 * gv) / gd;
      sv = -(g11 * gv - g12 * gu) / gd;
    }
    // Damped, clamped step: never increase the distance.
    const double d0 = distance(eval(u, v), p);
    double factor = 1.0;
    double nu = u;
    double nv = v;
    for (int halving = 0; halving < 30; ++halving) {
      nu = std::clamp(u + factor * su, 0.0, 1.0);
      nv = std::clamp(v + factor * sv, 0.0, 1.0);
      if (distance(eval(nu, nv), p) <= d0) {
        break;
      }
      factor *= 0.5;
    }
    const double step = std::fabs(nu - u) + std::fabs(nv - v);
    u = nu;
    v = nv;
    if (step < 1e-15) {
      break;
    }
  }
  SurfaceProjection out;
  out.u = u;
  out.v = v;
  out.point = eval(u, v);
  out.distance = distance(out.point, p);
  return out;
}

// ------------------------------------------------------------------ plane

PlaneSurface::PlaneSurface(const Vec3& origin, const Vec3& a, const Vec3& b)
    : origin_(origin), a_(a), b_(b) {
  if (norm(cross(a, b)) <= 1e-14 * norm(a) * norm(b) || norm(a) == 0.0 || norm(b) == 0.0) {
    throw std::invalid_argument("plane: the two edge vectors must be non-zero and not parallel");
  }
}

Vec3 PlaneSurface::eval(double u, double v) const {
  return origin_ + u * a_ + v * b_;
}
Vec3 PlaneSurface::du(double /*u*/, double /*v*/) const {
  return a_;
}
Vec3 PlaneSurface::dv(double /*u*/, double /*v*/) const {
  return b_;
}

// ------------------------------------------------------------------ ruled

RuledSurface::RuledSurface(std::shared_ptr<const Curve> c0, std::shared_ptr<const Curve> c1)
    : c0_(std::move(c0)), c1_(std::move(c1)) {
  if (!c0_ || !c1_) {
    throw std::invalid_argument("ruled surface: null curve");
  }
}

Vec3 RuledSurface::eval(double u, double v) const {
  return (1.0 - v) * c0_->eval(u) + v * c1_->eval(u);
}
Vec3 RuledSurface::du(double u, double v) const {
  return (1.0 - v) * c0_->d1(u) + v * c1_->d1(u);
}
Vec3 RuledSurface::dv(double u, double /*v*/) const {
  return c1_->eval(u) - c0_->eval(u);
}

// ------------------------------------------------------------------ extruded

ExtrudedSurface::ExtrudedSurface(std::shared_ptr<const Curve> curve, const Vec3& direction)
    : curve_(std::move(curve)), direction_(direction) {
  if (!curve_ || norm(direction) == 0.0) {
    throw std::invalid_argument("extruded surface: null curve or zero direction");
  }
}

Vec3 ExtrudedSurface::eval(double u, double v) const {
  return curve_->eval(u) + v * direction_;
}
Vec3 ExtrudedSurface::du(double u, double /*v*/) const {
  return curve_->d1(u);
}
Vec3 ExtrudedSurface::dv(double /*u*/, double /*v*/) const {
  return direction_;
}

// ------------------------------------------------------------------ revolved

RevolvedSurface::RevolvedSurface(std::shared_ptr<const Curve> curve, const Vec3& axis_point,
                                 const Vec3& axis, double angle0, double angle1)
    : curve_(std::move(curve)), point_(axis_point), angle0_(angle0), angle1_(angle1) {
  if (!curve_ || norm(axis) == 0.0 || angle0 == angle1) {
    throw std::invalid_argument("revolved surface: null curve, zero axis or equal angles");
  }
  axis_ = axis / norm(axis);
}

Vec3 RevolvedSurface::rotate(const Vec3& p, double angle) const {
  const Vec3 r = p - point_;
  const double c = std::cos(angle);
  const double s = std::sin(angle);
  return point_ + c * r + s * cross(axis_, r) + (1.0 - c) * dot(axis_, r) * axis_;
}

Vec3 RevolvedSurface::eval(double u, double v) const {
  return rotate(curve_->eval(u), angle0_ + v * (angle1_ - angle0_));
}

Vec3 RevolvedSurface::du(double u, double v) const {
  // Rotation is linear: rotate the tangent (as a vector, about the origin).
  const double angle = angle0_ + v * (angle1_ - angle0_);
  const Vec3 t = curve_->d1(u);
  const double c = std::cos(angle);
  const double s = std::sin(angle);
  return c * t + s * cross(axis_, t) + (1.0 - c) * dot(axis_, t) * axis_;
}

Vec3 RevolvedSurface::dv(double u, double v) const {
  const double angle = angle0_ + v * (angle1_ - angle0_);
  const Vec3 r = curve_->eval(u) - point_;
  const double c = std::cos(angle);
  const double s = std::sin(angle);
  // d/dangle of the Rodrigues formula.
  const Vec3 d = -s * r + c * cross(axis_, r) + s * dot(axis_, r) * axis_;
  return (angle1_ - angle0_) * d;
}

// ------------------------------------------------------------------ NURBS

NurbsSurface::NurbsSurface(int degree_u, int degree_v, std::vector<double> knots_u,
                           std::vector<double> knots_v, int nu, int nv, std::vector<Vec3> control,
                           std::vector<double> weights)
    : pu_(degree_u),
      pv_(degree_v),
      ku_(std::move(knots_u)),
      kv_(std::move(knots_v)),
      nu_(nu),
      nv_(nv),
      control_(std::move(control)),
      weights_(std::move(weights)) {
  check_knots(ku_, pu_, nu_, "u");
  check_knots(kv_, pv_, nv_, "v");
  const auto n = static_cast<std::size_t>(nu_) * static_cast<std::size_t>(nv_);
  if (control_.size() != n || weights_.size() != n) {
    throw std::invalid_argument("nurbs surface: need nu * nv control points and weights");
  }
  for (const double w : weights_) {
    if (!(w > 0.0)) {
      throw std::invalid_argument("nurbs surface: weights must be positive");
    }
  }
}

NurbsSurface::Values NurbsSurface::evaluate(double u, double v) const {
  const double lo_u = ku_[idx(pu_)];
  const double hi_u = ku_[idx(nu_)];
  const double lo_v = kv_[idx(pv_)];
  const double hi_v = kv_[idx(nv_)];
  const double uu = lo_u + std::clamp(u, 0.0, 1.0) * (hi_u - lo_u);
  const double vv = lo_v + std::clamp(v, 0.0, 1.0) * (hi_v - lo_v);
  const int su = find_span(ku_, pu_, nu_, uu);
  const int sv = find_span(kv_, pv_, nv_, vv);
  std::vector<double> nu_val;
  std::vector<double> nu_der;
  std::vector<double> nv_val;
  std::vector<double> nv_der;
  basis(ku_, pu_, su, uu, nu_val, nu_der);
  basis(kv_, pv_, sv, vv, nv_val, nv_der);
  // Homogeneous sums: A = sum N M w P, W = sum N M w, and their u / v derivatives.
  Vec3 a;
  Vec3 au;
  Vec3 av;
  double w = 0.0;
  double wu = 0.0;
  double wv = 0.0;
  for (int i = 0; i <= pu_; ++i) {
    for (int j = 0; j <= pv_; ++j) {
      const std::size_t k = idx(su - pu_ + i) * static_cast<std::size_t>(nv_) + idx(sv - pv_ + j);
      const double wk = weights_[k];
      const Vec3& pk = control_[k];
      const double b = nu_val[idx(i)] * nv_val[idx(j)] * wk;
      const double bu = nu_der[idx(i)] * nv_val[idx(j)] * wk;
      const double bv = nu_val[idx(i)] * nv_der[idx(j)] * wk;
      a += b * pk;
      au += bu * pk;
      av += bv * pk;
      w += b;
      wu += bu;
      wv += bv;
    }
  }
  Values out;
  out.s = a / w;
  // Chain rule for the mapping of [0, 1] onto the knot domain.
  out.su = (au - wu * out.s) / w * (hi_u - lo_u);
  out.sv = (av - wv * out.s) / w * (hi_v - lo_v);
  return out;
}

Vec3 NurbsSurface::eval(double u, double v) const {
  return evaluate(u, v).s;
}
Vec3 NurbsSurface::du(double u, double v) const {
  return evaluate(u, v).su;
}
Vec3 NurbsSurface::dv(double u, double v) const {
  return evaluate(u, v).sv;
}

}  // namespace og
