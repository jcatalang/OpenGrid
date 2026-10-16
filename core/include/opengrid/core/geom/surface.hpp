// Parametric surfaces x(u, v), (u, v) in [0, 1]^2, and closest-point projection.
//
// Types:
//  * PlaneSurface      x = o + u a + v b (a parallelogram patch);
//  * RuledSurface      x = (1 - v) c0(u) + v c1(u);
//  * ExtrudedSurface   x = c(u) + v d;
//  * RevolvedSurface   c(u) rotated about an axis by angle0 + v (angle1 - angle0) (Rodrigues);
//  * NurbsSurface      tensor-product NURBS of degrees (p, q) (Piegl & Tiller, The NURBS Book,
//                      2nd ed. (1997), algorithms A2.2 / A2.3 and eq. 4.20); each knot vector's
//                      valid domain is mapped linearly onto [0, 1].
// Projection: the closest point to p is searched on a grid of samples, then refined by Newton's
// method on the stationarity conditions (x - p) . x_u = (x - p) . x_v = 0 with a finite-difference
// Hessian of the squared distance (the converged point depends only on the analytic gradient),
// clamped to the parameter square, until the step is below 1e-15. Points on the surface are
// recovered to about 1e-14 relative to the surface size.
#pragma once

#include <memory>
#include <vector>

#include "opengrid/core/geom/curve.hpp"
#include "opengrid/core/geom/vec3.hpp"

namespace og {

struct SurfaceProjection {
  double u = 0.0;
  double v = 0.0;
  Vec3 point;
  double distance = 0.0;
};

// A point of a surface with its partial derivatives.
struct SurfacePoint {
  Vec3 s;
  Vec3 su;
  Vec3 sv;
};

class Surface {
 public:
  Surface() = default;
  Surface(const Surface&) = delete;
  Surface& operator=(const Surface&) = delete;
  virtual ~Surface() = default;

  [[nodiscard]] virtual Vec3 eval(double u, double v) const = 0;
  // Partial derivatives; the default uses central differences of eval.
  [[nodiscard]] virtual Vec3 du(double u, double v) const;
  [[nodiscard]] virtual Vec3 dv(double u, double v) const;
  // eval, du and dv together (one evaluation where the surface can share the work).
  [[nodiscard]] virtual SurfacePoint values(double u, double v) const;
  // Unit normal x_u x x_v (zero where the surface is degenerate).
  [[nodiscard]] Vec3 normal(double u, double v) const;

  // Closest point; `samples` per direction for the global search.
  [[nodiscard]] SurfaceProjection project(const Vec3& p, int samples = 24) const;
  // Projection refined from a given starting parameter (no global search).
  [[nodiscard]] SurfaceProjection project_from(const Vec3& p, double u0, double v0) const;
  // Size of the surface: diagonal of the bounding box of a sample grid.
  [[nodiscard]] double size() const;
};

class PlaneSurface final : public Surface {
 public:
  // Throws std::invalid_argument if a and b are parallel or zero.
  PlaneSurface(const Vec3& origin, const Vec3& a, const Vec3& b);
  [[nodiscard]] Vec3 eval(double u, double v) const override;
  [[nodiscard]] Vec3 du(double u, double v) const override;
  [[nodiscard]] Vec3 dv(double u, double v) const override;

 private:
  Vec3 origin_;
  Vec3 a_;
  Vec3 b_;
};

class RuledSurface final : public Surface {
 public:
  RuledSurface(std::shared_ptr<const Curve> c0, std::shared_ptr<const Curve> c1);
  [[nodiscard]] Vec3 eval(double u, double v) const override;
  [[nodiscard]] Vec3 du(double u, double v) const override;
  [[nodiscard]] Vec3 dv(double u, double v) const override;

 private:
  std::shared_ptr<const Curve> c0_;
  std::shared_ptr<const Curve> c1_;
};

class ExtrudedSurface final : public Surface {
 public:
  // Throws std::invalid_argument for a zero direction.
  ExtrudedSurface(std::shared_ptr<const Curve> curve, const Vec3& direction);
  [[nodiscard]] Vec3 eval(double u, double v) const override;
  [[nodiscard]] Vec3 du(double u, double v) const override;
  [[nodiscard]] Vec3 dv(double u, double v) const override;

 private:
  std::shared_ptr<const Curve> curve_;
  Vec3 direction_;
};

class RevolvedSurface final : public Surface {
 public:
  // Angles in radians. Throws std::invalid_argument for a zero axis or equal angles.
  RevolvedSurface(std::shared_ptr<const Curve> curve, const Vec3& axis_point, const Vec3& axis,
                  double angle0, double angle1);
  [[nodiscard]] Vec3 eval(double u, double v) const override;
  [[nodiscard]] Vec3 du(double u, double v) const override;
  [[nodiscard]] Vec3 dv(double u, double v) const override;

 private:
  [[nodiscard]] Vec3 rotate(const Vec3& p, double angle) const;
  std::shared_ptr<const Curve> curve_;
  Vec3 point_;
  Vec3 axis_;
  double angle0_;
  double angle1_;
};

class NurbsSurface final : public Surface {
 public:
  // control[i * nv + j] with i along u (nu points) and j along v (nv points). Throws
  // std::invalid_argument on inconsistent sizes, decreasing knots or non-positive weights.
  NurbsSurface(int degree_u, int degree_v, std::vector<double> knots_u, std::vector<double> knots_v,
               int nu, int nv, std::vector<Vec3> control, std::vector<double> weights);
  [[nodiscard]] Vec3 eval(double u, double v) const override;
  [[nodiscard]] Vec3 du(double u, double v) const override;
  [[nodiscard]] Vec3 dv(double u, double v) const override;
  [[nodiscard]] SurfacePoint values(double u, double v) const override;

 private:
  [[nodiscard]] SurfacePoint evaluate(double u, double v) const;
  int pu_;
  int pv_;
  std::vector<double> ku_;
  std::vector<double> kv_;
  int nu_;
  int nv_;
  std::vector<Vec3> control_;
  std::vector<double> weights_;
};

}  // namespace og
