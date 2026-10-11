#include "opengrid/core/quality/quality_tri.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numbers>
#include <stdexcept>
#include <utility>

namespace og {
namespace {

constexpr double kDeg = 180.0 / std::numbers::pi;

std::size_t u(std::int64_t i) {
  return static_cast<std::size_t>(i);
}

double angle_between(double ax, double ay, double bx, double by) {
  const double na = std::hypot(ax, ay);
  const double nb = std::hypot(bx, by);
  if (na == 0.0 || nb == 0.0) {
    return 90.0;
  }
  return std::acos(std::clamp(std::fabs(ax * bx + ay * by) / (na * nb), 0.0, 1.0)) * kDeg;
}

}  // namespace

Quality2D quality_triangles(const std::vector<std::array<double, 2>>& points,
                            const std::vector<std::array<std::int64_t, 3>>& triangles) {
  const std::size_t n = triangles.size();
  Quality2D q;
  for (auto* field : {&q.scaled_jacobian, &q.skewness, &q.aspect_ratio, &q.min_angle, &q.max_angle,
                      &q.non_orthogonality, &q.area_ratio, &q.area}) {
    field->assign(n, 0.0);
  }
  std::vector<std::array<double, 2>> centroid(n);
  const double sin60 = std::numbers::sqrt3 / 2.0;
  for (std::size_t t = 0; t < n; ++t) {
    for (const auto v : triangles[t]) {
      if (v < 0 || u(v) >= points.size()) {
        throw std::invalid_argument("quality: triangle refers to a missing point");
      }
    }
    const auto& a = points[u(triangles[t][0])];
    const auto& b = points[u(triangles[t][1])];
    const auto& c = points[u(triangles[t][2])];
    const std::array<const std::array<double, 2>*, 3> p = {&a, &b, &c};
    double min_jac = 1e300;
    double min_angle = 360.0;
    double max_angle = 0.0;
    for (std::size_t k = 0; k < 3; ++k) {
      const auto& o = *p[k];
      const auto& e1 = *p[(k + 1) % 3];
      const auto& e2 = *p[(k + 2) % 3];
      const double ux = e1[0] - o[0];
      const double uy = e1[1] - o[1];
      const double vx = e2[0] - o[0];
      const double vy = e2[1] - o[1];
      const double lu = std::hypot(ux, uy);
      const double lv = std::hypot(vx, vy);
      double jac = 0.0;
      double angle = 0.0;
      if (lu > 0.0 && lv > 0.0) {
        const double cz = ux * vy - uy * vx;
        jac = cz / (lu * lv) / sin60;
        angle = std::atan2(cz, ux * vx + uy * vy) * kDeg;
        if (angle < 0.0) {
          angle += 360.0;
        }
      }
      min_jac = std::min(min_jac, jac);
      min_angle = std::min(min_angle, angle);
      max_angle = std::max(max_angle, angle);
    }
    const double la = std::hypot(b[0] - c[0], b[1] - c[1]);
    const double lb = std::hypot(c[0] - a[0], c[1] - a[1]);
    const double lc = std::hypot(a[0] - b[0], a[1] - b[1]);
    const double area = 0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]));
    q.area[t] = area;
    q.scaled_jacobian[t] = min_jac;
    q.min_angle[t] = min_angle;
    q.max_angle[t] = max_angle;
    q.skewness[t] = std::max((max_angle - 60.0) / 120.0, (60.0 - min_angle) / 60.0);
    if (area > 0.0) {
      const double circumradius = la * lb * lc / (4.0 * area);
      const double inradius = 2.0 * area / (la + lb + lc);
      q.aspect_ratio[t] = circumradius / (2.0 * inradius);
    } else {
      q.aspect_ratio[t] = std::numeric_limits<double>::infinity();
    }
    centroid[t] = {(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0};
  }
  // Edge -> adjacent triangles.
  std::map<std::pair<std::int64_t, std::int64_t>, std::array<std::int64_t, 2>> edges;
  for (std::size_t t = 0; t < n; ++t) {
    for (std::size_t k = 0; k < 3; ++k) {
      auto a = triangles[t][k];
      auto b = triangles[t][(k + 1) % 3];
      if (b < a) {
        std::swap(a, b);
      }
      auto [it, inserted] = edges.try_emplace({a, b}, std::array<std::int64_t, 2>{-1, -1});
      if (it->second[0] < 0) {
        it->second[0] = static_cast<std::int64_t>(t);
      } else {
        it->second[1] = static_cast<std::int64_t>(t);
      }
    }
  }
  for (std::size_t t = 0; t < n; ++t) {
    double worst_angle = 0.0;
    double worst_ratio = 1.0;
    for (std::size_t k = 0; k < 3; ++k) {
      auto a = triangles[t][k];
      auto b = triangles[t][(k + 1) % 3];
      const auto& pa = points[u(a)];
      const auto& pb = points[u(b)];
      // Edge normal (in-plane perpendicular).
      const double nx = pb[1] - pa[1];
      const double ny = -(pb[0] - pa[0]);
      if (b < a) {
        std::swap(a, b);
      }
      const auto& pair = edges.at({a, b});
      std::int64_t other = pair[0];
      if (std::cmp_equal(other, t)) {
        other = pair[1];
      }
      double dx = 0.0;
      double dy = 0.0;
      if (other >= 0) {
        dx = centroid[u(other)][0] - centroid[t][0];
        dy = centroid[u(other)][1] - centroid[t][1];
        const double ao = std::fabs(q.area[u(other)]);
        const double at = std::fabs(q.area[t]);
        if (ao > 0.0 && at > 0.0) {
          worst_ratio = std::max({worst_ratio, ao / at, at / ao});
        }
      } else {
        dx = 0.5 * (pa[0] + pb[0]) - centroid[t][0];
        dy = 0.5 * (pa[1] + pb[1]) - centroid[t][1];
      }
      worst_angle = std::max(worst_angle, angle_between(nx, ny, dx, dy));
    }
    q.non_orthogonality[t] = worst_angle;
    q.area_ratio[t] = worst_ratio;
  }
  return q;
}

}  // namespace og
