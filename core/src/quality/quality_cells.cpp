#include "opengrid/core/quality/quality_cells.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
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

// Corner indices of cell c (triangles first, then quads); n = 3 or 4.
struct Cells {
  const std::vector<std::array<std::int64_t, 3>>& triangles;
  const std::vector<std::array<std::int64_t, 4>>& quads;

  [[nodiscard]] std::size_t size() const {
    return triangles.size() + quads.size();
  }
  [[nodiscard]] std::size_t corners(std::size_t c) const {
    if (c < triangles.size()) {
      return 3;
    }
    return 4;
  }
  [[nodiscard]] std::int64_t at(std::size_t c, std::size_t k) const {
    if (c < triangles.size()) {
      return triangles[c][k];
    }
    return quads[c - triangles.size()][k];
  }
};

}  // namespace

Quality2D quality_cells(const std::vector<std::array<double, 2>>& points,
                        const std::vector<std::array<std::int64_t, 3>>& triangles,
                        const std::vector<std::array<std::int64_t, 4>>& quads) {
  const Cells cells{triangles, quads};
  const std::size_t n = cells.size();
  Quality2D q;
  for (auto* field : {&q.scaled_jacobian, &q.skewness, &q.aspect_ratio, &q.min_angle, &q.max_angle,
                      &q.non_orthogonality, &q.area_ratio, &q.area}) {
    field->assign(n, 0.0);
  }
  std::vector<std::array<double, 2>> centroid(n);
  const double sin60 = std::numbers::sqrt3 / 2.0;
  for (std::size_t c = 0; c < n; ++c) {
    const std::size_t m = cells.corners(c);
    std::array<std::array<double, 2>, 4> p{};
    for (std::size_t k = 0; k < m; ++k) {
      const std::int64_t v = cells.at(c, k);
      if (v < 0 || u(v) >= points.size()) {
        throw std::invalid_argument("quality: cell refers to a missing point");
      }
      p[k] = points[u(v)];
    }
    double ideal = 90.0;
    double jac_scale = 1.0;
    if (m == 3) {
      ideal = 60.0;
      jac_scale = 1.0 / sin60;
    }
    double min_jac = 1e300;
    double min_angle = 360.0;
    double max_angle = 0.0;
    std::array<double, 4> length{};
    for (std::size_t k = 0; k < m; ++k) {
      const auto& o = p[k];
      const auto& e1 = p[(k + 1) % m];
      const auto& e2 = p[(k + m - 1) % m];
      const double ux = e1[0] - o[0];
      const double uy = e1[1] - o[1];
      const double vx = e2[0] - o[0];
      const double vy = e2[1] - o[1];
      const double lu = std::hypot(ux, uy);
      const double lv = std::hypot(vx, vy);
      length[k] = lu;
      double jac = 0.0;
      double angle = 0.0;
      if (lu > 0.0 && lv > 0.0) {
        const double cz = ux * vy - uy * vx;
        jac = cz / (lu * lv) * jac_scale;
        angle = std::atan2(cz, ux * vx + uy * vy) * kDeg;
        if (angle < 0.0) {
          angle += 360.0;
        }
      }
      min_jac = std::min(min_jac, jac);
      min_angle = std::min(min_angle, angle);
      max_angle = std::max(max_angle, angle);
    }
    // Shoelace area and area centroid.
    double area2 = 0.0;
    double cx = 0.0;
    double cy = 0.0;
    for (std::size_t k = 0; k < m; ++k) {
      const auto& a = p[k];
      const auto& b = p[(k + 1) % m];
      const double cross = a[0] * b[1] - b[0] * a[1];
      area2 += cross;
      cx += (a[0] + b[0]) * cross;
      cy += (a[1] + b[1]) * cross;
    }
    const double area = 0.5 * area2;
    if (area2 != 0.0) {
      centroid[c] = {cx / (3.0 * area2), cy / (3.0 * area2)};
    } else {
      double sx = 0.0;
      double sy = 0.0;
      for (std::size_t k = 0; k < m; ++k) {
        sx += p[k][0];
        sy += p[k][1];
      }
      centroid[c] = {sx / static_cast<double>(m), sy / static_cast<double>(m)};
    }
    q.area[c] = area;
    q.scaled_jacobian[c] = min_jac;
    q.min_angle[c] = min_angle;
    q.max_angle[c] = max_angle;
    q.skewness[c] = std::max((max_angle - ideal) / (180.0 - ideal), (ideal - min_angle) / ideal);
    if (m == 3) {
      const double la = length[1];
      const double lb = length[2];
      const double lc = length[0];
      if (area > 0.0) {
        const double circumradius = la * lb * lc / (4.0 * area);
        const double inradius = 2.0 * area / (la + lb + lc);
        q.aspect_ratio[c] = circumradius / (2.0 * inradius);
      } else {
        q.aspect_ratio[c] = std::numeric_limits<double>::infinity();
      }
    } else {
      const double l0 = 0.5 * (length[0] + length[2]);
      const double l1 = 0.5 * (length[1] + length[3]);
      const double small = std::min(l0, l1);
      if (small > 0.0) {
        q.aspect_ratio[c] = std::max(l0, l1) / small;
      } else {
        q.aspect_ratio[c] = std::numeric_limits<double>::infinity();
      }
    }
  }
  // Edge neighbours from sorted half-edges: (min, max, cell, corner).
  struct Half {
    std::int64_t a;
    std::int64_t b;
    std::size_t cell;
    std::size_t corner;
  };
  std::vector<Half> half;
  half.reserve(3 * triangles.size() + 4 * quads.size());
  for (std::size_t c = 0; c < n; ++c) {
    const std::size_t m = cells.corners(c);
    for (std::size_t k = 0; k < m; ++k) {
      std::int64_t a = cells.at(c, k);
      std::int64_t b = cells.at(c, (k + 1) % m);
      if (b < a) {
        std::swap(a, b);
      }
      half.push_back({a, b, c, k});
    }
  }
  std::sort(half.begin(), half.end(), [](const Half& x, const Half& y) {
    if (x.a != y.a) {
      return x.a < y.a;
    }
    if (x.b != y.b) {
      return x.b < y.b;
    }
    return x.cell < y.cell;
  });
  // neighbour[cell][corner] = cell across the edge (corner, corner + 1), or n on the boundary.
  std::vector<std::array<std::size_t, 4>> neighbour(n, {n, n, n, n});
  for (std::size_t i = 0; i + 1 < half.size(); ++i) {
    if (half[i].a == half[i + 1].a && half[i].b == half[i + 1].b) {
      neighbour[half[i].cell][half[i].corner] = half[i + 1].cell;
      neighbour[half[i + 1].cell][half[i + 1].corner] = half[i].cell;
      ++i;
    }
  }
  for (std::size_t c = 0; c < n; ++c) {
    const std::size_t m = cells.corners(c);
    double worst_angle = 0.0;
    double worst_ratio = 1.0;
    for (std::size_t k = 0; k < m; ++k) {
      const auto& pa = points[u(cells.at(c, k))];
      const auto& pb = points[u(cells.at(c, (k + 1) % m))];
      const double nx = pb[1] - pa[1];
      const double ny = -(pb[0] - pa[0]);
      const std::size_t other = neighbour[c][k];
      double dx = 0.0;
      double dy = 0.0;
      if (other < n) {
        dx = centroid[other][0] - centroid[c][0];
        dy = centroid[other][1] - centroid[c][1];
        const double ao = std::fabs(q.area[other]);
        const double ac = std::fabs(q.area[c]);
        if (ao > 0.0 && ac > 0.0) {
          worst_ratio = std::max({worst_ratio, ao / ac, ac / ao});
        }
      } else {
        dx = 0.5 * (pa[0] + pb[0]) - centroid[c][0];
        dy = 0.5 * (pa[1] + pb[1]) - centroid[c][1];
      }
      worst_angle = std::max(worst_angle, angle_between(nx, ny, dx, dy));
    }
    q.non_orthogonality[c] = worst_angle;
    q.area_ratio[c] = worst_ratio;
  }
  return q;
}

}  // namespace og
