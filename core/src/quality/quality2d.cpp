#include "opengrid/core/quality/quality2d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace og {
namespace {

constexpr double kDeg = 180.0 / std::numbers::pi;

double cross_z(const Vec3& a, const Vec3& b) {
  return a.x * b.y - a.y * b.x;
}

std::array<Vec3, 4> corners(const Grid2D& g, int i, int j) {
  return {g.at(i, j), g.at(i + 1, j), g.at(i + 1, j + 1), g.at(i, j + 1)};
}

// Signed area (shoelace) and area centroid of a quadrilateral.
void area_centroid(const std::array<Vec3, 4>& p, double& area, Vec3& centroid) {
  double a = 0.0;
  double cx = 0.0;
  double cy = 0.0;
  for (std::size_t k = 0; k < 4; ++k) {
    const Vec3& u = p[k];
    const Vec3& v = p[(k + 1) % 4];
    const double c = u.x * v.y - v.x * u.y;
    a += c;
    cx += (u.x + v.x) * c;
    cy += (u.y + v.y) * c;
  }
  a *= 0.5;
  area = a;
  const double z = 0.25 * (p[0].z + p[1].z + p[2].z + p[3].z);
  if (std::fabs(a) > 0.0) {
    centroid = {cx / (6.0 * a), cy / (6.0 * a), z};
  } else {
    centroid = 0.25 * (p[0] + p[1] + p[2] + p[3]);
  }
}

double angle_between(const Vec3& a, const Vec3& b) {
  const double na = norm(a);
  const double nb = norm(b);
  if (na == 0.0 || nb == 0.0) {
    return 90.0;
  }
  const double c = std::clamp(std::fabs(dot(a, b)) / (na * nb), 0.0, 1.0);
  return std::acos(c) * kDeg;
}

}  // namespace

Quality2D quality_2d(const Grid2D& g) {
  if (g.ni < 2 || g.nj < 2) {
    throw std::invalid_argument("quality: the grid needs at least 2 x 2 points");
  }
  const int ci = g.ni - 1;
  const int cj = g.nj - 1;
  const auto n = static_cast<std::size_t>(ci) * static_cast<std::size_t>(cj);
  Quality2D q;
  for (auto* field : {&q.scaled_jacobian, &q.skewness, &q.aspect_ratio, &q.min_angle, &q.max_angle,
                      &q.non_orthogonality, &q.area_ratio, &q.area}) {
    field->assign(n, 0.0);
  }
  std::vector<Vec3> centroid(n);
  auto cell = [ci](int i, int j) {
    return static_cast<std::size_t>(j) * static_cast<std::size_t>(ci) + static_cast<std::size_t>(i);
  };

  for (int j = 0; j < cj; ++j) {
    for (int i = 0; i < ci; ++i) {
      const std::size_t c = cell(i, j);
      const auto p = corners(g, i, j);
      double min_jac = 1.0;
      double min_angle = 360.0;
      double max_angle = 0.0;
      for (std::size_t k = 0; k < 4; ++k) {
        const Vec3 e1 = p[(k + 1) % 4] - p[k];
        const Vec3 e2 = p[(k + 3) % 4] - p[k];
        const double l1 = norm(e1);
        const double l2 = norm(e2);
        double jac = 0.0;
        double angle = 0.0;
        if (l1 > 0.0 && l2 > 0.0) {
          const double cz = cross_z(e1, e2);
          jac = cz / (l1 * l2);
          angle = std::atan2(cz, dot(e1, e2)) * kDeg;
          if (angle < 0.0) {
            angle += 360.0;
          }
        }
        min_jac = std::min(min_jac, jac);
        min_angle = std::min(min_angle, angle);
        max_angle = std::max(max_angle, angle);
      }
      q.scaled_jacobian[c] = min_jac;
      q.min_angle[c] = min_angle;
      q.max_angle[c] = max_angle;
      q.skewness[c] = std::max((max_angle - 90.0) / 90.0, (90.0 - min_angle) / 90.0);
      const double l0 = 0.5 * (distance(p[0], p[1]) + distance(p[3], p[2]));
      const double l1 = 0.5 * (distance(p[1], p[2]) + distance(p[0], p[3]));
      const double small = std::min(l0, l1);
      if (small > 0.0) {
        q.aspect_ratio[c] = std::max(l0, l1) / small;
      } else {
        q.aspect_ratio[c] = std::numeric_limits<double>::infinity();
      }
      area_centroid(p, q.area[c], centroid[c]);
    }
  }

  for (int j = 0; j < cj; ++j) {
    for (int i = 0; i < ci; ++i) {
      const std::size_t c = cell(i, j);
      const auto p = corners(g, i, j);
      // Faces: 0 bottom (p0 p1, neighbour j-1), 1 right (p1 p2, i+1), 2 top (p2 p3, j+1),
      // 3 left (p3 p0, i-1).
      const std::array<std::array<int, 2>, 4> offsets = {{{0, -1}, {1, 0}, {0, 1}, {-1, 0}}};
      double worst_angle = 0.0;
      double worst_ratio = 1.0;
      for (std::size_t f = 0; f < 4; ++f) {
        const Vec3& a = p[f];
        const Vec3& b = p[(f + 1) % 4];
        const Vec3 face = b - a;
        const Vec3 normal{face.y, -face.x, 0.0};
        const int nbr_i = i + offsets[f][0];
        const int nbr_j = j + offsets[f][1];
        Vec3 d;
        if (nbr_i >= 0 && nbr_i < ci && nbr_j >= 0 && nbr_j < cj) {
          const std::size_t nc = cell(nbr_i, nbr_j);
          d = centroid[nc] - centroid[c];
          const double an = std::fabs(q.area[nc]);
          const double ac = std::fabs(q.area[c]);
          if (an > 0.0 && ac > 0.0) {
            worst_ratio = std::max({worst_ratio, an / ac, ac / an});
          }
        } else {
          d = 0.5 * (a + b) - centroid[c];
        }
        worst_angle = std::max(worst_angle, angle_between(normal, d));
      }
      q.non_orthogonality[c] = worst_angle;
      q.area_ratio[c] = worst_ratio;
    }
  }
  return q;
}

}  // namespace og
