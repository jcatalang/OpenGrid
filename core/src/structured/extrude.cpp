#include "opengrid/core/structured/extrude.hpp"

#include <cmath>
#include <stdexcept>

namespace og {
namespace {

void check(const Grid2D& base, const std::vector<double>& s) {
  if (base.ni < 2 || base.nj < 2) {
    throw std::invalid_argument("extrusion: the base grid needs at least 2 x 2 points");
  }
  if (s.size() < 2 || s.front() != 0.0 || s.back() != 1.0) {
    throw std::invalid_argument("extrusion: layer fractions must run from 0 to 1");
  }
  for (std::size_t k = 1; k < s.size(); ++k) {
    if (!(s[k] > s[k - 1])) {
      throw std::invalid_argument("extrusion: layer fractions must increase");
    }
  }
}

// Sum over the base cells of (x_i x x_j) . direction(cell centre).
template <class F>
double orientation(const Grid2D& base, F&& direction) {
  double sum = 0.0;
  for (int j = 0; j + 1 < base.nj; ++j) {
    for (int i = 0; i + 1 < base.ni; ++i) {
      const Vec3 a = base.at(i + 1, j + 1) - base.at(i, j);
      const Vec3 b = base.at(i, j + 1) - base.at(i + 1, j);
      const Vec3 centre =
          0.25 * (base.at(i, j) + base.at(i + 1, j) + base.at(i, j + 1) + base.at(i + 1, j + 1));
      sum += dot(cross(a, b), direction(centre));
    }
  }
  return sum;
}

Grid2D mirror_i(const Grid2D& base) {
  Grid2D out(base.ni, base.nj);
  for (int j = 0; j < base.nj; ++j) {
    for (int i = 0; i < base.ni; ++i) {
      out.at(i, j) = base.at(base.ni - 1 - i, j);
    }
  }
  return out;
}

template <class F>
Extruded build(const Grid2D& base, const std::vector<double>& s, bool mirror, F&& place) {
  Extruded out;
  out.mirrored = mirror;
  Grid2D b = base;
  if (mirror) {
    b = mirror_i(base);
  }
  const auto nk = static_cast<int>(s.size());
  out.grid = Grid3D(b.ni, b.nj, nk);
  for (int k = 0; k < nk; ++k) {
    for (int j = 0; j < b.nj; ++j) {
      for (int i = 0; i < b.ni; ++i) {
        out.grid.at(i, j, k) = place(b, i, j, s[static_cast<std::size_t>(k)]);
      }
    }
  }
  return out;
}

}  // namespace

Extruded extrude_translate(const Grid2D& base, const Vec3& vector,
                           const std::vector<double>& fractions) {
  check(base, fractions);
  if (norm(vector) == 0.0) {
    throw std::invalid_argument("extrusion: zero translation vector");
  }
  const bool mirror = orientation(base, [&](const Vec3&) { return vector; }) < 0.0;
  return build(base, fractions, mirror,
               [&](const Grid2D& b, int i, int j, double s) { return b.at(i, j) + s * vector; });
}

Extruded extrude_rotate(const Grid2D& base, const Vec3& axis_point, const Vec3& axis, double angle,
                        const std::vector<double>& fractions) {
  check(base, fractions);
  if (norm(axis) == 0.0 || angle == 0.0) {
    throw std::invalid_argument("extrusion: zero rotation axis or angle");
  }
  const Vec3 n = axis / norm(axis);
  auto rotate = [&](const Vec3& p, double a) {
    const Vec3 r = p - axis_point;
    const double c = std::cos(a);
    const double si = std::sin(a);
    return axis_point + c * r + si * cross(n, r) + (1.0 - c) * dot(n, r) * n;
  };
  const bool mirror =
      orientation(base, [&](const Vec3& p) { return angle * cross(n, p - axis_point); }) < 0.0;
  return build(base, fractions, mirror, [&](const Grid2D& b, int i, int j, double s) {
    return rotate(b.at(i, j), s * angle);
  });
}

Extruded extrude_normal(const Grid2D& base, double distance, const std::vector<double>& fractions) {
  check(base, fractions);
  if (distance == 0.0) {
    throw std::invalid_argument("extrusion: zero distance");
  }
  const bool mirror = distance < 0.0;
  Grid2D b = base;
  if (mirror) {
    b = mirror_i(base);
  }
  // Node normals of the (possibly mirrored) base: sum of adjacent cell area vectors.
  std::vector<Vec3> normals(b.x.size());
  for (int j = 0; j + 1 < b.nj; ++j) {
    for (int i = 0; i + 1 < b.ni; ++i) {
      const Vec3 area = cross(b.at(i + 1, j + 1) - b.at(i, j), b.at(i, j + 1) - b.at(i + 1, j));
      for (const auto& [a, c] :
           {std::pair{i, j}, std::pair{i + 1, j}, std::pair{i, j + 1}, std::pair{i + 1, j + 1}}) {
        normals[b.index(a, c)] += area;
      }
    }
  }
  for (Vec3& n : normals) {
    const double length = norm(n);
    if (length == 0.0) {
      throw std::invalid_argument("extrusion: degenerate base grid (zero normal)");
    }
    n = n / length;
  }
  const double d = std::fabs(distance);
  Extruded out;
  out.mirrored = mirror;
  const auto nk = static_cast<int>(fractions.size());
  out.grid = Grid3D(b.ni, b.nj, nk);
  for (int k = 0; k < nk; ++k) {
    for (int j = 0; j < b.nj; ++j) {
      for (int i = 0; i < b.ni; ++i) {
        out.grid.at(i, j, k) =
            b.at(i, j) + fractions[static_cast<std::size_t>(k)] * d * normals[b.index(i, j)];
      }
    }
  }
  return out;
}

}  // namespace og
