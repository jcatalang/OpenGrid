// Structured 2D grid of ni x nj points stored with i fastest: x[j * ni + i].
#pragma once

#include <cstddef>
#include <vector>

#include "opengrid/core/geom/vec3.hpp"

namespace og {

struct Grid2D {
  int ni = 0;
  int nj = 0;
  std::vector<Vec3> x;

  Grid2D() = default;
  Grid2D(int ni_, int nj_)
      : ni(ni_), nj(nj_), x(static_cast<std::size_t>(ni_) * static_cast<std::size_t>(nj_)) {}

  [[nodiscard]] std::size_t index(int i, int j) const {
    return static_cast<std::size_t>(j) * static_cast<std::size_t>(ni) + static_cast<std::size_t>(i);
  }
  Vec3& at(int i, int j) {
    return x[index(i, j)];
  }
  [[nodiscard]] const Vec3& at(int i, int j) const {
    return x[index(i, j)];
  }
};

// The four boundaries of a structured domain, all ordered by increasing index:
// bottom (j = 0) and top (j = nj - 1) have ni points, left (i = 0) and right (i = ni - 1)
// have nj points; corners are shared (bottom[0] == left[0], bottom[ni-1] == right[0],
// top[0] == left[nj-1], top[ni-1] == right[nj-1]).
struct Boundary2D {
  std::vector<Vec3> bottom;
  std::vector<Vec3> right;
  std::vector<Vec3> top;
  std::vector<Vec3> left;
};

}  // namespace og
