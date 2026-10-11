// Structured 3D grid of ni x nj x nk points stored with i fastest, then j: x[(k nj + j) ni + i].
// Blocks are right-handed: (x_i x x_j) . x_k > 0.
#pragma once

#include <cstddef>
#include <vector>

#include "opengrid/core/geom/vec3.hpp"
#include "opengrid/core/structured/grid2d.hpp"

namespace og {

struct Grid3D {
  int ni = 0;
  int nj = 0;
  int nk = 0;
  std::vector<Vec3> x;

  Grid3D() = default;
  Grid3D(int ni_, int nj_, int nk_)
      : ni(ni_),
        nj(nj_),
        nk(nk_),
        x(static_cast<std::size_t>(ni_) * static_cast<std::size_t>(nj_) *
          static_cast<std::size_t>(nk_)) {}

  [[nodiscard]] std::size_t index(int i, int j, int k) const {
    return (static_cast<std::size_t>(k) * static_cast<std::size_t>(nj) +
            static_cast<std::size_t>(j)) *
               static_cast<std::size_t>(ni) +
           static_cast<std::size_t>(i);
  }
  Vec3& at(int i, int j, int k) {
    return x[index(i, j, k)];
  }
  [[nodiscard]] const Vec3& at(int i, int j, int k) const {
    return x[index(i, j, k)];
  }
};

// The six faces of a block, each a Grid2D with its first index running along the first listed
// block index and its second along the second:
//   imin, imax: (j, k)   jmin, jmax: (i, k)   kmin, kmax: (i, j)
// Shared edges and corners must coincide.
struct BlockFaces {
  Grid2D imin;
  Grid2D imax;
  Grid2D jmin;
  Grid2D jmax;
  Grid2D kmin;
  Grid2D kmax;
};

}  // namespace og
