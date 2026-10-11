#include "opengrid/core/structured/surface_tfi.hpp"

#include <algorithm>

#include "opengrid/core/structured/tfi.hpp"

namespace og {

SurfaceTfi tfi_on_surface(const Surface& surface, const Boundary2D& boundary) {
  SurfaceTfi out;
  auto to_uv = [&](const std::vector<Vec3>& side) {
    std::vector<Vec3> uv;
    uv.reserve(side.size());
    for (const Vec3& p : side) {
      const SurfaceProjection q = surface.project(p);
      out.max_boundary_distance = std::max(out.max_boundary_distance, q.distance);
      uv.emplace_back(q.u, q.v, 0.0);
    }
    return uv;
  };
  Boundary2D uv;
  uv.bottom = to_uv(boundary.bottom);
  uv.right = to_uv(boundary.right);
  uv.top = to_uv(boundary.top);
  uv.left = to_uv(boundary.left);
  // Shared corners must have identical parameters for tfi_linear.
  uv.left.front() = uv.bottom.front();
  uv.right.front() = uv.bottom.back();
  uv.left.back() = uv.top.front();
  uv.right.back() = uv.top.back();
  const Grid2D param = tfi_linear(uv);
  out.grid = Grid2D(param.ni, param.nj);
  for (std::size_t k = 0; k < param.x.size(); ++k) {
    out.grid.x[k] = surface.eval(param.x[k].x, param.x[k].y);
  }
  const int ni = param.ni;
  const int nj = param.nj;
  for (int i = 0; i < ni; ++i) {
    out.grid.at(i, 0) = boundary.bottom[static_cast<std::size_t>(i)];
    out.grid.at(i, nj - 1) = boundary.top[static_cast<std::size_t>(i)];
  }
  for (int j = 0; j < nj; ++j) {
    out.grid.at(0, j) = boundary.left[static_cast<std::size_t>(j)];
    out.grid.at(ni - 1, j) = boundary.right[static_cast<std::size_t>(j)];
  }
  return out;
}

}  // namespace og
