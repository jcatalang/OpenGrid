// Structured domains on a parametric surface: the boundary nodes are projected onto the surface
// (their parameters (u, v)), the interior is the linear TFI of those parameters (Soni's
// arc-length blending, as tfi_linear) and is evaluated on the surface; the boundary nodes are
// kept exactly. The domain must not cross a parameter seam of a closed surface.
#pragma once

#include "opengrid/core/geom/surface.hpp"
#include "opengrid/core/structured/grid2d.hpp"

namespace og {

struct SurfaceTfi {
  Grid2D grid;
  double max_boundary_distance = 0.0;  // largest distance of a boundary node from the surface
};

SurfaceTfi tfi_on_surface(const Surface& surface, const Boundary2D& boundary);

}  // namespace og
