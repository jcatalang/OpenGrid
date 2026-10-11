// Cell quality metrics of planar unstructured meshes made of triangles and quadrilaterals
// (counter-clockwise about +z). Same names and definitions as quality2d.hpp (quadrilaterals) and
// quality_tri.hpp (triangles); the ideal angle is 60 degrees for triangles, 90 for quads:
//  scaled_jacobian    min over corners of sin(angle), divided by sin(60 deg) for triangles
//  skewness           equiangle skewness max((theta_max - te) / (180 - te), (te - theta_min) / te)
//  aspect_ratio       triangles: circumradius / (2 inradius); quads: ratio of the larger to the
//                     smaller mean length of opposite edges
//  min_angle, max_angle, non_orthogonality (FV, degrees), area_ratio (edge neighbours), area
// Cells are numbered triangles first, then quads.
// References: P. M. Knupp, SIAM J. Sci. Comput. 23 (2001); J. H. Ferziger, M. Peric (2002).
#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include "opengrid/core/quality/quality2d.hpp"

namespace og {

Quality2D quality_cells(const std::vector<std::array<double, 2>>& points,
                        const std::vector<std::array<std::int64_t, 3>>& triangles,
                        const std::vector<std::array<std::int64_t, 4>>& quads);

}  // namespace og
