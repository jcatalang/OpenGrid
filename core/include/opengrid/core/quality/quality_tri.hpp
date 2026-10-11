// Cell quality metrics of planar triangle meshes (counter-clockwise about +z), with the same
// names as the quadrilateral metrics in quality2d.hpp:
//  scaled_jacobian    min over corners of sin(angle) / sin(60 deg) (1 for equilateral, <= 0
//                     degenerate or inverted)
//  skewness           equiangle skewness max((theta_max - 60) / 120, (60 - theta_min) / 60)
//  aspect_ratio       circumradius / (2 inradius) (1 for equilateral)
//  min_angle, max_angle, non_orthogonality (FV, degrees), area_ratio, area
// References: P. M. Knupp, SIAM J. Sci. Comput. 23 (2001); J. H. Ferziger, M. Peric (2002).
#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include "opengrid/core/quality/quality2d.hpp"

namespace og {

Quality2D quality_triangles(const std::vector<std::array<double, 2>>& points,
                            const std::vector<std::array<std::int64_t, 3>>& triangles);

}  // namespace og
