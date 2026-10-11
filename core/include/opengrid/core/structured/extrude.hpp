// Blocks by extrusion of a structured surface grid (the k = 0 face), with layer fractions
// s_0 = 0 < s_1 < ... < s_{nk-1} = 1 (from a connector-style distribution):
//  * translate: x_k = x + s_k d;
//  * rotate:    x_k = x rotated about the axis (point, direction) by s_k * angle (Rodrigues);
//  * normal:    x_k = x + s_k * distance * n(x), with n the unit node normal of the surface grid
//               (normalised sum of the adjacent cells' area vectors; the base grid's
//               orientation, x_i x x_j, sets the side).
// The block is right-handed when the extrusion goes to the side of x_i x x_j; otherwise the i
// direction is reversed so that it is (the base face then appears mirrored in i).
#pragma once

#include <vector>

#include "opengrid/core/structured/grid3d.hpp"

namespace og {

struct Extruded {
  Grid3D grid;
  bool mirrored = false;  // i reversed to make the block right-handed
};

// Throw std::invalid_argument for invalid fractions, a zero vector / axis / distance.
Extruded extrude_translate(const Grid2D& base, const Vec3& vector,
                           const std::vector<double>& fractions);
Extruded extrude_rotate(const Grid2D& base, const Vec3& axis_point, const Vec3& axis, double angle,
                        const std::vector<double>& fractions);
Extruded extrude_normal(const Grid2D& base, double distance, const std::vector<double>& fractions);

}  // namespace og
