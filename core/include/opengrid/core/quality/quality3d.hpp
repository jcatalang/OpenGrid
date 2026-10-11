// Cell quality metrics of structured 3D blocks (hexahedra, cells (i, j, k) -> (i+1, j+1, k+1),
// i fastest):
//  scaled_jacobian    min over the 8 corners of det(e1, e2, e3) / (|e1| |e2| |e3|) with the three
//                     edges leaving the corner in right-handed order; <= 0: degenerate/inverted
//  skewness           largest equiangle skewness of the six faces, max((t_max - 90) / 90,
//                     (90 - t_min) / 90)
//  aspect_ratio       largest / smallest mean edge length of the i, j and k edge families
//  min_angle,
//  max_angle          smallest / largest face corner angle (degrees)
//  non_orthogonality  finite-volume non-orthogonality (degrees): max over the faces of the angle
//                     between the face area vector and the centroid-to-neighbour-centroid vector
//                     (centroid to face centre on the block boundary)
//  volume_ratio       max over face neighbours of max(V / V_n, V_n / V)
//  volume             exact volume of the trilinear cell (2 x 2 x 2 Gauss quadrature of det J,
//                     which is of degree 2 in each variable)
// References: P. M. Knupp, SIAM J. Sci. Comput. 23 (2001); Ferziger & Peric (2002).
#pragma once

#include <vector>

#include "opengrid/core/structured/grid3d.hpp"

namespace og {

struct Quality3D {
  std::vector<double> scaled_jacobian;
  std::vector<double> skewness;
  std::vector<double> aspect_ratio;
  std::vector<double> min_angle;
  std::vector<double> max_angle;
  std::vector<double> non_orthogonality;
  std::vector<double> volume_ratio;
  std::vector<double> volume;
};

Quality3D quality_3d(const Grid3D& grid);

}  // namespace og
