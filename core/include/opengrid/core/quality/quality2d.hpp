// Cell quality metrics of planar structured 2D grids (cells (i, j) -> (i+1, j+1), stored i
// fastest; the grid is assumed counter-clockwise about +z, i.e. r_xi x r_eta . z > 0).
//
//  scaled_jacobian    min over corners of (e1 x e2).z / (|e1| |e2|), in [-1, 1]; <= 0 means
//                     a degenerate or inverted cell
//  skewness           equiangle skewness max((theta_max - 90) / 90, (90 - theta_min) / 90)
//  aspect_ratio       ratio of the larger to the smaller mean length of opposite edges
//  min_angle,
//  max_angle          included angles in degrees (reflex angles > 180 for non-convex cells)
//  non_orthogonality  finite-volume non-orthogonality in degrees: max over the cell's faces of
//                     the angle between the face normal and the vector joining the centroids
//                     of the adjacent cells (cell centroid to face midpoint on the boundary)
//  area_ratio         max over edge neighbours of max(A / A_n, A_n / A)
// References: P. M. Knupp, "Algebraic mesh quality metrics", SIAM J. Sci. Comput. 23 (2001);
// J. H. Ferziger, M. Peric, Computational Methods for Fluid Dynamics (2002), §8.6.
#pragma once

#include <vector>

#include "opengrid/core/structured/grid2d.hpp"

namespace og {

struct Quality2D {
  std::vector<double> scaled_jacobian;
  std::vector<double> skewness;
  std::vector<double> aspect_ratio;
  std::vector<double> min_angle;
  std::vector<double> max_angle;
  std::vector<double> non_orthogonality;
  std::vector<double> area_ratio;
  std::vector<double> area;
};

Quality2D quality_2d(const Grid2D& grid);

}  // namespace og
