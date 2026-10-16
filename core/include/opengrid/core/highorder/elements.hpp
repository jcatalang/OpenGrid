// Curved high-order quadrilaterals and hexahedra from a structured grid.
//
// Method: a grid with p x p (x p) cells per element (every interval count a multiple of the
// order p) is the equispaced node set of order-p Lagrange elements: the elliptic smoothing or
// TFI that made the grid placed every high-order node, its boundary nodes lie on the geometry,
// and neighbouring elements share their edge and face nodes. Each element's mapping x(xi) of
// the reference square (cube) [-1, 1]^d is the tensor-product Lagrange interpolant of these
// nodes; other node sets (Gauss-Lobatto-Legendre, Gauss-Legendre) are that mapping evaluated
// at their points. Validity: the Jacobian determinant of the mapping at the (p + 2)^d
// Gauss-Lobatto points of each element (2D: the in-plane determinant, positive for the grid's
// orientation); the element's scaled Jacobian is min / max of it (an element is invalid when
// the minimum is not positive).
// References: S. Sherwin, J. Peiro, "Mesh generation in curvilinear domains using high-order
// elements", Int. J. Numer. Meth. Engng 53 (2002) 207-223; P.-O. Persson, J. Peraire,
// "Curved mesh generation and mesh refinement using Lagrangian solid mechanics", AIAA
// 2009-949 (validity of curved elements through the Jacobian).
#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include "opengrid/core/geom/vec3.hpp"
#include "opengrid/core/highorder/nodes.hpp"

namespace og {

struct HighOrderElements {
  int dim = 2;
  int order = 1;
  std::int64_t count = 0;            // elements
  std::array<std::int64_t, 3> ne{};  // elements per grid direction (i, j, k)
  // Element nodes: count x (order + 1)^dim x 3, tensor-product order (i fastest), elements in
  // grid order (i fastest).
  std::vector<double> xyz;
  std::vector<double> min_jacobian;  // per element
  std::vector<double> max_jacobian;
};

// grid: ni x nj (x nk) points, i fastest (nk = 1 for a 2D grid); throws std::invalid_argument
// unless every interval count is a multiple of `order`.
HighOrderElements high_order_from_grid(const std::vector<Vec3>& grid, int ni, int nj, int nk,
                                       int order, NodeSet set);

}  // namespace og
