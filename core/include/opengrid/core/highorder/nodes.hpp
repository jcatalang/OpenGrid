// One-dimensional node sets and Lagrange interpolation for high-order elements.
//
// Method: on the reference interval [-1, 1],
//  * equispaced: x_k = -1 + 2 k / (n - 1) (the node sets of Lagrange cells in VTK, Gmsh and
//    CGNS files);
//  * Gauss-Lobatto-Legendre (GLL): the end points and the roots of P'_{n-1} (Newton from the
//    Chebyshev-Gauss-Lobatto points), weights 2 / (n (n - 1) P_{n-1}(x)^2); exact for degree
//    2 n - 3; the nodes of spectral-element (continuous Galerkin) methods and of DGSEM with
//    boundary points (shared nodes on element faces);
//  * Gauss-Legendre (GL): the roots of P_n (Newton from cos(pi (k + 3/4) / (n + 1/2))),
//    weights 2 / ((1 - x^2) P'_n(x)^2); exact for degree 2 n - 1; interior points only (DG
//    solution and quadrature points).
// Interpolation and differentiation matrices use the barycentric form of the Lagrange
// polynomials (stable for these point sets).
// References: G. E. Karniadakis, S. J. Sherwin, Spectral/hp Element Methods for CFD (2005),
// app. B; J. S. Hesthaven, T. Warburton, Nodal Discontinuous Galerkin Methods (2008), ch. 3;
// J.-P. Berrut, L. N. Trefethen, "Barycentric Lagrange interpolation", SIAM Rev. 46 (2004)
// 501-517.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace og {

enum class NodeSet : std::uint8_t { Equispaced, GaussLobatto, GaussLegendre };

// Parses "equispaced", "gll" / "gauss_lobatto" or "gauss" / "gauss_legendre"; throws
// std::invalid_argument otherwise.
NodeSet node_set_from(const std::string& name);

// `count` points on [-1, 1], increasing (count >= 1; GLL needs count >= 2).
std::vector<double> nodes_1d(NodeSet set, int count);

// Quadrature weights of the GLL or GL points (sum 2); throws for equispaced points.
std::vector<double> weights_1d(NodeSet set, int count);

// Row-major (to.size() x from.size()): values of the Lagrange polynomials through `from` at
// `to` (interpolation of nodal values).
std::vector<double> interpolation_matrix(const std::vector<double>& from,
                                         const std::vector<double>& to);

// Row-major (to.size() x from.size()): derivatives of the Lagrange polynomials through `from`
// at `to`.
std::vector<double> differentiation_matrix(const std::vector<double>& from,
                                           const std::vector<double>& to);

}  // namespace og
