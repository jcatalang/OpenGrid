// Transfinite interpolation (TFI) for structured 2D domains.
//
// Method: Boolean sum of Lagrange (linear) projectors, x = P_xi + P_eta - P_xi P_eta, with
// blending parameters taken from the normalised arc length of the boundary distributions
// (Soni's form): at node (i, j) the parameters (xi, eta) solve
//     xi = u_b(i) + eta (u_t(i) - u_b(i)),   eta = v_l(j) + xi (v_r(j) - v_l(j)),
// where u_b, u_t, v_l, v_r are the arc-length fractions of the bottom, top, left and right
// boundary nodes. TFI with these parameters reproduces the boundaries exactly, and is exact
// for grids such as an annular sector whose opposite boundaries carry the same distribution.
// Hermite TFI replaces the eta projector by cubic Hermite interpolation between bottom and top
// with prescribed derivatives dx/deta: normal to the boundary (in the domain plane) on the
// selected walls, scaled by strength * |top - bottom|, and the chord (top - bottom) elsewhere.
// References: J. F. Thompson, Z. U. A. Warsi, C. W. Mastin, Numerical Grid Generation (1985),
// ch. VIII; B. K. Soni, "Two- and three-dimensional grid generation for internal flow
// applications of CFD", AIAA 85-1526 (1985); L. E. Eriksson, AIAA J. 20 (1982) 1313-1320.
#pragma once

#include "opengrid/core/structured/grid2d.hpp"

namespace og {

// Throws std::invalid_argument for fewer than 2 points per direction, inconsistent sizes, or
// corners that do not coincide (tolerance relative to the boundary size).
Grid2D tfi_linear(const Boundary2D& boundary);

struct HermiteOptions {
  bool bottom = true;     // orthogonal at j = 0
  bool top = false;       // orthogonal at j = nj - 1
  double strength = 1.0;  // derivative magnitude relative to |top - bottom|
};

Grid2D tfi_hermite(const Boundary2D& boundary, const HermiteOptions& options);

}  // namespace og
