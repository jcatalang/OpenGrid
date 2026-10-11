// Elliptic smoothing of a structured grid on a parametric surface (boundary nodes fixed).
//
// Method: the 2D Poisson system written with the 3D tangent vectors (alpha = |r_eta|^2,
// beta = r_xi . r_eta, gamma = |r_xi|^2), Laplace or Thomas-Middlecoff control functions computed
// from the boundary distributions with 3D dot products, point SOR; every updated node is
// projected back onto the surface (Newton from its previous parameters), so the grid stays on
// the surface. Reference: J. F. Thompson, B. K. Soni, N. P. Weatherill (eds.), Handbook of Grid
// Generation (1999), ch. 9 (surface grid generation by projection).
#pragma once

#include "opengrid/core/geom/surface.hpp"
#include "opengrid/core/structured/elliptic.hpp"

namespace og {

// Uses options.control (Laplace or ThomasMiddlecoff), omega, max_iterations, tolerance; wall
// control and attractions are not supported on surfaces (std::invalid_argument).
EllipticResult smooth_on_surface(const Surface& surface, const Grid2D& initial,
                                 const EllipticOptions& options,
                                 SmoothingProgress* progress = nullptr);

}  // namespace og
