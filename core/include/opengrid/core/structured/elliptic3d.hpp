// Elliptic smoothing of structured 3D blocks (boundary faces fixed).
//
// Method: the 3D Poisson system in computational space (xi^1..3 = i, j, k),
//   sum_mn a_mn r_{mn} + sum_m a_mm P_m r_m = 0,
// where a_mn are the cofactors of the covariant metric g_mn = r_m . r_n (the contravariant
// metric times J^2), discretised with second-order central differences (the source terms in the
// neighbour coefficients, as in the 2D smoother) and solved by point SOR with coefficients frozen
// per node update. Planes k = const are relaxed in zebra order (odd k, then even k; the stencil
// reaches k +/- 1 only), the planes of one colour in parallel, nodes within a plane in order -
// the result does not depend on the thread count. Optional Anderson acceleration of the
// iteration (util/anderson.hpp), as in the 2D smoother. Control functions:
//  * Laplace (Winslow-type smoothing): P_m = 0.
//  * Thomas-Middlecoff (3D extension): on the four faces containing the m-lines,
//    P_m = -(r_m . r_mm) / |r_m|^2, blended into the block by transfinite interpolation in index
//    space.
// Wall control on any of the six faces: decaying sources P_m = p_m exp(-a d) (d = grid lines from
// the wall) with direct feedback (Hilgenstock): the source of the wall-normal index from the
// relative error of the first-cell height, the two tangential sources from the deviation of the
// first grid line from the wall normal, under-relaxed and clamped - the 3D form of the 2D wall
// control (structured/elliptic.hpp).
// References: J. F. Thompson, Z. U. A. Warsi, C. W. Mastin, Numerical Grid Generation (1985),
// ch. VI; P. D. Thomas, J. F. Middlecoff, AIAA J. 18 (1980) 652-656; A. Hilgenstock, Numerical
// Grid Generation in CFD '88 (1988) 137-146.
#pragma once

#include <string>
#include <vector>

#include "opengrid/core/structured/elliptic.hpp"
#include "opengrid/core/structured/grid3d.hpp"

namespace og {

struct WallControl3D {
  int face = 4;  // 0 imin, 1 imax, 2 jmin, 3 jmax, 4 kmin, 5 kmax
  WallControl control;
};

struct Elliptic3DOptions {
  ControlType control = ControlType::Laplace;  // Laplace or ThomasMiddlecoff
  std::vector<WallControl3D> walls;
  double omega = 1.0;
  int max_iterations = 500;
  double tolerance = 1e-6;  // on the residual relative to the first iteration
  int acceleration = 0;     // Anderson acceleration depth (0: plain relaxation)
};

struct Elliptic3DResult {
  Grid3D grid;
  std::vector<double> residuals;
  int iterations = 0;
  bool converged = false;
  bool cancelled = false;
  std::vector<std::string> warnings;
};

// Throws std::invalid_argument for blocks smaller than 3 x 3 x 3 or unsupported controls.
Elliptic3DResult smooth_elliptic_3d(const Grid3D& initial, const Elliptic3DOptions& options,
                                    SmoothingProgress* progress = nullptr);

}  // namespace og
