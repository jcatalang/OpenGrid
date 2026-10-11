// Elliptic smoothing of planar structured 2D grids (boundary nodes fixed).
//
// Method: the Poisson system in computational space (xi = i, eta = j),
//   alpha x_xixi - 2 beta x_xieta + gamma x_etaeta + alpha P x_xi + gamma Q x_eta = 0
// (same for y), alpha = |r_eta|^2, beta = r_xi . r_eta, gamma = |r_xi|^2, discretised with
// second-order central differences and solved by point SOR or alternating-direction line SOR
// (an eta-line sweep then a xi-line sweep per iteration; source terms implicit along the
// lines) with coefficients frozen per sweep. The control terms are written in the
// alpha/gamma-scaled form of Thomas & Middlecoff, which keeps P and Q dimensionless and O(1)
// even at very small wall spacings. Control functions:
//  * Laplace (Winslow-type smoothing): P = Q = 0.
//  * Thomas-Middlecoff: P from the bottom/top boundary distributions, -(r_xi.r_xixi)/|r_xi|^2,
//    interpolated linearly in eta; Q likewise from left/right, interpolated in xi.
//  * TTM attraction (Thompson, Thames & Mastin): P = -sum a sgn(xi - xi_l) exp(-c |xi - xi_l|)
//    for xi-lines and -sum b sgn(xi - xi_p) exp(-d r_p) for points (Q likewise with eta).
// Wall control (bottom and/or top boundary) uses Steger-Sorenson source terms
// P = p(xi) exp(-a eta'), Q = q(xi) exp(-a eta') (eta' = grid lines from the wall). p and q
// are updated every iteration by direct feedback on the achieved geometry (Hilgenstock):
// q from the relative error of the first-cell height, p from the deviation of the first grid
// line from the wall normal, under-relaxed and clamped. At convergence the requested spacing
// and orthogonality are met. (Sorenson's GRAPE update, which infers p and q from derivative
// estimates, was tried first and left ~6% spacing and ~4 degree angle errors.)
// References: J. F. Thompson, F. C. Thames, C. W. Mastin, J. Comput. Phys. 15 (1974) 299-319;
// P. D. Thomas, J. F. Middlecoff, AIAA J. 18 (1980) 652-656; J. L. Steger, R. L. Sorenson,
// J. Comput. Phys. 33 (1979) 405-410; R. L. Sorenson, NASA TM-81198 (1980); A. Hilgenstock,
// "A fast method for the elliptic generation of three-dimensional grids with full boundary
// control", Numerical Grid Generation in CFD '88, Pineridge (1988) 137-146.
#pragma once

#include <atomic>
#include <mutex>
#include <optional>
#include <string>
#include <vector>

#include "opengrid/core/structured/grid2d.hpp"

namespace og {

enum class ControlType { Laplace, ThomasMiddlecoff, Ttm };
enum class EllipticSolver { PointSor, LineSor };

struct LineAttraction {
  bool xi_line = true;  // attraction to the line xi = position (else eta = position)
  double position = 0.0;
  double amplitude = 0.0;
  double decay = 0.5;
};

struct PointAttraction {
  double xi = 0.0;
  double eta = 0.0;
  double amplitude = 0.0;
  double decay = 0.5;
};

struct WallControl {
  bool enabled = false;
  std::optional<double> spacing;  // first-cell height; default: from the side boundaries
  double decay = 0.1;             // a in exp(-a eta'), per grid line
  double relaxation = 0.1;        // feedback gain of the p, q updates
  double max_source = 100.0;      // clamp on |p|, |q| (reaching it is reported as a warning)
};

struct EllipticOptions {
  ControlType control = ControlType::Laplace;
  std::vector<LineAttraction> lines;
  std::vector<PointAttraction> points;
  WallControl bottom;
  WallControl top;
  EllipticSolver solver = EllipticSolver::LineSor;
  double omega = 1.0;
  int max_iterations = 1000;
  double tolerance = 1e-6;  // on the residual relative to the first iteration (floor: 1e-13 L)
};

// Shared with a monitoring thread: iteration count, residual history and a cancel flag.
class SmoothingProgress {
 public:
  void cancel() {
    cancel_.store(true);
  }
  [[nodiscard]] bool cancelled() const {
    return cancel_.load();
  }
  [[nodiscard]] int iteration() const {
    return iteration_.load();
  }
  [[nodiscard]] std::vector<double> residuals() const {
    const std::scoped_lock lock(mutex_);
    return residuals_;
  }
  void record(double residual) {
    const std::scoped_lock lock(mutex_);
    residuals_.push_back(residual);
    iteration_.store(static_cast<int>(residuals_.size()));
  }

 private:
  std::atomic<bool> cancel_{false};
  std::atomic<int> iteration_{0};
  mutable std::mutex mutex_;
  std::vector<double> residuals_;
};

struct EllipticResult {
  Grid2D grid;
  std::vector<double> residuals;  // RMS node correction per iteration
  int iterations = 0;
  bool converged = false;
  bool cancelled = false;
  std::vector<std::string> warnings;
};

// Throws std::invalid_argument for grids smaller than 3 x 3 or not planar in z = const.
EllipticResult smooth_elliptic(const Grid2D& initial, const EllipticOptions& options,
                               SmoothingProgress* progress = nullptr);

}  // namespace og
