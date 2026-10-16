// Anderson acceleration of a fixed-point iteration u <- G(u).
//
// Method: type-II Anderson mixing without damping (H. F. Walker, P. Ni, "Anderson acceleration
// for fixed-point iterations", SIAM J. Numer. Anal. 49 (2011) 1715-1735). With the last m
// differences of the maps dG_k = G(u_k) - G(u_{k-1}) and of the residuals dF_k = f_k - f_{k-1}
// (f = G(u) - u), the next iterate is G(u) - dG gamma, where gamma minimises |f - dF gamma|_2.
// For a linear map this is equivalent to GMRES on (I - G) u = b; for the relaxation sweeps of
// elliptic grid generation it removes the slowly decaying smooth error modes. The least-squares
// problem is solved through its normal equations (Gram matrix updated by one row per step, a
// relative ridge of 1e-10, Gaussian elimination with partial pivoting); a singular or
// non-finite system clears the history. Dot products sum fixed-size chunks in index order, so
// the result does not depend on the number of threads.
//
// Safeguard (guarded_step): undamped Anderson mixing is not globally convergent (Walker & Ni
// §5; it can diverge where the plain iteration contracts, e.g. elliptic smoothing with stiff
// wall sources at wall spacings of 1e-5 body lengths). The history restarts when the residual
// grows above 4x its best value since the last restart; the iterate of lowest residual is
// kept, and when the residual exceeds kBlowup times that lowest value, is not finite, or has
// not fallen by 10 % for kStall steps (mixing that stagnates while the plain iteration would
// converge, e.g. with fixed wall layers), the iteration restarts from it with an empty history.
// After kMaxRestarts such restarts acceleration is switched off and the plain iteration
// continues, which converges where its map contracts.
#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace og {

class AndersonAccelerator {
 public:
  AndersonAccelerator(std::size_t size, int depth);

  // u: the iterate before the map; g: G(u), replaced by the accelerated iterate. Returns false
  // when no extrapolation was made (first step after a reset): g is then left unchanged.
  bool step(const std::vector<double>& u, std::vector<double>& g);

  // Forget the history (after a residual jump, or when the map changes).
  void reset();

  // The safeguarded step for an iteration whose residual (the size of its last update) is
  // ``residual``: g is G(u), replaced by the accelerated iterate or, after a blow-up, by the
  // best iterate so far. Returns true when g was replaced (the caller then unpacks it).
  bool guarded_step(const std::vector<double>& u, std::vector<double>& g, double residual);

  // Blow-ups that restarted from the best iterate, and whether acceleration is still on.
  [[nodiscard]] int restarts() const {
    return restarts_;
  }
  [[nodiscard]] bool enabled() const {
    return enabled_;
  }
  // The warning for a run with restarts ("" without).
  [[nodiscard]] std::string restart_note() const;

  static constexpr double kBlowup = 100.0;
  static constexpr int kStall = 200;
  static constexpr int kMaxRestarts = 2;

 private:
  bool solve(std::vector<double>& gamma) const;

  std::size_t size_;
  std::size_t depth_;
  std::size_t count_ = 0;  // columns in the history (<= depth)
  std::size_t next_ = 0;   // ring position of the next column
  bool primed_ = false;    // last_f_ / last_g_ hold the previous step
  // The history rings in single precision: they are differences between consecutive
  // iterates and residuals (small, and only ever combined linearly into the next iterate), and
  // reading them is what a step costs on a large grid - 2 x depth x size values per step,
  // halved by the float (11b/12d: 34 % of the 1M-cell group iteration). Everything else is
  // double; the chunked sums keep the result independent of the thread count.
  std::vector<float> df_;     // depth x size, ring of residual differences
  std::vector<float> dg_;     // depth x size, ring of map differences
  std::vector<double> gram_;  // depth x depth, by ring slot
  std::vector<double> last_f_;
  std::vector<double> last_g_;
  std::vector<double> f_;
  bool started_ = false;  // a residual has been seen
  double best_ = 0.0;     // best residual since the last history restart
  double lowest_ = 0.0;   // lowest residual of the run
  double mark_ = 0.0;     // residual of the last 10 % improvement
  int since_mark_ = 0;
  std::vector<double> best_state_;
  int restarts_ = 0;
  bool enabled_ = true;
};

}  // namespace og
