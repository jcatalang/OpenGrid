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
#pragma once

#include <cstddef>
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

 private:
  [[nodiscard]] double dot(const double* a, const double* b) const;
  bool solve(std::vector<double>& gamma) const;

  std::size_t size_;
  std::size_t depth_;
  std::size_t count_ = 0;     // columns in the history (<= depth)
  std::size_t next_ = 0;      // ring position of the next column
  bool primed_ = false;       // last_f_ / last_g_ hold the previous step
  std::vector<double> df_;    // depth x size, ring of residual differences
  std::vector<double> dg_;    // depth x size, ring of map differences
  std::vector<double> gram_;  // depth x depth, by ring slot
  std::vector<double> last_f_;
  std::vector<double> last_g_;
  std::vector<double> f_;
};

}  // namespace og
