// A chain of curve pieces traversed in sequence: piece k runs on curve_k from t0_k to t1_k
// (t1 < t0 traverses the curve backwards). Consecutive pieces must meet within a tolerance.
//
// Two parametrisations of the whole chain:
//  * arc length s in [0, L];
//  * a chain parameter u in [0, 1] in which piece k occupies [l_k / L, l_{k+1} / L] (its
//    cumulative length fractions) and the curve parameter varies linearly inside the piece.
#pragma once

#include <memory>
#include <utility>
#include <vector>

#include "opengrid/core/geom/curve.hpp"

namespace og {

struct ChainSegment {
  std::shared_ptr<const Curve> curve;
  double t0 = 0.0;
  double t1 = 1.0;
};

class CurveChain {
 public:
  // Throws std::invalid_argument (naming the segment index) for an empty chain, a null curve,
  // parameters outside [0, 1], t0 == t1, or a gap larger than `tolerance` between pieces.
  CurveChain(std::vector<ChainSegment> segments, double tolerance);

  [[nodiscard]] double length() const {
    return cumulative_.back();
  }
  [[nodiscard]] std::size_t size() const {
    return segments_.size();
  }
  [[nodiscard]] Vec3 start() const;
  [[nodiscard]] Vec3 end() const;

  // Point at arc length s from the start (clamped to [0, L]).
  [[nodiscard]] Vec3 point_at_length(double s) const;
  // Piece index and curve parameter of the point at arc length s (clamped to [0, L]).
  [[nodiscard]] std::pair<std::size_t, double> locate_length(double s) const;
  // Point at chain parameter u in [0, 1] (clamped).
  [[nodiscard]] Vec3 point_at_param(double u) const;
  // Arc length from the start to chain parameter u.
  [[nodiscard]] double length_at_param(double u) const;

 private:
  struct Piece {
    ChainSegment segment;
    double offset;  // arc length of the curve from t = 0 to t0
    double length;  // arc length of the piece
  };
  [[nodiscard]] std::size_t piece_at_length(double s) const;
  [[nodiscard]] double piece_param(std::size_t k, double u) const;

  std::vector<Piece> segments_;
  std::vector<double> cumulative_;  // arc length at the start of each piece, plus total
};

}  // namespace og
