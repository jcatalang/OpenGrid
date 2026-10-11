#include "opengrid/core/geom/curve_chain.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace og {

CurveChain::CurveChain(std::vector<ChainSegment> segments, double tolerance) {
  if (segments.empty()) {
    throw std::invalid_argument("curve chain: at least one segment is required");
  }
  cumulative_.push_back(0.0);
  for (std::size_t k = 0; k < segments.size(); ++k) {
    const ChainSegment& seg = segments[k];
    const std::string where = "segments[" + std::to_string(k) + "]";
    if (!seg.curve) {
      throw std::invalid_argument("curve chain: " + where + " has no curve");
    }
    if (!(seg.t0 >= 0.0 && seg.t0 <= 1.0 && seg.t1 >= 0.0 && seg.t1 <= 1.0)) {
      throw std::invalid_argument("curve chain: " + where + " parameters must lie in [0, 1]");
    }
    if (seg.t0 == seg.t1) {
      throw std::invalid_argument("curve chain: " + where + " has t0 == t1");
    }
    if (k > 0) {
      const ChainSegment& prev = segments[k - 1];
      const double gap = distance(prev.curve->eval(prev.t1), seg.curve->eval(seg.t0));
      if (gap > tolerance) {
        throw std::invalid_argument("curve chain: " + where + " does not start where segments[" +
                                    std::to_string(k - 1) + "] ends (gap " + std::to_string(gap) +
                                    ")");
      }
    }
    const double offset = seg.curve->length_between(0.0, seg.t0);
    const double length = std::fabs(seg.curve->length_between(seg.t0, seg.t1));
    segments_.push_back({seg, offset, length});
    cumulative_.push_back(cumulative_.back() + length);
  }
}

Vec3 CurveChain::start() const {
  const ChainSegment& s = segments_.front().segment;
  return s.curve->eval(s.t0);
}

Vec3 CurveChain::end() const {
  const ChainSegment& s = segments_.back().segment;
  return s.curve->eval(s.t1);
}

std::size_t CurveChain::piece_at_length(double s) const {
  const auto it = std::upper_bound(cumulative_.begin(), cumulative_.end(), s);
  auto k = static_cast<std::size_t>(it - cumulative_.begin());
  if (k == 0) {
    return 0;
  }
  k -= 1;
  return std::min(k, segments_.size() - 1);
}

Vec3 CurveChain::point_at_length(double s) const {
  s = std::clamp(s, 0.0, length());
  const std::size_t k = piece_at_length(s);
  const Piece& piece = segments_[k];
  const double local = s - cumulative_[k];
  if (piece.segment.t1 > piece.segment.t0) {
    return piece.segment.curve->eval(piece.segment.curve->t_at_length(piece.offset + local));
  }
  return piece.segment.curve->eval(piece.segment.curve->t_at_length(piece.offset - local));
}

double CurveChain::piece_param(std::size_t k, double u) const {
  const Piece& piece = segments_[k];
  const double u0 = cumulative_[k] / length();
  const double u1 = cumulative_[k + 1] / length();
  double frac = 0.0;
  if (u1 > u0) {
    frac = std::clamp((u - u0) / (u1 - u0), 0.0, 1.0);
  }
  return piece.segment.t0 + frac * (piece.segment.t1 - piece.segment.t0);
}

Vec3 CurveChain::point_at_param(double u) const {
  u = std::clamp(u, 0.0, 1.0);
  const std::size_t k = piece_at_length(u * length());
  return segments_[k].segment.curve->eval(piece_param(k, u));
}

double CurveChain::length_at_param(double u) const {
  u = std::clamp(u, 0.0, 1.0);
  const std::size_t k = piece_at_length(u * length());
  const ChainSegment& seg = segments_[k].segment;
  return cumulative_[k] + std::fabs(seg.curve->length_between(seg.t0, piece_param(k, u)));
}

}  // namespace og
