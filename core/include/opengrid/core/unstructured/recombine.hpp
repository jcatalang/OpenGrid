// Quad-dominant recombination: pairs of adjacent triangles are merged into quadrilaterals.
//
// Method: a matching in the dual graph of the triangulation (triangles = nodes, shared edges =
// edges), weighted by the quality of the quad each pair would form,
//   eta = max(0, 1 - (2 / pi) max_k |pi / 2 - alpha_k|)
// (J.-F. Remacle et al., "Blossom-Quad", Int. J. Numer. Meth. Eng. 89 (2012) 1102-1119; eta = 1
// for a rectangle, 0 when a corner reaches 0 or 180 degrees). Only strictly convex quads with
// eta >= min_quality are candidates; boundary segments are never removed.
//  1. Greedy matching in order of decreasing eta (ties by triangle indices), which keeps the
//     best quads.
//  2. Augmenting paths: from every unmatched triangle a breadth-first search over alternating
//     paths (candidate edge, matched pair, candidate edge, ...) of at most max_path matched pairs
//     looks for another unmatched triangle; the path is flipped, gaining one quad. This is the
//     augmenting-path step of maximum-cardinality matching (Edmonds), bounded in length, so the
//     result is quad-dominant rather than guaranteed all-quad.
// Blossom-Quad's minimum-cost perfect matching is not used: it needs Edmonds' weighted blossom
// algorithm, and an all-quad mesh is produced instead by subdividing a coarser quad-dominant mesh
// (domain_mesher.hpp). Output is deterministic.
#pragma once

#include <array>
#include <cstdint>
#include <vector>

namespace og {

struct RecombineOptions {
  double min_quality = 0.3;  // eta of candidate quads (corner angles within 90 +- 63 degrees)
  int max_path = 4;          // matched pairs on an augmenting path
};

struct Recombined {
  std::vector<std::array<std::int64_t, 3>> triangles;  // unmatched, in input order
  std::vector<std::array<std::int64_t, 4>> quads;      // counter-clockwise
};

double quad_quality(const std::array<double, 2>& a, const std::array<double, 2>& b,
                    const std::array<double, 2>& c, const std::array<double, 2>& d);

Recombined recombine(const std::vector<std::array<double, 2>>& points,
                     const std::vector<std::array<std::int64_t, 3>>& triangles,
                     const RecombineOptions& options);

}  // namespace og
