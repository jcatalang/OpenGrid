#include "opengrid/core/unstructured/recombine.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <tuple>
#include <utility>

#include "opengrid/core/util/predicates.hpp"

namespace og {
namespace {

using Index = std::int64_t;

std::size_t u(Index i) {
  return static_cast<std::size_t>(i);
}

struct Candidate {
  double quality;
  Index t1;
  Index t2;
  std::array<Index, 4> quad;
};

}  // namespace

double quad_quality(const std::array<double, 2>& a, const std::array<double, 2>& b,
                    const std::array<double, 2>& c, const std::array<double, 2>& d) {
  const std::array<const std::array<double, 2>*, 4> p = {&a, &b, &c, &d};
  double worst = 0.0;
  for (std::size_t k = 0; k < 4; ++k) {
    const auto& o = *p[k];
    const auto& e1 = *p[(k + 1) % 4];
    const auto& e2 = *p[(k + 3) % 4];
    const double ux = e1[0] - o[0];
    const double uy = e1[1] - o[1];
    const double vx = e2[0] - o[0];
    const double vy = e2[1] - o[1];
    const double angle = std::atan2(ux * vy - uy * vx, ux * vx + uy * vy);
    worst = std::max(worst, std::fabs(0.5 * std::numbers::pi - angle));
  }
  return std::max(0.0, 1.0 - 2.0 / std::numbers::pi * worst);
}

Recombined recombine(const std::vector<std::array<double, 2>>& points,
                     const std::vector<std::array<std::int64_t, 3>>& triangles,
                     const RecombineOptions& options) {
  const std::size_t nt = triangles.size();
  // Shared edges from sorted half-edges.
  std::vector<std::tuple<Index, Index, Index>> half;  // (min, max, triangle)
  half.reserve(3 * nt);
  for (std::size_t t = 0; t < nt; ++t) {
    for (std::size_t k = 0; k < 3; ++k) {
      Index a = triangles[t][k];
      Index b = triangles[t][(k + 1) % 3];
      if (b < a) {
        std::swap(a, b);
      }
      half.emplace_back(a, b, static_cast<Index>(t));
    }
  }
  std::sort(half.begin(), half.end());
  std::vector<Candidate> candidates;
  for (std::size_t i = 0; i + 1 < half.size(); ++i) {
    const auto [a, b, t1] = half[i];
    const auto [c, d, t2] = half[i + 1];
    if (a != c || b != d) {
      continue;
    }
    ++i;
    // t1 = (p0, p1, p2) with the shared edge (p1, p2) opposite p0; q opposite in t2.
    const auto& tri1 = triangles[u(t1)];
    std::size_t k1 = 0;
    while (tri1[k1] == a || tri1[k1] == b) {
      ++k1;
    }
    Index q = -1;
    for (const Index v : triangles[u(t2)]) {
      if (v != a && v != b) {
        q = v;
      }
    }
    const std::array<Index, 4> quad = {tri1[k1], tri1[(k1 + 1) % 3], q, tri1[(k1 + 2) % 3]};
    bool convex = true;
    for (std::size_t k = 0; k < 4 && convex; ++k) {
      convex = orient2d(points[u(quad[k])].data(), points[u(quad[(k + 1) % 4])].data(),
                        points[u(quad[(k + 2) % 4])].data()) > 0;
    }
    if (!convex) {
      continue;
    }
    const double eta = quad_quality(points[u(quad[0])], points[u(quad[1])], points[u(quad[2])],
                                    points[u(quad[3])]);
    if (eta >= options.min_quality) {
      candidates.push_back({eta, std::min(t1, t2), std::max(t1, t2), quad});
    }
  }
  std::sort(candidates.begin(), candidates.end(), [](const Candidate& x, const Candidate& y) {
    if (x.quality != y.quality) {
      return x.quality > y.quality;
    }
    if (x.t1 != y.t1) {
      return x.t1 < y.t1;
    }
    return x.t2 < y.t2;
  });
  // Candidate edges per triangle (indices into candidates, best first).
  std::vector<std::vector<Index>> options_of(nt);
  for (Index e = 0; std::cmp_less(e, candidates.size()); ++e) {
    options_of[u(candidates[u(e)].t1)].push_back(e);
    options_of[u(candidates[u(e)].t2)].push_back(e);
  }
  std::vector<Index> mate(nt, -1);  // matched candidate per triangle
  for (Index e = 0; std::cmp_less(e, candidates.size()); ++e) {
    const Candidate& cand = candidates[u(e)];
    if (mate[u(cand.t1)] < 0 && mate[u(cand.t2)] < 0) {
      mate[u(cand.t1)] = e;
      mate[u(cand.t2)] = e;
    }
  }
  auto other = [&](Index e, Index t) {
    if (candidates[u(e)].t1 == t) {
      return candidates[u(e)].t2;
    }
    return candidates[u(e)].t1;
  };
  // Bounded augmenting paths: BFS over "even" triangles reached by an alternating path.
  std::vector<Index> parent_edge(nt, -1);  // candidate edge used to reach an odd triangle
  std::vector<Index> parent_even(nt, -1);  // even triangle we came from
  std::vector<int> depth(nt, 0);
  std::vector<std::uint64_t> seen(nt, 0);
  std::uint64_t stamp = 0;
  bool improved = true;
  while (improved) {
    improved = false;
    for (Index start = 0; std::cmp_less(start, nt); ++start) {
      if (mate[u(start)] >= 0 || options_of[u(start)].empty()) {
        continue;
      }
      stamp += 1;
      std::vector<Index> frontier = {start};
      seen[u(start)] = stamp;
      depth[u(start)] = 0;
      Index found = -1;
      Index found_edge = -1;
      Index found_from = -1;
      while (!frontier.empty() && found < 0) {
        std::vector<Index> next;
        for (const Index even : frontier) {
          for (const Index e : options_of[u(even)]) {
            const Index odd = other(e, even);
            if (seen[u(odd)] == stamp) {
              continue;
            }
            if (mate[u(odd)] < 0) {
              found = odd;
              found_edge = e;
              found_from = even;
              break;
            }
            if (depth[u(even)] >= options.max_path) {
              continue;
            }
            const Index partner = other(mate[u(odd)], odd);
            if (seen[u(partner)] == stamp) {
              continue;
            }
            seen[u(odd)] = stamp;
            seen[u(partner)] = stamp;
            parent_edge[u(odd)] = e;
            parent_even[u(odd)] = even;
            depth[u(partner)] = depth[u(even)] + 1;
            // partner is reached through odd's matched edge.
            parent_even[u(partner)] = odd;
            next.push_back(partner);
          }
          if (found >= 0) {
            break;
          }
        }
        frontier = std::move(next);
      }
      if (found < 0) {
        continue;
      }
      // Flip the path: found -(found_edge)- found_from, then back along the parents.
      Index t = found_from;
      Index edge = found_edge;
      mate[u(found)] = edge;
      while (true) {
        mate[u(t)] = edge;
        if (t == start) {
          break;
        }
        // t (even) was matched to odd = parent_even[t]; odd was reached by parent_edge[odd].
        const Index odd = parent_even[u(t)];
        edge = parent_edge[u(odd)];
        mate[u(odd)] = edge;
        t = parent_even[u(odd)];
      }
      improved = true;
    }
  }
  Recombined out;
  std::vector<bool> emitted(candidates.size(), false);
  for (std::size_t t = 0; t < nt; ++t) {
    const Index e = mate[t];
    if (e < 0) {
      out.triangles.push_back(triangles[t]);
    } else if (!emitted[u(e)]) {
      emitted[u(e)] = true;
      out.quads.push_back(candidates[u(e)].quad);
    }
  }
  return out;
}

}  // namespace og
