#include "opengrid/core/unstructured/smoothing.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <utility>

#include "opengrid/core/unstructured/advancing_front.hpp"

namespace og {
namespace {

using Index = std::int64_t;

std::size_t u(Index i) {
  return static_cast<std::size_t>(i);
}

// Cells as (corners, count); triangles first, then quads.
struct CellRef {
  std::array<Index, 4> v{};
  int n = 3;
};

class Topology {
 public:
  Topology(std::size_t n_points, const std::vector<std::array<Index, 3>>& triangles,
           const std::vector<std::array<Index, 4>>& quads) {
    cells_.reserve(triangles.size() + quads.size());
    for (const auto& t : triangles) {
      cells_.push_back({{t[0], t[1], t[2], -1}, 3});
    }
    for (const auto& q : quads) {
      cells_.push_back({q, 4});
    }
    // Node -> incident cells and edge neighbours (CSR).
    std::vector<std::vector<Index>> incident(n_points);
    std::vector<std::vector<Index>> neighbours(n_points);
    for (Index c = 0; std::cmp_less(c, cells_.size()); ++c) {
      const CellRef& cell = cells_[u(c)];
      for (int k = 0; k < cell.n; ++k) {
        const Index v = cell.v[u(k)];
        incident[u(v)].push_back(c);
        neighbours[u(v)].push_back(cell.v[u((k + 1) % cell.n)]);
        neighbours[u(v)].push_back(cell.v[u((k + cell.n - 1) % cell.n)]);
      }
    }
    cell_start_.assign(n_points + 1, 0);
    nb_start_.assign(n_points + 1, 0);
    for (std::size_t v = 0; v < n_points; ++v) {
      auto& list = neighbours[v];
      std::sort(list.begin(), list.end());
      list.erase(std::unique(list.begin(), list.end()), list.end());
      cell_start_[v + 1] = cell_start_[v] + incident[v].size();
      nb_start_[v + 1] = nb_start_[v] + list.size();
    }
    cell_list_.reserve(cell_start_.back());
    nb_list_.reserve(nb_start_.back());
    for (std::size_t v = 0; v < n_points; ++v) {
      cell_list_.insert(cell_list_.end(), incident[v].begin(), incident[v].end());
      nb_list_.insert(nb_list_.end(), neighbours[v].begin(), neighbours[v].end());
    }
  }

  [[nodiscard]] const CellRef& cell(Index c) const {
    return cells_[u(c)];
  }
  [[nodiscard]] std::size_t cell_count() const {
    return cells_.size();
  }
  [[nodiscard]] std::pair<const Index*, const Index*> cells_of(std::size_t v) const {
    return {cell_list_.data() + cell_start_[v], cell_list_.data() + cell_start_[v + 1]};
  }
  [[nodiscard]] std::pair<const Index*, const Index*> neighbours_of(std::size_t v) const {
    return {nb_list_.data() + nb_start_[v], nb_list_.data() + nb_start_[v + 1]};
  }

 private:
  std::vector<CellRef> cells_;
  std::vector<std::size_t> cell_start_;
  std::vector<Index> cell_list_;
  std::vector<std::size_t> nb_start_;
  std::vector<Index> nb_list_;
};

double worst_around(const std::vector<std::array<double, 2>>& points, const Topology& topo,
                    std::size_t v) {
  double worst = 2.0;
  const auto [begin, end] = topo.cells_of(v);
  for (const Index* c = begin; c != end; ++c) {
    const CellRef& cell = topo.cell(*c);
    worst = std::min(worst, cell_quality(points, cell.v, cell.n));
  }
  return worst;
}

std::array<double, 2> proposal(const std::vector<std::array<double, 2>>& points,
                               const Topology& topo, std::size_t v, SmoothingMethod method,
                               bool& ok) {
  const auto [begin, end] = topo.neighbours_of(v);
  const auto& pv = points[v];
  ok = end != begin;
  if (!ok) {
    return pv;
  }
  if (method == SmoothingMethod::laplacian || method == SmoothingMethod::length_weighted) {
    double sx = 0.0;
    double sy = 0.0;
    double sw = 0.0;
    for (const Index* w = begin; w != end; ++w) {
      const auto& pw = points[u(*w)];
      double weight = 1.0;
      if (method == SmoothingMethod::length_weighted) {
        weight = std::hypot(pw[0] - pv[0], pw[1] - pv[1]);
      }
      sx += weight * pw[0];
      sy += weight * pw[1];
      sw += weight;
    }
    if (sw <= 0.0) {
      ok = false;
      return pv;
    }
    return {sx / sw, sy / sw};
  }
  // Angle based: ring of edge neighbours ordered by angle around v.
  std::vector<std::pair<double, Index>> ring;
  for (const Index* w = begin; w != end; ++w) {
    const auto& pw = points[u(*w)];
    ring.emplace_back(std::atan2(pw[1] - pv[1], pw[0] - pv[0]), *w);
  }
  std::sort(ring.begin(), ring.end());
  const std::size_t m = ring.size();
  if (m < 3) {
    ok = false;
    return pv;
  }
  double sx = 0.0;
  double sy = 0.0;
  double count = 0.0;
  for (std::size_t j = 0; j < m; ++j) {
    const auto& pj = points[u(ring[j].second)];
    const auto& prev = points[u(ring[(j + m - 1) % m].second)];
    const auto& next = points[u(ring[(j + 1) % m].second)];
    const double ax = prev[0] - pj[0];
    const double ay = prev[1] - pj[1];
    const double bx = next[0] - pj[0];
    const double by = next[1] - pj[1];
    const double la = std::hypot(ax, ay);
    const double lb = std::hypot(bx, by);
    const double vx = pv[0] - pj[0];
    const double vy = pv[1] - pj[1];
    const double lv = std::hypot(vx, vy);
    if (la == 0.0 || lb == 0.0 || lv == 0.0) {
      continue;
    }
    double dx = ax / la + bx / lb;
    double dy = ay / la + by / lb;
    const double ld = std::hypot(dx, dy);
    if (ld < 1e-12) {
      continue;
    }
    dx /= ld;
    dy /= ld;
    if (dx * vx + dy * vy < 0.0) {
      dx = -dx;
      dy = -dy;
    }
    sx += pj[0] + lv * dx;
    sy += pj[1] + lv * dy;
    count += 1.0;
  }
  if (count == 0.0) {
    ok = false;
    return pv;
  }
  return {sx / count, sy / count};
}

}  // namespace

double cell_quality(const std::vector<std::array<double, 2>>& points,
                    const std::array<std::int64_t, 4>& corners, int n) {
  double worst = 2.0;
  double scale = 1.0;
  if (n == 3) {
    scale = 2.0 / std::numbers::sqrt3;
  }
  for (int k = 0; k < n; ++k) {
    const auto& o = points[u(corners[u(k)])];
    const auto& a = points[u(corners[u((k + 1) % n)])];
    const auto& b = points[u(corners[u((k + n - 1) % n)])];
    const double ux = a[0] - o[0];
    const double uy = a[1] - o[1];
    const double vx = b[0] - o[0];
    const double vy = b[1] - o[1];
    const double lu = std::hypot(ux, uy);
    const double lv = std::hypot(vx, vy);
    if (lu == 0.0 || lv == 0.0) {
      return -1.0;
    }
    worst = std::min(worst, (ux * vy - uy * vx) / (lu * lv) * scale);
  }
  return worst;
}

std::int64_t smooth_mesh(std::vector<std::array<double, 2>>& points,
                         std::vector<std::array<std::int64_t, 3>>& triangles,
                         const std::vector<std::array<std::int64_t, 4>>& quads,
                         const std::vector<std::array<std::int64_t, 2>>& segments,
                         std::size_t n_fixed, const SmoothingOptions& options) {
  if (options.method == SmoothingMethod::none || options.iterations <= 0) {
    return 0;
  }
  std::int64_t moved = 0;
  for (int pass = 0; pass < options.iterations; ++pass) {
    // Edge swaps change the topology; rebuild it every sweep.
    const Topology topo(points.size(), triangles, quads);
    for (std::size_t v = n_fixed; v < points.size(); ++v) {
      bool ok = false;
      const std::array<double, 2> target = proposal(points, topo, v, options.method, ok);
      if (!ok) {
        continue;
      }
      const std::array<double, 2> old = points[v];
      const double before = worst_around(points, topo, v);
      points[v] = target;
      if (worst_around(points, topo, v) < before) {
        points[v] = old;
      } else {
        moved += 1;
      }
    }
    if (options.edge_swap && !triangles.empty()) {
      flip_to_delaunay(points, triangles, segments);
    }
  }
  return moved;
}

std::int64_t untangle(std::vector<std::array<double, 2>>& points,
                      const std::vector<std::array<std::int64_t, 3>>& triangles,
                      const std::vector<std::array<std::int64_t, 4>>& quads, std::size_t n_fixed) {
  // Cheap check first: most meshes have nothing to untangle.
  bool any = false;
  for (const auto& t : triangles) {
    any = any || cell_quality(points, {t[0], t[1], t[2], -1}, 3) <= 0.0;
  }
  for (const auto& q : quads) {
    any = any || cell_quality(points, q, 4) <= 0.0;
  }
  if (!any) {
    return 0;
  }
  const Topology topo(points.size(), triangles, quads);
  auto count_bad = [&]() {
    std::int64_t bad = 0;
    for (Index c = 0; std::cmp_less(c, topo.cell_count()); ++c) {
      const CellRef& cell = topo.cell(c);
      if (cell_quality(points, cell.v, cell.n) <= 0.0) {
        ++bad;
      }
    }
    return bad;
  };
  constexpr std::array<std::array<double, 2>, 8> kDirections = {
      {{1, 0},
       {0.7071067811865476, 0.7071067811865476},
       {0, 1},
       {-0.7071067811865476, 0.7071067811865476},
       {-1, 0},
       {-0.7071067811865476, -0.7071067811865476},
       {0, -1},
       {0.7071067811865476, -0.7071067811865476}}};
  for (int sweep = 0; sweep < 20; ++sweep) {
    if (count_bad() == 0) {
      return 0;
    }
    bool changed = false;
    for (std::size_t v = n_fixed; v < points.size(); ++v) {
      double current = worst_around(points, topo, v);
      if (current > 0.0) {
        continue;
      }
      // Step: a quarter of the mean edge length around the node.
      const auto [begin, end] = topo.neighbours_of(v);
      double mean = 0.0;
      for (const Index* w = begin; w != end; ++w) {
        mean += std::hypot(points[u(*w)][0] - points[v][0], points[u(*w)][1] - points[v][1]);
      }
      if (begin == end) {
        continue;
      }
      // Start from the neighbours' centre if that is better (it repairs a node thrown out of
      // its ring in one step; the compass search alone can stall at a kink of the worst-cell
      // objective).
      bool ok = false;
      const std::array<double, 2> centre =
          proposal(points, topo, v, SmoothingMethod::laplacian, ok);
      if (ok) {
        const std::array<double, 2> origin = points[v];
        points[v] = centre;
        const double value = worst_around(points, topo, v);
        if (value > current) {
          current = value;
          changed = true;
        } else {
          points[v] = origin;
        }
      }
      double step = 0.25 * mean / static_cast<double>(end - begin);
      for (int iteration = 0; iteration < 60 && step > 1e-9 * mean; ++iteration) {
        const std::array<double, 2> origin = points[v];
        std::array<double, 2> best_point = origin;
        double best = current;
        for (const auto& d : kDirections) {
          points[v] = {origin[0] + step * d[0], origin[1] + step * d[1]};
          const double value = worst_around(points, topo, v);
          if (value > best) {
            best = value;
            best_point = points[v];
          }
        }
        points[v] = best_point;
        if (best > current) {
          current = best;
          changed = true;
        } else {
          step *= 0.5;
        }
      }
    }
    if (!changed) {
      break;
    }
  }
  return count_bad();
}

}  // namespace og
