#include "opengrid/core/unstructured/optimize3d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numbers>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include "opengrid/core/geom/vec3.hpp"

namespace og {

namespace {

using Index = std::int64_t;

std::size_t u(Index i) {
  return static_cast<std::size_t>(i);
}

Index sz(std::size_t n) {
  return static_cast<Index>(n);
}

Vec3 at(const std::vector<std::array<double, 3>>& points, Index i) {
  const auto& p = points[u(i)];
  return {p[0], p[1], p[2]};
}

// Corners (node, three neighbours, positive triple product for a valid cell) and the factor
// scaling the ideal cell's corner value to 1, per cell size (as CELL_CORNERS in volumes.py).
struct CornerTable {
  std::vector<std::array<int, 4>> corners;
  double factor = 1.0;
};

const CornerTable& corner_table(int n) {
  static const CornerTable tet{{{0, 1, 2, 3}, {1, 2, 0, 3}, {2, 0, 1, 3}, {3, 0, 2, 1}},
                               std::numbers::sqrt2};
  static const CornerTable pyramid{{{0, 1, 3, 4}, {1, 2, 0, 4}, {2, 3, 1, 4}, {3, 0, 2, 4}},
                                   std::numbers::sqrt2};
  static const CornerTable prism{
      {{0, 1, 2, 3}, {1, 2, 0, 4}, {2, 0, 1, 5}, {3, 5, 4, 0}, {4, 3, 5, 1}, {5, 4, 3, 2}},
      2.0 / std::numbers::sqrt3};
  static const CornerTable hex{{{0, 1, 3, 4},
                                {1, 2, 0, 5},
                                {2, 3, 1, 6},
                                {3, 0, 2, 7},
                                {4, 7, 5, 0},
                                {5, 4, 6, 1},
                                {6, 5, 7, 2},
                                {7, 6, 4, 3}},
                               1.0};
  if (n == 4) {
    return tet;
  }
  if (n == 5) {
    return pyramid;
  }
  if (n == 6) {
    return prism;
  }
  if (n == 8) {
    return hex;
  }
  throw std::invalid_argument("cell_quality_3d: a cell has 4, 5, 6 or 8 nodes");
}

// Edges of each cell size (local node pairs).
const std::vector<std::array<int, 2>>& cell_edges(int n) {
  static const std::vector<std::array<int, 2>> tet = {{0, 1}, {1, 2}, {2, 0},
                                                      {0, 3}, {1, 3}, {2, 3}};
  static const std::vector<std::array<int, 2>> pyramid = {{0, 1}, {1, 2}, {2, 3}, {3, 0},
                                                          {0, 4}, {1, 4}, {2, 4}, {3, 4}};
  static const std::vector<std::array<int, 2>> prism = {{0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5},
                                                        {5, 3}, {0, 3}, {1, 4}, {2, 5}};
  static const std::vector<std::array<int, 2>> hex = {{0, 1}, {1, 2}, {2, 3}, {3, 0},
                                                      {4, 5}, {5, 6}, {6, 7}, {7, 4},
                                                      {0, 4}, {1, 5}, {2, 6}, {3, 7}};
  if (n == 4) {
    return tet;
  }
  if (n == 5) {
    return pyramid;
  }
  if (n == 6) {
    return prism;
  }
  return hex;
}

// Triangular faces of each cell size (local nodes; orientation irrelevant here).
const std::vector<std::array<int, 3>>& cell_triangles(int n) {
  static const std::vector<std::array<int, 3>> tet = {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {0, 3, 2}};
  static const std::vector<std::array<int, 3>> pyramid = {
      {0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}};
  static const std::vector<std::array<int, 3>> prism = {{0, 2, 1}, {3, 4, 5}};
  static const std::vector<std::array<int, 3>> none;
  if (n == 4) {
    return tet;
  }
  if (n == 5) {
    return pyramid;
  }
  if (n == 6) {
    return prism;
  }
  return none;
}

struct Cell {
  int n = 4;
  std::array<Index, 8> v{};
  bool alive = true;
  bool swappable = false;
  bool touched = false;  // changed during the current swap pass
};

constexpr double kInvSqrt3 = std::numbers::inv_sqrt3;
constexpr std::array<std::array<double, 3>, 14> kDirections = {
    {{1, 0, 0},
     {-1, 0, 0},
     {0, 1, 0},
     {0, -1, 0},
     {0, 0, 1},
     {0, 0, -1},
     {kInvSqrt3, kInvSqrt3, kInvSqrt3},
     {kInvSqrt3, kInvSqrt3, -kInvSqrt3},
     {kInvSqrt3, -kInvSqrt3, kInvSqrt3},
     {kInvSqrt3, -kInvSqrt3, -kInvSqrt3},
     {-kInvSqrt3, kInvSqrt3, kInvSqrt3},
     {-kInvSqrt3, kInvSqrt3, -kInvSqrt3},
     {-kInvSqrt3, -kInvSqrt3, kInvSqrt3},
     {-kInvSqrt3, -kInvSqrt3, -kInvSqrt3}}};

class Optimizer {
 public:
  Optimizer(const Optimize3DInput& in, const Optimize3DOptions& opt)
      : points_(in.points), movable_(in.movable), opt_(opt) {
    auto add = [&](const auto& list, int n, bool swappable_from_zero, Index first_swappable) {
      for (std::size_t c = 0; c < list.size(); ++c) {
        Cell cell;
        cell.n = n;
        for (int k = 0; k < n; ++k) {
          const Index v = list[c][static_cast<std::size_t>(k)];
          if (v < 0 || std::cmp_greater_equal(v, points_.size())) {
            throw std::invalid_argument("optimize_3d: a cell refers to a missing point");
          }
          cell.v[static_cast<std::size_t>(k)] = v;
        }
        cell.swappable = swappable_from_zero && sz(c) >= first_swappable;
        cells_.push_back(cell);
      }
    };
    add(in.tets, 4, true, in.n_fixed_tets);
    add(in.pyramids, 5, false, 0);
    add(in.prisms, 6, false, 0);
    add(in.hexes, 8, false, 0);
  }

  Optimize3DResult run() {
    Optimize3DResult r;
    std::tie(r.worst_before, r.below_before, std::ignore) = summary();
    for (int pass = 0; pass < opt_.passes; ++pass) {
      Index changes = 0;
      if (opt_.move_nodes) {
        const Index moved = move_nodes();
        r.moved += moved;
        changes += moved;
      }
      if (opt_.swap) {
        const auto [s23, s32] = swap();
        r.swaps_23 += s23;
        r.swaps_32 += s32;
        changes += s23 + s32;
      }
      if (changes == 0) {
        break;
      }
    }
    std::tie(r.worst_after, r.below_after, r.inverted_after) = summary();
    r.points = points_;
    for (const Cell& cell : cells_) {
      if (cell.alive && cell.n == 4) {
        r.tets.push_back({cell.v[0], cell.v[1], cell.v[2], cell.v[3]});
      }
    }
    // The fixed tetrahedra are the first cells and never die, so they stay first and in order.
    return r;
  }

 private:
  [[nodiscard]] double quality(const Cell& cell) const {
    return cell_quality_3d(points_, cell.v.data(), cell.n);
  }

  [[nodiscard]] std::tuple<double, Index, Index> summary() const {
    double worst = 1.0;
    Index below = 0;
    Index inverted = 0;
    for (const Cell& cell : cells_) {
      if (!cell.alive) {
        continue;
      }
      const double q = quality(cell);
      worst = std::min(worst, q);
      if (q < opt_.threshold) {
        ++below;
      }
      if (q <= 0.0) {
        ++inverted;
      }
    }
    return {worst, below, inverted};
  }

  // Node -> alive cells (CSR).
  void build_node_cells() {
    start_.assign(points_.size() + 1, 0);
    for (const Cell& cell : cells_) {
      if (cell.alive) {
        for (int k = 0; k < cell.n; ++k) {
          ++start_[u(cell.v[static_cast<std::size_t>(k)]) + 1];
        }
      }
    }
    for (std::size_t i = 0; i < points_.size(); ++i) {
      start_[i + 1] += start_[i];
    }
    node_cells_.assign(u(start_.back()), 0);
    std::vector<Index> fill(start_.begin(), start_.end() - 1);
    for (std::size_t c = 0; c < cells_.size(); ++c) {
      const Cell& cell = cells_[c];
      if (cell.alive) {
        for (int k = 0; k < cell.n; ++k) {
          node_cells_[u(fill[u(cell.v[static_cast<std::size_t>(k)])]++)] = sz(c);
        }
      }
    }
  }

  [[nodiscard]] double worst_around(std::size_t v) const {
    double worst = std::numeric_limits<double>::infinity();
    for (Index k = start_[v]; k < start_[v + 1]; ++k) {
      worst = std::min(worst, quality(cells_[u(node_cells_[u(k)])]));
    }
    return worst;
  }

  // Objective of a node position: the worst quality of its cells and the sum of squared
  // deficits below the threshold. A position is better when the deficit sum is lower and the
  // worst quality is not.
  struct Score {
    double worst = 0.0;
    double deficit = 0.0;
    [[nodiscard]] bool better_than(const Score& o) const {
      return deficit < o.deficit - 1e-15 && worst >= o.worst;
    }
  };

  [[nodiscard]] Score score(std::size_t v) const {
    Score s{std::numeric_limits<double>::infinity(), 0.0};
    for (Index k = start_[v]; k < start_[v + 1]; ++k) {
      const double q = quality(cells_[u(node_cells_[u(k)])]);
      s.worst = std::min(s.worst, q);
      if (q < opt_.threshold) {
        s.deficit += (opt_.threshold - q) * (opt_.threshold - q);
      }
    }
    return s;
  }

  Index move_nodes() {
    build_node_cells();
    // Movable nodes of cells below the threshold, worst first.
    std::vector<std::pair<double, std::size_t>> order;
    for (std::size_t v = 0; v < points_.size(); ++v) {
      if (!movable_[v] || start_[v] == start_[v + 1]) {
        continue;
      }
      const double worst = worst_around(v);
      if (worst < opt_.threshold) {
        order.emplace_back(worst, v);
      }
    }
    std::sort(order.begin(), order.end());
    Index moved = 0;
    for (const auto& [initial, v] : order) {
      Score current = score(v);
      if (current.worst >= opt_.threshold) {
        continue;
      }
      // The other nodes of the node's cells: their centre and mean distance.
      std::set<Index> others;
      for (Index k = start_[v]; k < start_[v + 1]; ++k) {
        const Cell& cell = cells_[u(node_cells_[u(k)])];
        for (int j = 0; j < cell.n; ++j) {
          const Index w = cell.v[static_cast<std::size_t>(j)];
          if (u(w) != v) {
            others.insert(w);
          }
        }
      }
      const Vec3 origin = at(points_, sz(v));
      Vec3 centre;
      double mean = 0.0;
      for (const Index w : others) {
        centre += at(points_, w);
        mean += distance(at(points_, w), origin);
      }
      centre = centre / static_cast<double>(others.size());
      mean /= static_cast<double>(others.size());
      bool improved = false;
      auto place = [&](const Vec3& p) { points_[v] = {p.x, p.y, p.z}; };
      place(centre);
      const Score at_centre = score(v);
      if (at_centre.better_than(current)) {
        current = at_centre;
        improved = true;
      } else {
        place(origin);
      }
      double step = 0.2 * mean;
      for (int iteration = 0; iteration < 60 && step > 1e-6 * mean; ++iteration) {
        const Vec3 here = at(points_, sz(v));
        Vec3 best_point = here;
        Score best = current;
        for (const auto& d : kDirections) {
          place(here + Vec3{d[0], d[1], d[2]} * step);
          const Score value = score(v);
          if (value.better_than(best)) {
            best = value;
            best_point = at(points_, sz(v));
          }
        }
        place(best_point);
        if (best_point.x != here.x || best_point.y != here.y || best_point.z != here.z) {
          current = best;
          improved = true;
        } else {
          step *= 0.5;
        }
      }
      if (improved) {
        ++moved;
      }
    }
    return moved;
  }

  // Positive orientation of a tetrahedron's node tuple (a, b, c, d) relative to cell t: true when
  // the tuple is an even permutation of the cell's nodes.
  [[nodiscard]] bool even(const Cell& t, const std::array<Index, 4>& tuple) const {
    std::array<int, 4> pos{};
    for (std::size_t i = 0; i < 4; ++i) {
      for (std::size_t k = 0; k < 4; ++k) {
        if (t.v[k] == tuple[i]) {
          pos[i] = static_cast<int>(k);
        }
      }
    }
    int inversions = 0;
    for (std::size_t i = 0; i < 4; ++i) {
      for (std::size_t j = i + 1; j < 4; ++j) {
        inversions += static_cast<int>(pos[i] > pos[j]);
      }
    }
    return inversions % 2 == 0;
  }

  [[nodiscard]] double tet_quality(const std::array<Index, 4>& t) const {
    return cell_quality_3d(points_, t.data(), 4);
  }

  std::pair<Index, Index> swap() {
    // Faces (sorted triangles) and edges (sorted pairs) of all alive cells, with their cells.
    std::vector<std::pair<std::array<Index, 3>, Index>> faces;
    std::vector<std::pair<std::array<Index, 2>, Index>> edges;
    for (std::size_t c = 0; c < cells_.size(); ++c) {
      Cell& cell = cells_[c];
      cell.touched = false;
      if (!cell.alive) {
        continue;
      }
      for (const auto& f : cell_triangles(cell.n)) {
        std::array<Index, 3> key = {cell.v[u(f[0])], cell.v[u(f[1])], cell.v[u(f[2])]};
        std::sort(key.begin(), key.end());
        faces.emplace_back(key, sz(c));
      }
      for (const auto& e : cell_edges(cell.n)) {
        std::array<Index, 2> key = {cell.v[u(e[0])], cell.v[u(e[1])]};
        std::sort(key.begin(), key.end());
        edges.emplace_back(key, sz(c));
      }
    }
    std::sort(faces.begin(), faces.end());
    std::sort(edges.begin(), edges.end());
    std::vector<double> q(cells_.size(), 1.0);
    for (std::size_t c = 0; c < cells_.size(); ++c) {
      if (cells_[c].alive && cells_[c].n == 4) {
        q[c] = quality(cells_[c]);
      }
    }
    // Candidates: (worst quality of the cells involved, kind 0 face / 1 edge, first entry).
    struct Candidate {
      double worst;
      int kind;
      std::size_t first;
      std::size_t count;
      bool operator<(const Candidate& o) const {
        if (worst != o.worst) {
          return worst < o.worst;
        }
        if (kind != o.kind) {
          return kind < o.kind;
        }
        return first < o.first;
      }
    };
    std::vector<Candidate> candidates;
    for (std::size_t i = 0; i < faces.size();) {
      std::size_t j = i;
      while (j < faces.size() && faces[j].first == faces[i].first) {
        ++j;
      }
      if (j - i == 2) {
        const Cell& a = cells_[u(faces[i].second)];
        const Cell& b = cells_[u(faces[i + 1].second)];
        const double worst = std::min(q[u(faces[i].second)], q[u(faces[i + 1].second)]);
        if (a.swappable && b.swappable && worst < opt_.threshold) {
          candidates.push_back({worst, 0, i, 2});
        }
      }
      i = j;
    }
    for (std::size_t i = 0; i < edges.size();) {
      std::size_t j = i;
      while (j < edges.size() && edges[j].first == edges[i].first) {
        ++j;
      }
      if (j - i == 3) {
        bool ok = true;
        double worst = 1.0;
        for (std::size_t k = i; k < j; ++k) {
          ok = ok && cells_[u(edges[k].second)].swappable;
          worst = std::min(worst, q[u(edges[k].second)]);
        }
        if (ok && worst < opt_.threshold) {
          candidates.push_back({worst, 1, i, 3});
        }
      }
      i = j;
    }
    std::sort(candidates.begin(), candidates.end());
    // Edges created during this pass (a 2-3 swap must not duplicate an existing edge).
    std::set<std::array<Index, 2>> new_edges;
    auto edge_exists = [&](Index a, Index b) {
      std::array<Index, 2> key = {std::min(a, b), std::max(a, b)};
      const auto it = std::lower_bound(edges.begin(), edges.end(), key,
                                       [](const std::pair<std::array<Index, 2>, Index>& e,
                                          const std::array<Index, 2>& k) { return e.first < k; });
      return (it != edges.end() && it->first == key) || new_edges.count(key) > 0;
    };
    Index s23 = 0;
    Index s32 = 0;
    auto replace = [&](const std::vector<Index>& old,
                       const std::vector<std::array<Index, 4>>& add) {
      for (const Index c : old) {
        cells_[u(c)].alive = false;
        cells_[u(c)].touched = true;
      }
      for (const auto& t : add) {
        Cell cell;
        cell.n = 4;
        std::copy(t.begin(), t.end(), cell.v.begin());
        cell.swappable = true;
        cell.touched = true;
        cells_.push_back(cell);
        for (const auto& e : cell_edges(4)) {
          new_edges.insert({std::min(t[u(e[0])], t[u(e[1])]), std::max(t[u(e[0])], t[u(e[1])])});
        }
      }
    };
    for (const Candidate& cand : candidates) {
      if (cand.kind == 0) {
        const Index ca = faces[cand.first].second;
        const Index cb = faces[cand.first + 1].second;
        if (cells_[u(ca)].touched || cells_[u(cb)].touched) {
          continue;
        }
        const auto& key = faces[cand.first].first;
        const Cell& ta = cells_[u(ca)];
        const Cell& tb = cells_[u(cb)];
        Index d = -1;
        Index e = -1;
        for (std::size_t k = 0; k < 4; ++k) {
          if (std::find(key.begin(), key.end(), ta.v[k]) == key.end()) {
            d = ta.v[k];
          }
          if (std::find(key.begin(), key.end(), tb.v[k]) == key.end()) {
            e = tb.v[k];
          }
        }
        if (d == e || edge_exists(d, e)) {
          continue;
        }
        Index a = key[0];
        Index b = key[1];
        const Index c = key[2];
        if (!even(ta, {a, b, c, d})) {
          std::swap(a, b);
        }
        const std::vector<std::array<Index, 4>> add = {{a, b, e, d}, {b, c, e, d}, {c, a, e, d}};
        double worst_new = 1.0;
        for (const auto& t : add) {
          worst_new = std::min(worst_new, tet_quality(t));
        }
        if (worst_new > cand.worst + 1e-12) {
          replace({ca, cb}, add);
          ++s23;
        }
        continue;
      }
      std::array<Index, 3> around{};
      bool skip = false;
      for (std::size_t k = 0; k < 3; ++k) {
        around[k] = edges[cand.first + k].second;
        skip = skip || cells_[u(around[k])].touched;
      }
      if (skip) {
        continue;
      }
      const auto& key = edges[cand.first].first;
      const Index d = key[0];
      const Index e = key[1];
      // The ring: each tetrahedron's two other nodes; closed ring a-b, b-c, c-a.
      std::vector<std::array<Index, 2>> pairs;
      for (const Index t : around) {
        std::array<Index, 2> pair{};
        std::size_t m = 0;
        for (std::size_t k = 0; k < 4; ++k) {
          const Index w = cells_[u(t)].v[k];
          if (w != d && w != e && m < 2) {
            pair[m++] = w;
          }
        }
        pairs.push_back(pair);
      }
      std::set<Index> ring;
      std::map<Index, int> degree;
      for (const auto& p : pairs) {
        ring.insert(p.begin(), p.end());
        ++degree[p[0]];
        ++degree[p[1]];
      }
      bool closed = ring.size() == 3;
      for (const auto& [node, count] : degree) {
        closed = closed && count == 2;
      }
      if (!closed) {
        continue;
      }
      // Orientation from the first tetrahedron: (a, b, e, d) is positive (as in a 2-3 swap).
      Index a = pairs[0][0];
      Index b = pairs[0][1];
      if (!even(cells_[u(around[0])], {a, b, e, d})) {
        std::swap(a, b);
      }
      Index c = -1;
      for (const Index w : ring) {
        if (w != a && w != b) {
          c = w;
        }
      }
      const std::vector<std::array<Index, 4>> add = {{a, b, c, d}, {b, a, c, e}};
      double worst_new = 1.0;
      for (const auto& t : add) {
        worst_new = std::min(worst_new, tet_quality(t));
      }
      if (worst_new > cand.worst + 1e-12) {
        replace({around[0], around[1], around[2]}, add);
        ++s32;
      }
    }
    // Drop dead cells, keeping the order (fixed tetrahedra first).
    std::vector<Cell> kept;
    kept.reserve(cells_.size());
    for (const Cell& cell : cells_) {
      if (cell.alive) {
        kept.push_back(cell);
      }
    }
    std::stable_partition(kept.begin(), kept.end(), [](const Cell& c) { return !c.swappable; });
    cells_ = std::move(kept);
    return {s23, s32};
  }

  std::vector<std::array<double, 3>> points_;
  std::vector<bool> movable_;
  Optimize3DOptions opt_;
  std::vector<Cell> cells_;
  std::vector<Index> start_;
  std::vector<Index> node_cells_;
};

}  // namespace

double cell_quality_3d(const std::vector<std::array<double, 3>>& points, const std::int64_t* nodes,
                       int n) {
  const CornerTable& table = corner_table(n);
  double worst = 1.0;
  for (const auto& c : table.corners) {
    const Vec3 o = at(points, nodes[c[0]]);
    const Vec3 e1 = at(points, nodes[c[1]]) - o;
    const Vec3 e2 = at(points, nodes[c[2]]) - o;
    const Vec3 e3 = at(points, nodes[c[3]]) - o;
    const double lengths = norm(e1) * norm(e2) * norm(e3);
    if (lengths <= 0.0) {
      return -1.0;
    }
    worst = std::min(worst, dot(cross(e1, e2), e3) / lengths * table.factor);
  }
  return worst;
}

Optimize3DResult optimize_3d(const Optimize3DInput& input, const Optimize3DOptions& options) {
  if (input.movable.size() != input.points.size()) {
    throw std::invalid_argument("optimize_3d: movable must have one flag per point");
  }
  if (options.passes < 0 || !(options.threshold <= 1.0)) {
    throw std::invalid_argument("optimize_3d: passes >= 0 and threshold <= 1 are required");
  }
  Optimizer optimizer(input, options);
  return optimizer.run();
}

}  // namespace og
