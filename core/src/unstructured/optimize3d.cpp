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
#include "opengrid/core/quality/quality_cells3d.hpp"
#include "opengrid/core/util/parallel.hpp"

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
// Point insertion tries at most this many cells (the worst first).
constexpr std::size_t kMaxInsertions = 1000;
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
    r.below_floor_before = below_floor();
    if (opt_.move_nodes && r.below_floor_before > 0) {
      r.moved += rescue();
    }
    if (opt_.insert && opt_.swap && below_floor() > 0) {
      r.inserted = insert_points();
      if (r.inserted > 0 && opt_.move_nodes) {
        r.moved += rescue();
      }
    }
    r.below_floor_after = below_floor();
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

  // Floor margin of a cell: >= 0 when it meets the quality floor.
  [[nodiscard]] double margin(const Cell& cell) const {
    const CellQuality3D q = cell_quality_3d_full(points_, cell.v.data(), cell.n);
    return std::min(q.scaled_jacobian - opt_.floor.min_jacobian,
                    opt_.floor.max_skewness - q.skewness);
  }

  [[nodiscard]] Index below_floor() const {
    std::vector<char> below(cells_.size(), 0);
    parallel_for(0, sz(cells_.size()), [&](Index c) {
      const Cell& cell = cells_[u(c)];
      if (cell.alive && margin(cell) < 0.0) {
        below[u(c)] = 1;
      }
    });
    return static_cast<Index>(std::count(below.begin(), below.end(), 1));
  }

  [[nodiscard]] double worst_margin(std::size_t v) const {
    double worst = std::numeric_limits<double>::infinity();
    for (Index k = start_[v]; k < start_[v + 1]; ++k) {
      worst = std::min(worst, margin(cells_[u(node_cells_[u(k)])]));
    }
    return worst;
  }

  Index rescue() {
    build_node_cells();
    const std::vector<double> worst = per_node([&](std::size_t v) { return worst_margin(v); });
    std::vector<std::pair<double, std::size_t>> order;
    for (std::size_t v = 0; v < points_.size(); ++v) {
      if (movable_[v] && start_[v] != start_[v + 1] && worst[v] < 0.0) {
        order.emplace_back(worst[v], v);
      }
    }
    std::sort(order.begin(), order.end());
    Index moved = 0;
    for (const auto& [initial, v] : order) {
      double current = worst_margin(v);
      if (current >= 0.0) {
        continue;
      }
      double mean = 0.0;
      int count = 0;
      const Vec3 origin = at(points_, sz(v));
      for (Index k = start_[v]; k < start_[v + 1]; ++k) {
        const Cell& cell = cells_[u(node_cells_[u(k)])];
        for (int j = 0; j < cell.n; ++j) {
          const Index w = cell.v[static_cast<std::size_t>(j)];
          if (u(w) != v) {
            mean += distance(at(points_, w), origin);
            ++count;
          }
        }
      }
      mean /= std::max(count, 1);
      bool improved = false;
      double step = 0.2 * mean;
      for (int iteration = 0; iteration < 80 && step > 1e-7 * mean; ++iteration) {
        const Vec3 here = at(points_, sz(v));
        Vec3 best_point = here;
        double best = current;
        for (const auto& d : kDirections) {
          const Vec3 p = here + Vec3{d[0], d[1], d[2]} * step;
          points_[v] = {p.x, p.y, p.z};
          const double value = worst_margin(v);
          if (value > best) {
            best = value;
            best_point = p;
          }
        }
        points_[v] = {best_point.x, best_point.y, best_point.z};
        if (best > current) {
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

  // Point insertion (optimize3d.hpp): the number of points inserted.
  Index insert_points() {
    build_node_cells();
    std::vector<std::pair<double, Index>> order;
    for (std::size_t c = 0; c < cells_.size(); ++c) {
      Cell& cell = cells_[c];
      cell.touched = false;
      bool fixed = true;
      for (std::size_t k = 0; k < 4; ++k) {
        fixed = fixed && !movable_[u(cell.v[k])];
      }
      if (cell.alive && cell.n == 4 && cell.swappable && fixed) {
        const double m = margin(cell);
        if (m < 0.0) {
          order.emplace_back(m, sz(c));
        }
      }
    }
    std::sort(order.begin(), order.end());
    if (order.size() > kMaxInsertions) {
      order.resize(kMaxInsertions);  // (the worst first)
    }
    Index inserted = 0;
    for (const auto& [initial, c] : order) {
      if (!cells_[u(c)].alive || cells_[u(c)].touched) {
        continue;
      }
      std::vector<Index> cavity = {c};
      for (int grow = 0; grow < 3; ++grow) {
        cavity = grown(cavity);
        if (try_star(cavity)) {
          ++inserted;
          break;
        }
      }
    }
    compact();
    return inserted;
  }

  // The alive cell other than `t` with triangle `key` (-1: none), from the node lists of the
  // pass (cells created during it are not listed: their faces count as boundary).
  [[nodiscard]] Index across(Index t, const std::array<Index, 3>& key) const {
    for (Index k = start_[u(key[0])]; k < start_[u(key[0]) + 1]; ++k) {
      const Index w = node_cells_[u(k)];
      if (w != t && cells_[u(w)].alive && has_triangle(cells_[u(w)], key)) {
        return w;
      }
    }
    return -1;
  }

  [[nodiscard]] std::array<Index, 3> face_key(const Cell& cell, const std::array<int, 3>& f) const {
    std::array<Index, 3> key = {cell.v[u(f[0])], cell.v[u(f[1])], cell.v[u(f[2])]};
    std::sort(key.begin(), key.end());
    return key;
  }

  // The cavity with the swappable tetrahedra across its faces (untouched in this pass).
  [[nodiscard]] std::vector<Index> grown(const std::vector<Index>& cavity) const {
    std::vector<Index> out = cavity;
    for (const Index t : cavity) {
      for (const auto& f : cell_triangles(4)) {
        const Index w = across(t, face_key(cells_[u(t)], f));
        if (w >= 0 && cells_[u(w)].n == 4 && cells_[u(w)].swappable && !cells_[u(w)].touched &&
            std::find(out.begin(), out.end(), w) == out.end()) {
          out.push_back(w);
        }
      }
    }
    return out;
  }

  [[nodiscard]] double volume6(const std::array<Index, 4>& t) const {
    const Vec3 a = at(points_, t[0]);
    return std::fabs(
        dot(at(points_, t[1]) - a, cross(at(points_, t[2]) - a, at(points_, t[3]) - a)));
  }

  // Replaces the cavity by the star of a new node when that is better (optimize3d.hpp).
  bool try_star(const std::vector<Index>& cavity) {
    std::vector<std::array<Index, 4>> star;
    std::vector<Index> nodes;
    double old_worst = std::numeric_limits<double>::infinity();
    double old_volume = 0.0;
    const Index p = sz(points_.size());
    for (const Index t : cavity) {
      const Cell& cell = cells_[u(t)];
      old_worst = std::min(old_worst, margin(cell));
      old_volume += volume6({cell.v[0], cell.v[1], cell.v[2], cell.v[3]});
      for (std::size_t k = 0; k < 4; ++k) {
        if (std::find(nodes.begin(), nodes.end(), cell.v[k]) == nodes.end()) {
          nodes.push_back(cell.v[k]);
        }
      }
      for (const auto& f : cell_triangles(4)) {
        const Index w = across(t, face_key(cell, f));
        if (w >= 0 && std::find(cavity.begin(), cavity.end(), w) != cavity.end()) {
          continue;  // an inner face of the cavity
        }
        // The cell with its node opposite this face replaced by the new node: same
        // orientation, positive when the node is on the same side of the face.
        std::array<Index, 4> t4 = {cell.v[0], cell.v[1], cell.v[2], cell.v[3]};
        t4[u(6 - f[0] - f[1] - f[2])] = p;
        star.push_back(t4);
      }
    }
    Vec3 centre{0.0, 0.0, 0.0};
    for (const Index v : nodes) {
      centre += at(points_, v);
    }
    centre = centre / static_cast<double>(nodes.size());
    double spread = 0.0;
    for (const Index v : nodes) {
      spread += distance(at(points_, v), centre) / static_cast<double>(nodes.size());
    }
    // Starts: the cavity's centroid and each of its cells' (a flat cavity's kernel can be
    // small); the star's worst floor margin, finite even with inverted cells (it rises as they
    // straighten), is maximised from each, the best kept.
    std::vector<Vec3> starts = {centre};
    // ... and the apex of an ideal tetrahedron on each face of the bad cell (the first), on the
    // side of its fourth node: above a layer front for a sliver lying on it.
    const Cell& bad = cells_[u(cavity.front())];
    for (const auto& f : cell_triangles(4)) {
      const Vec3 a = at(points_, bad.v[u(f[0])]);
      const Vec3 b = at(points_, bad.v[u(f[1])]);
      const Vec3 c = at(points_, bad.v[u(f[2])]);
      const Vec3 opposite = at(points_, bad.v[u(6 - f[0] - f[1] - f[2])]);
      Vec3 normal = cross(b - a, c - a);
      const double length = norm(normal);
      if (!(length > 0.0)) {
        continue;
      }
      normal = normal / length;
      const Vec3 mid = (a + b + c) / 3.0;
      if (dot(opposite - mid, normal) < 0.0) {
        normal = -normal;
      }
      const double edge = (distance(a, b) + distance(b, c) + distance(c, a)) / 3.0;
      starts.push_back(mid + normal * (std::sqrt(2.0 / 3.0) * edge));
    }
    starts.push_back((at(points_, bad.v[0]) + at(points_, bad.v[1]) + at(points_, bad.v[2]) +
                      at(points_, bad.v[3])) /
                     4.0);
    points_.push_back({centre.x, centre.y, centre.z});
    movable_.push_back(true);
    auto worst_star = [&]() {
      double worst = std::numeric_limits<double>::infinity();
      for (const auto& t4 : star) {
        Cell cell;
        std::copy(t4.begin(), t4.end(), cell.v.begin());
        worst = std::min(worst, margin(cell));
      }
      return worst;
    };
    double current = -std::numeric_limits<double>::infinity();
    Vec3 chosen = centre;
    for (const Vec3& start : starts) {
      points_[u(p)] = {start.x, start.y, start.z};
      double value = worst_star();
      double step = 0.2 * spread;
      for (int iteration = 0; iteration < 80 && step > 1e-7 * spread; ++iteration) {
        const Vec3 here = at(points_, p);
        Vec3 best_point = here;
        double best = value;
        for (const auto& d : kDirections) {
          const Vec3 q = here + Vec3{d[0], d[1], d[2]} * step;
          points_[u(p)] = {q.x, q.y, q.z};
          const double trial = worst_star();
          if (trial > best) {
            best = trial;
            best_point = q;
          }
        }
        points_[u(p)] = {best_point.x, best_point.y, best_point.z};
        if (best > value) {
          value = best;
        } else {
          step *= 0.5;
        }
      }
      if (value > current) {
        current = value;
        chosen = at(points_, p);
      }
    }
    points_[u(p)] = {chosen.x, chosen.y, chosen.z};
    double new_volume = 0.0;
    bool positive = true;
    for (const auto& t4 : star) {
      new_volume += volume6(t4);
      Cell cell;
      std::copy(t4.begin(), t4.end(), cell.v.begin());
      positive = positive && quality(cell) > 0.0;
    }
    const bool tiles =
        positive && std::fabs(new_volume - old_volume) <= 1e-9 * std::max(old_volume, 1e-300);
    // Kept only when it reaches the floor: otherwise the cell stays for the caller, which can
    // still remedy cells whose nodes are all fixed (the volume generator lowers the layers).
    if (!tiles || !(current >= 0.0) || !(current > old_worst)) {
      points_.pop_back();
      movable_.pop_back();
      return false;
    }
    for (const Index t : cavity) {
      cells_[u(t)].alive = false;
      cells_[u(t)].touched = true;
    }
    for (const auto& t4 : star) {
      Cell cell;
      std::copy(t4.begin(), t4.end(), cell.v.begin());
      cell.swappable = true;
      cell.touched = true;
      cells_.push_back(cell);
    }
    return true;
  }

  // Drops dead cells, keeping the order with the fixed ones first (as after a swap pass).
  void compact() {
    std::vector<Cell> kept;
    kept.reserve(cells_.size());
    for (const Cell& cell : cells_) {
      if (cell.alive) {
        kept.push_back(cell);
      }
    }
    std::stable_partition(kept.begin(), kept.end(), [](const Cell& c) { return !c.swappable; });
    cells_ = std::move(kept);
  }

  // f(v) for every movable node with cells (others: +infinity), in parallel.
  template <class F>
  [[nodiscard]] std::vector<double> per_node(F&& f) const {
    std::vector<double> out(points_.size(), std::numeric_limits<double>::infinity());
    parallel_for(0, sz(points_.size()), [&](Index i) {
      const std::size_t v = u(i);
      if (movable_[v] && start_[v] != start_[v + 1]) {
        out[v] = f(v);
      }
    });
    return out;
  }

  [[nodiscard]] std::tuple<double, Index, Index> summary() const {
    std::vector<double> values(cells_.size(), 1.0);
    parallel_for(0, sz(cells_.size()), [&](Index c) {
      if (cells_[u(c)].alive) {
        values[u(c)] = quality(cells_[u(c)]);
      }
    });
    double worst = 1.0;
    Index below = 0;
    Index inverted = 0;
    for (std::size_t c = 0; c < cells_.size(); ++c) {
      if (!cells_[c].alive) {
        continue;
      }
      const double q = values[c];
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
    const std::vector<double> worst = per_node([&](std::size_t v) { return worst_around(v); });
    std::vector<std::pair<double, std::size_t>> order;
    for (std::size_t v = 0; v < points_.size(); ++v) {
      if (movable_[v] && start_[v] != start_[v + 1] && worst[v] < opt_.threshold) {
        order.emplace_back(worst[v], v);
      }
    }
    std::sort(order.begin(), order.end());
    Index moved = 0;
    for (const auto& entry : order) {
      // Not a structured binding: lambdas below capture v (clang rejects capturing structured
      // bindings when OpenMP is enabled).
      const auto v = entry.second;
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

  // Whether cell w has the sorted triangle `key` among its triangular faces.
  [[nodiscard]] static bool has_triangle(const Cell& w, const std::array<Index, 3>& key) {
    for (const auto& f : cell_triangles(w.n)) {
      std::array<Index, 3> t = {w.v[u(f[0])], w.v[u(f[1])], w.v[u(f[2])]};
      std::sort(t.begin(), t.end());
      if (t == key) {
        return true;
      }
    }
    return false;
  }

  // Whether cell w has the edge (a, b), a < b.
  [[nodiscard]] static bool has_edge(const Cell& w, Index a, Index b) {
    for (const auto& e : cell_edges(w.n)) {
      const Index x = w.v[u(e[0])];
      const Index y = w.v[u(e[1])];
      if (std::min(x, y) == a && std::max(x, y) == b) {
        return true;
      }
    }
    return false;
  }

  std::pair<Index, Index> swap() {
    // Node -> cells at the start of the pass (cells killed during the pass stay listed, new ones
    // are tracked separately), ascending cell index.
    build_node_cells();
    for (Cell& cell : cells_) {
      cell.touched = false;
    }
    std::vector<double> q(cells_.size(), 1.0);
    parallel_for(0, sz(cells_.size()), [&](Index c) {
      if (cells_[u(c)].alive && cells_[u(c)].n == 4) {
        q[u(c)] = quality(cells_[u(c)]);
      }
    });
    // Candidates: a face of exactly two cells (2-3 swap) or an edge of exactly three cells (3-2
    // swap), all swappable tetrahedra, the worst below the threshold - so every candidate has a
    // tetrahedron below the threshold, and only those are searched. Order: worst quality, kind
    // (0 face, 1 edge), sorted node key; cells by ascending index.
    struct Candidate {
      double worst;
      int kind;
      std::array<Index, 3> key;
      std::array<Index, 3> cells;
      bool operator<(const Candidate& o) const {
        if (worst != o.worst) {
          return worst < o.worst;
        }
        if (kind != o.kind) {
          return kind < o.kind;
        }
        return key < o.key;
      }
    };
    std::vector<Candidate> candidates;
    const std::size_t n_start = cells_.size();
    for (std::size_t c = 0; c < n_start; ++c) {
      const Cell& cell = cells_[c];
      if (!cell.alive || cell.n != 4 || !cell.swappable || !(q[c] < opt_.threshold)) {
        continue;
      }
      for (const auto& f : cell_triangles(4)) {
        std::array<Index, 3> key = {cell.v[u(f[0])], cell.v[u(f[1])], cell.v[u(f[2])]};
        std::sort(key.begin(), key.end());
        std::vector<Index> around;
        for (Index k = start_[u(key[0])]; k < start_[u(key[0]) + 1]; ++k) {
          const Index w = node_cells_[u(k)];
          if (has_triangle(cells_[u(w)], key)) {
            around.push_back(w);
          }
        }
        if (around.size() != 2 || !cells_[u(around[0])].swappable ||
            !cells_[u(around[1])].swappable) {
          continue;
        }
        const double worst = std::min(q[u(around[0])], q[u(around[1])]);
        candidates.push_back({worst, 0, key, {around[0], around[1], -1}});
      }
      for (const auto& e : cell_edges(4)) {
        const Index a = std::min(cell.v[u(e[0])], cell.v[u(e[1])]);
        const Index b = std::max(cell.v[u(e[0])], cell.v[u(e[1])]);
        std::vector<Index> around;
        for (Index k = start_[u(a)]; k < start_[u(a) + 1]; ++k) {
          const Index w = node_cells_[u(k)];
          if (has_edge(cells_[u(w)], a, b)) {
            around.push_back(w);
          }
        }
        if (around.size() != 3) {
          continue;
        }
        bool ok = true;
        double worst = 1.0;
        for (const Index w : around) {
          ok = ok && cells_[u(w)].swappable;
          worst = std::min(worst, q[u(w)]);
        }
        if (ok) {
          candidates.push_back({worst, 1, {a, b, -1}, {around[0], around[1], around[2]}});
        }
      }
    }
    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end(),
                                 [](const Candidate& x, const Candidate& y) {
                                   return x.kind == y.kind && x.key == y.key;
                                 }),
                     candidates.end());
    // Edges created during this pass (a 2-3 swap must not duplicate an existing edge).
    std::set<std::array<Index, 2>> new_edges;
    auto edge_exists = [&](Index a, Index b) {
      const std::array<Index, 2> key = {std::min(a, b), std::max(a, b)};
      for (Index k = start_[u(key[0])]; k < start_[u(key[0]) + 1]; ++k) {
        if (has_edge(cells_[u(node_cells_[u(k)])], key[0], key[1])) {
          return true;
        }
      }
      return new_edges.count(key) > 0;
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
        const Index ca = cand.cells[0];
        const Index cb = cand.cells[1];
        if (cells_[u(ca)].touched || cells_[u(cb)].touched) {
          continue;
        }
        const auto& key = cand.key;
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
        around[k] = cand.cells[k];
        skip = skip || cells_[u(around[k])].touched;
      }
      if (skip) {
        continue;
      }
      const auto& key = cand.key;
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
    compact();
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
