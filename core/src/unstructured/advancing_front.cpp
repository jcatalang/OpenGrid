#include "opengrid/core/unstructured/advancing_front.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <queue>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "opengrid/core/util/predicates.hpp"

namespace og {
namespace {

using Index = MeshIndex;

std::size_t u(Index i) {
  return static_cast<std::size_t>(i);
}

std::uint64_t directed_key(Index a, Index b) {
  return (static_cast<std::uint64_t>(a) << 32U) ^ static_cast<std::uint64_t>(b);
}

std::uint64_t undirected_key(Index a, Index b) {
  if (b < a) {
    std::swap(a, b);
  }
  return directed_key(a, b);
}

// Uniform bucket grid over a fixed bounding box; ids may appear in several cells.
class Grid {
 public:
  Grid(double xmin, double ymin, double xmax, double ymax, double cell) : x0_(xmin), y0_(ymin) {
    cell_ = std::max(cell, 1e-300);
    nx_ = std::max(1L, static_cast<long>((xmax - xmin) / cell_) + 1);
    ny_ = std::max(1L, static_cast<long>((ymax - ymin) / cell_) + 1);
    while (nx_ * ny_ > 4'000'000) {
      cell_ *= 2.0;
      nx_ = std::max(1L, static_cast<long>((xmax - xmin) / cell_) + 1);
      ny_ = std::max(1L, static_cast<long>((ymax - ymin) / cell_) + 1);
    }
    cells_.assign(static_cast<std::size_t>(nx_ * ny_), {});
  }

  void insert(Index id, double xa, double ya, double xb, double yb) {
    for_cells(xa, ya, xb, yb, [&](std::vector<Index>& bucket) { bucket.push_back(id); });
  }

  void erase(Index id, double xa, double ya, double xb, double yb) {
    for_cells(xa, ya, xb, yb, [&](std::vector<Index>& bucket) {
      const auto it = std::find(bucket.begin(), bucket.end(), id);
      if (it != bucket.end()) {
        *it = bucket.back();
        bucket.pop_back();
      }
    });
  }

  template <class F>
  void query(double xa, double ya, double xb, double yb, F&& f) {
    for_cells(xa, ya, xb, yb, [&](std::vector<Index>& bucket) {
      for (const Index id : bucket) {
        f(id);
      }
    });
  }

 private:
  [[nodiscard]] long clamp_x(double x) const {
    return std::clamp(static_cast<long>(std::floor((x - x0_) / cell_)), 0L, nx_ - 1);
  }
  [[nodiscard]] long clamp_y(double y) const {
    return std::clamp(static_cast<long>(std::floor((y - y0_) / cell_)), 0L, ny_ - 1);
  }

  template <class F>
  void for_cells(double xa, double ya, double xb, double yb, F&& f) {
    const long i0 = clamp_x(std::min(xa, xb));
    const long i1 = clamp_x(std::max(xa, xb));
    const long j0 = clamp_y(std::min(ya, yb));
    const long j1 = clamp_y(std::max(ya, yb));
    for (long j = j0; j <= j1; ++j) {
      for (long i = i0; i <= i1; ++i) {
        f(cells_[static_cast<std::size_t>(j * nx_ + i)]);
      }
    }
  }

  double x0_;
  double y0_;
  double cell_ = 1.0;
  long nx_ = 1;
  long ny_ = 1;
  std::vector<std::vector<Index>> cells_;
};

struct FrontEdge {
  Index a;
  Index b;
  bool alive = true;
  int tries = 0;
};

constexpr int kMaxTries = 3;
// Size change per layer: the target size of a front edge of length l is clamped to
// [kShrink, kGrowth] * l, so the front follows the size field gradually.
constexpr double kShrink = 0.55;
constexpr double kGrowth = 1.3;
// Existing front nodes closer than kReuseRadius * h to the ideal apex are tried before a new
// point. Coarsening (e.g. around a small hole) needs the neighbouring apexes, which lie about
// one edge length l >= h / kGrowth from the ideal apex, to qualify: kReuseRadius must exceed
// 1 / kGrowth (0.77) with a margin, or the front stops coarsening and keeps the boundary spacing
// (measured: median edge / target size 0.48 at 0.8 vs 1.01 at 0.9 on a graded cylinder).
constexpr double kReuseRadius = 0.9;

class Front {
 public:
  Front(const TriangulationInput& input, const RefinementOptions& options)
      : p_(input.points), segments_(input.segments) {
    double xmin = p_[0][0];
    double xmax = xmin;
    double ymin = p_[0][1];
    double ymax = ymin;
    for (const auto& q : p_) {
      xmin = std::min(xmin, q[0]);
      xmax = std::max(xmax, q[0]);
      ymin = std::min(ymin, q[1]);
      ymax = std::max(ymax, q[1]);
    }
    std::vector<double> sum(p_.size(), 0.0);
    std::vector<int> count(p_.size(), 0);
    double total = 0.0;
    for (const auto& s : segments_) {
      const double len = dist(s[0], s[1]);
      total += len;
      for (const Index v : s) {
        sum[u(v)] += len;
        count[u(v)] += 1;
      }
    }
    std::vector<std::array<double, 2>> bp;
    std::vector<double> bs;
    for (std::size_t i = 0; i < p_.size(); ++i) {
      if (count[i] > 0) {
        bp.push_back(p_[i]);
        bs.push_back(sum[i] / count[i]);
      }
    }
    size_ = std::make_unique<SizeField>(bp, bs, options.size);
    const double mean = total / static_cast<double>(segments_.size());
    nodes_ = std::make_unique<Grid>(xmin, ymin, xmax, ymax, 2.0 * mean);
    edges_grid_ = std::make_unique<Grid>(xmin, ymin, xmax, ymax, 2.0 * mean);
    degree_.assign(p_.size(), 0);
    for (Index i = 0; std::cmp_less(i, p_.size()); ++i) {
      nodes_->insert(i, p_[u(i)][0], p_[u(i)][1], p_[u(i)][0], p_[u(i)][1]);
    }
    for (const auto& s : segments_) {
      add_edge(s[0], s[1]);
    }
    max_points_ = options.max_points;
  }

  void run(TriangleMesh& out) {
    while (!heap_.empty()) {
      const auto [length, tries, order, id] = heap_.top();
      heap_.pop();
      FrontEdge& e = edges_[u(id)];
      if (!e.alive || e.tries != tries) {
        continue;
      }
      if (std::cmp_greater_equal(p_.size(), max_points_)) {
        out.warnings.emplace_back("advancing front stopped at the maximum number of points");
        break;
      }
      if (advance(id)) {
        continue;
      }
      e.tries += 1;
      if (e.tries <= kMaxTries) {
        heap_.emplace(length, e.tries, sequence_++, id);
      }
    }
    close_remaining(out);
  }

  void output(TriangleMesh& out) {
    out.points = p_;
    out.triangles = triangles_;
  }

  std::vector<std::array<double, 2>>& points() {
    return p_;
  }
  std::vector<std::array<Index, 3>>& triangles() {
    return triangles_;
  }

 private:
  [[nodiscard]] double dist(Index a, Index b) const {
    return std::hypot(p_[u(a)][0] - p_[u(b)][0], p_[u(a)][1] - p_[u(b)][1]);
  }

  void add_edge(Index a, Index b) {
    const auto id = static_cast<Index>(edges_.size());
    edges_.push_back({a, b});
    map_[directed_key(a, b)] = id;
    degree_[u(a)] += 1;
    degree_[u(b)] += 1;
    edges_grid_->insert(id, p_[u(a)][0], p_[u(a)][1], p_[u(b)][0], p_[u(b)][1]);
    heap_.emplace(dist(a, b), 0, sequence_++, id);
  }

  void remove_edge(Index id) {
    FrontEdge& e = edges_[u(id)];
    e.alive = false;
    map_.erase(directed_key(e.a, e.b));
    edges_grid_->erase(id, p_[u(e.a)][0], p_[u(e.a)][1], p_[u(e.b)][0], p_[u(e.b)][1]);
    for (const Index v : {e.a, e.b}) {
      degree_[u(v)] -= 1;
      if (degree_[u(v)] == 0) {
        // The node is now interior to the mesh: it can never be on the front again.
        nodes_->erase(v, p_[u(v)][0], p_[u(v)][1], p_[u(v)][0], p_[u(v)][1]);
      }
    }
  }

  // Closed segments pq and rs intersect (sharing an endpoint counts only if they overlap).
  bool crosses(const double* pp, const double* q, const double* r, const double* s,
               bool share) const {
    if (std::max(pp[0], q[0]) < std::min(r[0], s[0]) ||
        std::max(r[0], s[0]) < std::min(pp[0], q[0]) ||
        std::max(pp[1], q[1]) < std::min(r[1], s[1]) ||
        std::max(r[1], s[1]) < std::min(pp[1], q[1])) {
      return false;  // disjoint bounding boxes
    }
    const int o1 = orient2d(pp, q, r);
    const int o2 = orient2d(pp, q, s);
    if (o1 * o2 > 0) {
      return false;  // r and s strictly on the same side of pq
    }
    const int o3 = orient2d(r, s, pp);
    const int o4 = orient2d(r, s, q);
    if (share) {
      // Shared endpoint: only a collinear overlap is a conflict.
      if (o1 == 0 && o2 == 0) {
        return overlap(pp, q, r, s);
      }
      return false;
    }
    if (o1 * o2 < 0 && o3 * o4 < 0) {
      return true;
    }
    if (o1 == 0 && on_segment(pp, q, r)) {
      return true;
    }
    if (o2 == 0 && on_segment(pp, q, s)) {
      return true;
    }
    if (o3 == 0 && on_segment(r, s, pp)) {
      return true;
    }
    if (o4 == 0 && on_segment(r, s, q)) {
      return true;
    }
    return false;
  }

  static bool on_segment(const double* a, const double* b, const double* c) {
    return std::min(a[0], b[0]) <= c[0] && c[0] <= std::max(a[0], b[0]) &&
           std::min(a[1], b[1]) <= c[1] && c[1] <= std::max(a[1], b[1]);
  }

  static bool overlap(const double* a, const double* b, const double* c, const double* d) {
    // Collinear segments sharing an endpoint overlap if the other ends point the same way.
    const double* shared = nullptr;
    const double* x = nullptr;
    const double* y = nullptr;
    if (a[0] == c[0] && a[1] == c[1]) {
      shared = a;
      x = b;
      y = d;
    } else if (a[0] == d[0] && a[1] == d[1]) {
      shared = a;
      x = b;
      y = c;
    } else if (b[0] == c[0] && b[1] == c[1]) {
      shared = b;
      x = a;
      y = d;
    } else {
      shared = b;
      x = a;
      y = c;
    }
    return (x[0] - shared[0]) * (y[0] - shared[0]) + (x[1] - shared[1]) * (y[1] - shared[1]) > 0.0;
  }

  static double segment_distance(const std::array<double, 2>& q, const std::array<double, 2>& a,
                                 const std::array<double, 2>& b) {
    const double dx = b[0] - a[0];
    const double dy = b[1] - a[1];
    const double len2 = dx * dx + dy * dy;
    double t = 0.0;
    if (len2 > 0.0) {
      t = std::clamp(((q[0] - a[0]) * dx + (q[1] - a[1]) * dy) / len2, 0.0, 1.0);
    }
    return std::hypot(q[0] - (a[0] + t * dx), q[1] - (a[1] + t * dy));
  }

  // Can triangle (a, b, c) be added? c < 0 denotes the new point `q`.
  bool valid(Index a, Index b, Index c, const std::array<double, 2>& q, double h) {
    const double* pa = p_[u(a)].data();
    const double* pb = p_[u(b)].data();
    const double* pc = q.data();
    if (orient2d(pa, pb, pc) <= 0) {
      return false;
    }
    const double xa = std::min({pa[0], pb[0], pc[0]});
    const double xb = std::max({pa[0], pb[0], pc[0]});
    const double ya = std::min({pa[1], pb[1], pc[1]});
    const double yb = std::max({pa[1], pb[1], pc[1]});
    const double margin = 0.5 * h;
    bool ok = true;
    stamp_ += 1;
    edges_grid_->query(xa - margin, ya - margin, xb + margin, yb + margin, [&](Index id) {
      if (!ok || seen_[u(id)] == stamp_) {
        return;
      }
      seen_[u(id)] = stamp_;
      const FrontEdge& e = edges_[u(id)];
      if (!e.alive || (e.a == a && e.b == b)) {
        return;
      }
      // Front edges the triangle closes against: (b, c) and (c, a), traversed as in the triangle.
      if (c >= 0 && ((e.a == b && e.b == c) || (e.a == c && e.b == a))) {
        return;
      }
      const double* pu = p_[u(e.a)].data();
      const double* pv = p_[u(e.b)].data();
      for (const Index end : {a, b}) {
        const bool share = e.a == end || e.b == end || (c >= 0 && (e.a == c || e.b == c));
        if (crosses(p_[u(end)].data(), pc, pu, pv, share)) {
          ok = false;
          return;
        }
      }
      if (c < 0 && segment_distance(q, p_[u(e.a)], p_[u(e.b)]) < 0.3 * h) {
        ok = false;
      }
    });
    if (!ok) {
      return false;
    }
    nodes_->query(xa - margin, ya - margin, xb + margin, yb + margin, [&](Index v) {
      if (!ok || degree_[u(v)] == 0 || v == a || v == b || v == c) {
        return;
      }
      const double* pv = p_[u(v)].data();
      if (orient2d(pa, pb, pv) >= 0 && orient2d(pb, pc, pv) >= 0 && orient2d(pc, pa, pv) >= 0) {
        ok = false;  // inside or on the triangle
        return;
      }
      if (c < 0 && std::hypot(pv[0] - q[0], pv[1] - q[1]) < 0.5 * h) {
        ok = false;  // a new point too close to the front
      }
    });
    return ok;
  }

  bool advance(Index id) {
    const Index a = edges_[u(id)].a;
    const Index b = edges_[u(id)].b;
    const auto& pa = p_[u(a)];
    const auto& pb = p_[u(b)];
    const double l = dist(a, b);
    const std::array<double, 2> m = {0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1])};
    const double nx = -(pb[1] - pa[1]) / l;
    const double ny = (pb[0] - pa[0]) / l;
    const double h = std::clamp((*size_)(m[0], m[1]), kShrink * l, kGrowth * l);
    const double height = std::sqrt(h * h - 0.25 * l * l);
    const std::array<double, 2> ideal = {m[0] + nx * height, m[1] + ny * height};
    // The search widens on every retry of this edge.
    const double radius = (1.3 + 0.7 * edges_[u(id)].tries) * h;

    std::vector<std::pair<double, Index>>& nearby = nearby_;
    nearby.clear();
    if (seen_.size() < edges_.size()) {
      seen_.resize(edges_.size() + 1024, 0);
    }
    if (node_seen_.size() < p_.size() + 1) {
      node_seen_.resize(p_.size() + 1024, 0);
    }
    node_stamp_ += 1;
    nodes_->query(ideal[0] - radius, ideal[1] - radius, ideal[0] + radius, ideal[1] + radius,
                  [&](Index v) {
                    if (degree_[u(v)] == 0 || v == a || v == b || node_seen_[u(v)] == node_stamp_) {
                      return;
                    }
                    node_seen_[u(v)] = node_stamp_;
                    const double d = std::hypot(p_[u(v)][0] - ideal[0], p_[u(v)][1] - ideal[1]);
                    if (d <= radius && orient2d(pa.data(), pb.data(), p_[u(v)].data()) > 0) {
                      nearby.emplace_back(d, v);
                    }
                  });
    std::sort(nearby.begin(), nearby.end());
    // Candidate order: close existing nodes, the ideal new point, then farther nodes.
    std::vector<Index>& order = order_;
    order.clear();
    bool new_added = false;
    for (const auto& [d, v] : nearby) {
      if (!new_added && d >= kReuseRadius * h) {
        order.push_back(-1);
        new_added = true;
      }
      order.push_back(v);
    }
    if (!new_added) {
      order.push_back(-1);
    }
    // The first valid candidate wins.
    Index best = -2;
    std::array<double, 2> best_point{};
    for (const Index c : order) {
      std::array<double, 2> q = ideal;
      if (c >= 0) {
        q = p_[u(c)];
      }
      if (valid(a, b, c, q, h)) {
        best = c;
        best_point = q;
        break;
      }
    }
    if (best == -2) {
      return false;
    }
    Index apex = best;
    if (best < 0) {
      apex = static_cast<Index>(p_.size());
      p_.push_back(best_point);
      degree_.push_back(0);
      nodes_->insert(apex, best_point[0], best_point[1], best_point[0], best_point[1]);
    }
    triangles_.push_back({a, b, apex});
    remove_edge(id);
    // Sides (b, apex) and (apex, a): close against the front or extend it.
    for (const auto& [x, y] : {std::pair{b, apex}, std::pair{apex, a}}) {
      const auto it = map_.find(directed_key(x, y));
      if (it != map_.end()) {
        remove_edge(it->second);
      } else {
        add_edge(y, x);
      }
    }
    return true;
  }

  void close_remaining(TriangleMesh& out) {
    std::vector<std::array<Index, 2>> open;
    for (const FrontEdge& e : edges_) {
      if (e.alive) {
        open.push_back({e.a, e.b});
      }
    }
    if (open.empty()) {
      return;
    }
    TriangulationInput rest;
    std::unordered_map<Index, Index> local;
    std::vector<Index> global;
    for (const auto& s : open) {
      std::array<Index, 2> seg{};
      for (std::size_t k = 0; k < 2; ++k) {
        const auto [it, inserted] = local.try_emplace(s[k], static_cast<Index>(global.size()));
        if (inserted) {
          global.push_back(s[k]);
          rest.points.push_back(p_[u(s[k])]);
        }
        seg[k] = it->second;
      }
      rest.segments.push_back(seg);
    }
    RefinementOptions options;
    options.size = size_options_;
    const TriangleMesh patch = triangulate_domain(rest, options);
    for (std::size_t i = global.size(); i < patch.points.size(); ++i) {
      global.push_back(static_cast<Index>(p_.size()));
      p_.push_back(patch.points[i]);
    }
    for (const auto& t : patch.triangles) {
      triangles_.push_back({global[u(t[0])], global[u(t[1])], global[u(t[2])]});
    }
    out.warnings.push_back("advancing front left " + std::to_string(open.size()) +
                           " edges open; closed them with the constrained Delaunay mesher");
  }

 public:
  SizeFieldOptions size_options_;

 private:
  std::vector<std::array<double, 2>> p_;
  std::vector<std::array<Index, 2>> segments_;
  std::vector<FrontEdge> edges_;
  std::unordered_map<std::uint64_t, Index> map_;
  std::vector<int> degree_;
  std::vector<std::array<Index, 3>> triangles_;
  std::unique_ptr<SizeField> size_;
  std::unique_ptr<Grid> nodes_;
  std::unique_ptr<Grid> edges_grid_;
  std::vector<std::uint64_t> seen_;
  std::uint64_t stamp_ = 0;
  std::vector<std::uint64_t> node_seen_;
  std::uint64_t node_stamp_ = 0;
  std::vector<std::pair<double, Index>> nearby_;
  std::vector<Index> order_;
  std::uint64_t sequence_ = 0;
  Index max_points_ = 0;
  std::priority_queue<std::tuple<double, int, std::uint64_t, Index>,
                      std::vector<std::tuple<double, int, std::uint64_t, Index>>, std::greater<>>
      heap_;
};

}  // namespace

std::int64_t flip_to_delaunay(const std::vector<std::array<double, 2>>& points,
                              std::vector<std::array<MeshIndex, 3>>& triangles,
                              const std::vector<std::array<MeshIndex, 2>>& segments) {
  const std::size_t nt = triangles.size();
  // Neighbour across the edge opposite vertex k (-1 on the domain boundary), from sorted
  // half-edges.
  std::vector<std::array<Index, 3>> nb(nt, {-1, -1, -1});
  {
    std::vector<std::pair<std::uint64_t, Index>> half;
    half.reserve(3 * nt);
    for (std::size_t t = 0; t < nt; ++t) {
      for (std::size_t k = 0; k < 3; ++k) {
        half.emplace_back(undirected_key(triangles[t][(k + 1) % 3], triangles[t][(k + 2) % 3]),
                          static_cast<Index>(3 * t + k));
      }
    }
    std::sort(half.begin(), half.end());
    for (std::size_t i = 0; i + 1 < half.size(); ++i) {
      if (half[i].first == half[i + 1].first) {
        const Index x = half[i].second;
        const Index y = half[i + 1].second;
        nb[u(x / 3)][u(x % 3)] = y / 3;
        nb[u(y / 3)][u(y % 3)] = x / 3;
        ++i;
      }
    }
  }
  // Segments are never flipped (in a single domain they are all on its boundary and have no
  // neighbour; the set also protects interior constraints if the caller has any).
  std::unordered_set<std::uint64_t> fixed;
  for (const auto& s : segments) {
    fixed.insert(undirected_key(s[0], s[1]));
  }
  auto replace = [&](Index t, Index old_nb, Index new_nb) {
    if (t < 0) {
      return;
    }
    for (Index& v : nb[u(t)]) {
      if (v == old_nb) {
        v = new_nb;
        return;
      }
    }
  };
  std::vector<Index> stack;
  stack.reserve(nt);
  for (Index t = static_cast<Index>(nt) - 1; t >= 0; --t) {
    stack.push_back(t);
  }
  std::int64_t flips = 0;
  const std::int64_t limit = 100 * static_cast<std::int64_t>(nt) + 1000;
  while (!stack.empty() && flips < limit) {
    const Index t1 = stack.back();
    stack.pop_back();
    for (std::size_t k1 = 0; k1 < 3; ++k1) {
      const Index t2 = nb[u(t1)][k1];
      if (t2 < 0) {
        continue;
      }
      const auto tri1 = triangles[u(t1)];
      const Index p0 = tri1[k1];
      const Index p1 = tri1[(k1 + 1) % 3];
      const Index p2 = tri1[(k1 + 2) % 3];
      if (fixed.count(undirected_key(p1, p2)) != 0) {
        continue;
      }
      const auto tri2 = triangles[u(t2)];
      std::size_t k2 = 0;
      while (tri2[k2] == p1 || tri2[k2] == p2) {
        ++k2;
      }
      const Index q = tri2[k2];
      const double* a = points[u(p0)].data();
      if (incircle(a, points[u(p1)].data(), points[u(p2)].data(), points[u(q)].data()) <= 0) {
        continue;
      }
      if (orient2d(a, points[u(q)].data(), points[u(p1)].data()) >= 0 ||
          orient2d(a, points[u(q)].data(), points[u(p2)].data()) <= 0) {
        continue;  // not a strictly convex quadrilateral
      }
      // t2 traverses the shared edge backwards: t2 = (q, p2, p1) starting at k2.
      const Index n1a = nb[u(t1)][(k1 + 1) % 3];         // across (p2, p0)
      const Index n1b = nb[u(t1)][(k1 + 2) % 3];         // across (p0, p1)
      const Index across_p1q = nb[u(t2)][(k2 + 1) % 3];  // opposite p2
      const Index across_qp2 = nb[u(t2)][(k2 + 2) % 3];  // opposite p1
      triangles[u(t1)] = {p0, p1, q};
      triangles[u(t2)] = {p0, q, p2};
      nb[u(t1)] = {across_p1q, t2, n1b};
      nb[u(t2)] = {across_qp2, n1a, t1};
      replace(across_p1q, t2, t1);
      replace(n1a, t1, t2);
      flips += 1;
      stack.push_back(t2);
      stack.push_back(t1);
      break;
    }
  }
  return flips;
}

std::int64_t smooth_interior(std::vector<std::array<double, 2>>& points,
                             const std::vector<std::array<MeshIndex, 3>>& triangles,
                             std::size_t n_fixed, int passes) {
  // Vertex -> incident triangles.
  std::vector<std::vector<Index>> incident(points.size());
  for (Index t = 0; std::cmp_less(t, triangles.size()); ++t) {
    for (const Index v : triangles[u(t)]) {
      incident[u(v)].push_back(t);
    }
  }
  // Sine of the smallest angle around v (the smallest angle of a triangle is at most 60
  // degrees, where the sine is monotone); -1 if a triangle is inverted.
  auto worst_angle = [&](Index v) {
    double worst = 2.0;
    for (const Index t : incident[u(v)]) {
      const auto& tri = triangles[u(t)];
      const auto& a = points[u(tri[0])];
      const auto& b = points[u(tri[1])];
      const auto& c = points[u(tri[2])];
      if (orient2d(a.data(), b.data(), c.data()) <= 0) {
        return -1.0;
      }
      const double la = std::hypot(b[0] - c[0], b[1] - c[1]);
      const double lb = std::hypot(c[0] - a[0], c[1] - a[1]);
      const double lc = std::hypot(a[0] - b[0], a[1] - b[1]);
      const double area2 = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
      // The smallest angle is opposite the shortest side, between the two longer ones.
      const double longest = std::max({la, lb, lc});
      const double middle = la + lb + lc - longest - std::min({la, lb, lc});
      worst = std::min(worst, area2 / (longest * middle));
    }
    return worst;
  };
  std::int64_t moved = 0;
  for (int pass = 0; pass < passes; ++pass) {
    for (std::size_t v = n_fixed; v < points.size(); ++v) {
      if (incident[v].empty()) {
        continue;
      }
      // Centroid of the neighbouring vertices.
      double sx = 0.0;
      double sy = 0.0;
      double count = 0.0;
      for (const Index t : incident[v]) {
        for (const Index w : triangles[u(t)]) {
          if (u(w) != v) {
            sx += points[u(w)][0];
            sy += points[u(w)][1];
            count += 1.0;
          }
        }
      }
      const std::array<double, 2> old = points[v];
      const double before = worst_angle(static_cast<Index>(v));
      points[v] = {sx / count, sy / count};
      if (worst_angle(static_cast<Index>(v)) < before) {
        points[v] = old;  // never make the neighbourhood worse
      } else {
        moved += 1;
      }
    }
  }
  return moved;
}

TriangleMesh advancing_front(const TriangulationInput& input, const RefinementOptions& options) {
  if (input.points.size() < 3 || input.segments.size() < 3) {
    throw std::invalid_argument("advancing front: at least 3 boundary points and segments");
  }
  for (const auto& s : input.segments) {
    for (const MeshIndex v : s) {
      if (v < 0 || std::cmp_greater_equal(v, input.points.size())) {
        throw std::invalid_argument("advancing front: segment refers to a missing point");
      }
    }
  }
  if (input.points.size() >= (std::size_t{1} << 31U)) {
    throw std::invalid_argument("advancing front: too many points");
  }
  // Validate the boundary (duplicates, closure, orientation, intersections) exactly as the
  // Delaunay mesher does, without refining.
  RefinementOptions check = options;
  check.refine = false;
  (void)triangulate_domain(input, check);

  Front front(input, options);
  front.size_options_ = options.size;
  TriangleMesh out;
  front.run(out);
  // Post-processing: Delaunay flips, guarded Laplacian smoothing of interior nodes, flips again.
  flip_to_delaunay(front.points(), front.triangles(), input.segments);
  smooth_interior(front.points(), front.triangles(), input.points.size(), 3);
  flip_to_delaunay(front.points(), front.triangles(), input.segments);
  front.output(out);
  return out;
}

}  // namespace og
