#include "opengrid/core/unstructured/triangulate.hpp"

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

#include "opengrid/core/util/predicates.hpp"

namespace og {
namespace {

using Index = MeshIndex;
constexpr Index kNone = -1;

std::size_t u(Index i) {
  return static_cast<std::size_t>(i);
}

int next(int k) {
  return (k + 1) % 3;
}

int prev(int k) {
  return (k + 2) % 3;
}

// Triangle with vertices v (counter-clockwise); n[k] and c[k] refer to the edge opposite v[k],
// i.e. (v[k+1], v[k+2]).
struct Tri {
  std::array<Index, 3> v{};
  std::array<Index, 3> n{kNone, kNone, kNone};
  std::array<bool, 3> c{false, false, false};
  bool alive = true;
  bool inside = false;
};

struct CavityEdge {
  Index a;
  Index b;
  Index outer;
  bool constrained;
};

class Mesher {
 public:
  explicit Mesher(const TriangulationInput& input)
      : n_input_(static_cast<Index>(input.points.size())) {
    if (input.points.size() < 3) {
      throw std::invalid_argument("triangulate: at least 3 boundary points are required");
    }
    p_ = input.points;
    check_duplicates();
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
    scale_ = std::max(xmax - xmin, ymax - ymin);
    const double cx = 0.5 * (xmin + xmax);
    const double cy = 0.5 * (ymin + ymax);
    const double m = std::max(scale_, 1e-300);
    // Super triangle far outside (exact predicates make its size harmless).
    p_.push_back({cx - 20.0 * m, cy - 10.0 * m});
    p_.push_back({cx + 20.0 * m, cy - 10.0 * m});
    p_.push_back({cx, cy + 20.0 * m});
    vt_.assign(p_.size(), kNone);
    Tri super;
    super.v = {n_input_, n_input_ + 1, n_input_ + 2};
    t_.push_back(super);
    for (Index k = 0; k < 3; ++k) {
      vt_[u(n_input_ + k)] = 0;
    }
  }

  void insert_boundary_points() {
    Index hint = 0;
    for (Index i = 0; i < n_input_; ++i) {
      bool crossed = false;
      const Index start = locate(p_[u(i)].data(), hint, crossed);
      std::vector<Index> tris;
      std::vector<CavityEdge> edges;
      if (start == kNone || !cavity(p_[u(i)].data(), start, tris, edges)) {
        throw std::invalid_argument("triangulate: degenerate boundary point " + std::to_string(i));
      }
      hint = commit(i, tris, edges, false);
    }
  }

  void recover_segments(const std::vector<std::array<Index, 2>>& segments) {
    for (const auto& s : segments) {
      recover(s[0], s[1]);
    }
    restore_delaunay();
  }

  void classify(const std::vector<std::array<Index, 2>>& segments) {
    std::vector<Index> stack;
    for (const auto& s : segments) {
      const auto [ti, k] = find_edge(s[0], s[1]);
      if (ti == kNone) {
        throw std::logic_error("triangulate: segment lost after recovery");
      }
      Index left = ti;
      if (t_[u(ti)].v[u(next(k))] != s[0]) {
        left = t_[u(ti)].n[u(k)];
      }
      if (left != kNone && !t_[u(left)].inside) {
        t_[u(left)].inside = true;
        stack.push_back(left);
      }
    }
    while (!stack.empty()) {
      const Index ti = stack.back();
      stack.pop_back();
      const Tri& tri = t_[u(ti)];
      for (const Index vertex : tri.v) {
        if (vertex >= n_input_ && vertex < n_input_ + 3) {
          throw std::invalid_argument(
              "triangulate: the boundary is not closed around the domain (check the loop "
              "orientation: outer counter-clockwise, holes clockwise)");
        }
      }
      for (int k = 0; k < 3; ++k) {
        const Index nb = tri.n[u(k)];
        if (tri.c[u(k)] || nb == kNone || t_[u(nb)].inside) {
          continue;
        }
        t_[u(nb)].inside = true;
        stack.push_back(nb);
      }
    }
  }

  void refine(const RefinementOptions& options, const std::vector<std::array<Index, 2>>& segments,
              TriangleMesh& out) {
    // Boundary sizes: mean length of the segments at each node.
    std::vector<double> sum(u(n_input_), 0.0);
    std::vector<int> count(u(n_input_), 0);
    for (const auto& s : segments) {
      const double len = dist(s[0], s[1]);
      for (const Index v : s) {
        sum[u(v)] += len;
        count[u(v)] += 1;
      }
    }
    std::vector<std::array<double, 2>> bp;
    std::vector<double> bs;
    for (Index i = 0; i < n_input_; ++i) {
      if (count[u(i)] > 0) {
        bp.push_back(p_[u(i)]);
        bs.push_back(sum[u(i)] / count[u(i)]);
      }
    }
    const SizeField size(bp, bs, options.size);
    const double ratio_limit =
        1.0 / (2.0 * std::sin(options.min_angle_deg * std::numbers::pi / 180.0));

    std::deque<Index> queue;
    for (Index ti = 0; std::cmp_less(ti, t_.size()); ++ti) {
      if (t_[u(ti)].alive && t_[u(ti)].inside) {
        queue.push_back(ti);
      }
    }
    bool capped = false;
    while (!queue.empty()) {
      const Index ti = queue.front();
      queue.pop_front();
      if (!t_[u(ti)].alive || !t_[u(ti)].inside) {
        continue;
      }
      const Tri& tri = t_[u(ti)];
      const auto& a = p_[u(tri.v[0])];
      const auto& b = p_[u(tri.v[1])];
      const auto& c = p_[u(tri.v[2])];
      std::array<double, 2> center{};
      double radius = 0.0;
      if (!circumcenter(a, b, c, center, radius)) {
        continue;
      }
      const double shortest = std::min({distance2(a, b), distance2(b, c), distance2(c, a)});
      const std::array<double, 2> centroid = {(a[0] + b[0] + c[0]) / 3.0,
                                              (a[1] + b[1] + c[1]) / 3.0};
      // A triangle standing on a boundary segment cannot become smaller than that segment
      // (boundary segments are never split), so the segment length bounds its target size.
      double boundary_edge = 0.0;
      bool on_boundary = false;
      for (int k = 0; k < 3; ++k) {
        if (tri.c[u(k)]) {
          on_boundary = true;
          boundary_edge = std::max(boundary_edge, dist(tri.v[u(next(k))], tri.v[u(prev(k))]));
        }
      }
      const double h = std::max(size(centroid[0], centroid[1]), boundary_edge);
      const double equivalent_edge = radius * std::numbers::sqrt3;
      const bool too_big = equivalent_edge > 1.2 * h;
      // Shape refinement stops well below the local size (small input angles cannot be fixed
      // without splitting the boundary and would otherwise refine forever).
      const bool bad_shape =
          radius > ratio_limit * std::sqrt(shortest) * (1.0 + 1e-12) && equivalent_edge > 0.02 * h;
      if (!too_big && !bad_shape) {
        continue;
      }
      if (radius <= 1e-12 * scale_) {
        continue;
      }
      if (static_cast<Index>(p_.size()) - 3 >= options.max_points) {
        capped = true;
        break;
      }
      // Shape refinement inserts the off-center (Ungor): on the bisector of the shortest edge,
      // where the new triangle on that edge has exactly the minimum angle, or the circumcenter
      // if that is closer. It is closer to the triangle than the circumcenter of a thin
      // triangle, so it is far more often inside the domain and away from the boundary.
      std::array<double, 2> target = center;
      if (bad_shape && !too_big) {
        target = off_center(tri, center, options.min_angle_deg);
      }
      Index created = try_insert(target, ti, queue);
      // Fallback for oversized interior triangles only: a centroid split does not improve
      // angles (it would cascade into ever thinner triangles, so thin ones are left as they
      // are: their size is held by a boundary segment that cannot be split), and on a boundary
      // triangle it leaves a triangle on the same segment.
      const bool thin = radius > ratio_limit * std::sqrt(shortest) * (1.0 + 1e-12);
      if (created == kNone && too_big && !on_boundary && !thin) {
        created = try_insert(centroid, ti, queue);
      }
      if (created == kNone) {
        out.rejected += 1;
      }
    }
    if (capped) {
      out.warnings.push_back("refinement stopped at the maximum number of points (" +
                             std::to_string(options.max_points) + ")");
    }
  }

  void output(TriangleMesh& out) const {
    std::vector<Index> map(p_.size(), kNone);
    for (Index i = 0; i < n_input_; ++i) {
      map[u(i)] = i;
      out.points.push_back(p_[u(i)]);
    }
    for (std::size_t i = u(n_input_ + 3); i < p_.size(); ++i) {
      map[i] = static_cast<Index>(out.points.size());
      out.points.push_back(p_[i]);
    }
    for (const Tri& tri : t_) {
      if (tri.alive && tri.inside) {
        out.triangles.push_back({map[u(tri.v[0])], map[u(tri.v[1])], map[u(tri.v[2])]});
      }
    }
  }

 private:
  [[nodiscard]] const double* point_ptr(Index i) const {
    return p_[u(i)].data();
  }

  [[nodiscard]] double dist(Index a, Index b) const {
    return std::hypot(p_[u(a)][0] - p_[u(b)][0], p_[u(a)][1] - p_[u(b)][1]);
  }

  static double distance2(const std::array<double, 2>& a, const std::array<double, 2>& b) {
    const double dx = a[0] - b[0];
    const double dy = a[1] - b[1];
    return dx * dx + dy * dy;
  }

  [[nodiscard]] std::array<double, 2> off_center(const Tri& tri,
                                                 const std::array<double, 2>& center,
                                                 double min_angle_deg) const {
    int shortest = 0;
    double length = std::numeric_limits<double>::infinity();
    for (int k = 0; k < 3; ++k) {
      const double l = dist(tri.v[u(next(k))], tri.v[u(prev(k))]);
      if (l < length) {
        length = l;
        shortest = k;
      }
    }
    const auto& p = p_[u(tri.v[u(next(shortest))])];
    const auto& q = p_[u(tri.v[u(prev(shortest))])];
    const std::array<double, 2> mid = {0.5 * (p[0] + q[0]), 0.5 * (p[1] + q[1])};
    const double dx = center[0] - mid[0];
    const double dy = center[1] - mid[1];
    const double d = std::hypot(dx, dy);
    const double wanted = 0.5 * length / std::tan(0.5 * min_angle_deg * std::numbers::pi / 180.0);
    if (d <= wanted || d == 0.0) {
      return center;
    }
    return {mid[0] + dx * wanted / d, mid[1] + dy * wanted / d};
  }

  static bool circumcenter(const std::array<double, 2>& a, const std::array<double, 2>& b,
                           const std::array<double, 2>& c, std::array<double, 2>& center,
                           double& radius) {
    const double bx = b[0] - a[0];
    const double by = b[1] - a[1];
    const double cx = c[0] - a[0];
    const double cy = c[1] - a[1];
    const double d = 2.0 * (bx * cy - by * cx);
    if (d == 0.0) {
      return false;
    }
    const double b2 = bx * bx + by * by;
    const double c2 = cx * cx + cy * cy;
    const double ux = (cy * b2 - by * c2) / d;
    const double uy = (bx * c2 - cx * b2) / d;
    center = {a[0] + ux, a[1] + uy};
    radius = std::hypot(ux, uy);
    return true;
  }

  void check_duplicates() const {
    std::vector<std::size_t> order(p_.size());
    for (std::size_t i = 0; i < order.size(); ++i) {
      order[i] = i;
    }
    std::sort(order.begin(), order.end(), [this](std::size_t a, std::size_t b) {
      if (p_[a][0] != p_[b][0]) {
        return p_[a][0] < p_[b][0];
      }
      if (p_[a][1] != p_[b][1]) {
        return p_[a][1] < p_[b][1];
      }
      return a < b;
    });
    for (std::size_t k = 1; k < order.size(); ++k) {
      if (p_[order[k]] == p_[order[k - 1]]) {
        throw std::invalid_argument("triangulate: duplicate points " +
                                    std::to_string(order[k - 1]) + " and " +
                                    std::to_string(order[k]));
      }
    }
  }

  static int edge_index(const Tri& tri, Index a, Index b) {
    for (int k = 0; k < 3; ++k) {
      const Index x = tri.v[u(next(k))];
      const Index y = tri.v[u(prev(k))];
      if ((x == a && y == b) || (x == b && y == a)) {
        return k;
      }
    }
    return -1;
  }

  static int vertex_index(const Tri& tri, Index a) {
    for (int k = 0; k < 3; ++k) {
      if (tri.v[u(k)] == a) {
        return k;
      }
    }
    return -1;
  }

  // Point location by a visibility walk from `start`; `crossed` reports whether the walk
  // stepped over a constrained edge (the point is then not visible from `start`).
  Index locate(const double* q, Index start, bool& crossed) const {
    crossed = false;
    Index current = start;
    if (current == kNone || !t_[u(current)].alive) {
      current = any_alive();
    }
    const std::size_t limit = 3 * t_.size() + 100;
    int rotation = 0;
    for (std::size_t step = 0; step < limit; ++step) {
      const Tri& tri = t_[u(current)];
      bool moved = false;
      for (int m = 0; m < 3; ++m) {
        const int k = (m + rotation) % 3;
        if (orient2d(point_ptr(tri.v[u(next(k))]), point_ptr(tri.v[u(prev(k))]), q) < 0) {
          if (tri.c[u(k)]) {
            crossed = true;
          }
          if (tri.n[u(k)] == kNone) {
            return kNone;
          }
          current = tri.n[u(k)];
          moved = true;
          break;
        }
      }
      if (!moved) {
        return current;
      }
      rotation = (rotation + 1) % 3;
    }
    // Pathological walk: exhaustive search (visibility unknown, so treat as crossed).
    crossed = true;
    for (Index ti = 0; std::cmp_less(ti, t_.size()); ++ti) {
      const Tri& tri = t_[u(ti)];
      if (tri.alive && orient2d(point_ptr(tri.v[0]), point_ptr(tri.v[1]), q) >= 0 &&
          orient2d(point_ptr(tri.v[1]), point_ptr(tri.v[2]), q) >= 0 &&
          orient2d(point_ptr(tri.v[2]), point_ptr(tri.v[0]), q) >= 0) {
        return ti;
      }
    }
    return kNone;
  }

  [[nodiscard]] Index any_alive() const {
    for (Index ti = static_cast<Index>(t_.size()) - 1; ti >= 0; --ti) {
      if (t_[u(ti)].alive) {
        return ti;
      }
    }
    return kNone;
  }

  // Bowyer-Watson cavity of q (never crossing constrained edges). Returns false if q is not
  // strictly inside the cavity polygon (degenerate insertion).
  bool cavity(const double* q, Index start, std::vector<Index>& tris,
              std::vector<CavityEdge>& edges) {
    generation_ += 1;
    if (stamp_.size() < t_.size()) {
      stamp_.resize(t_.size(), 0);
    }
    tris.clear();
    edges.clear();
    tris.push_back(start);
    stamp_[u(start)] = generation_;
    for (std::size_t head = 0; head < tris.size(); ++head) {
      const Tri& tri = t_[u(tris[head])];
      for (int k = 0; k < 3; ++k) {
        const Index nb = tri.n[u(k)];
        if (tri.c[u(k)] || nb == kNone || stamp_[u(nb)] == generation_) {
          continue;
        }
        const Tri& other = t_[u(nb)];
        if (incircle(point_ptr(other.v[0]), point_ptr(other.v[1]), point_ptr(other.v[2]), q) > 0) {
          stamp_[u(nb)] = generation_;
          tris.push_back(nb);
        }
      }
    }
    for (const Index ti : tris) {
      const Tri& tri = t_[u(ti)];
      for (int k = 0; k < 3; ++k) {
        const Index nb = tri.n[u(k)];
        if (nb != kNone && stamp_[u(nb)] == generation_ && !tri.c[u(k)]) {
          continue;
        }
        const Index a = tri.v[u(next(k))];
        const Index b = tri.v[u(prev(k))];
        if (orient2d(point_ptr(a), point_ptr(b), q) <= 0) {
          return false;
        }
        edges.push_back({a, b, nb, tri.c[u(k)]});
      }
    }
    return true;
  }

  Index commit(Index qi, const std::vector<Index>& tris, const std::vector<CavityEdge>& edges,
               bool inside) {
    for (const Index ti : tris) {
      t_[u(ti)].alive = false;
    }
    std::unordered_map<Index, Index> by_start;
    std::unordered_map<Index, Index> by_end;
    std::vector<Index> created;
    for (const CavityEdge& e : edges) {
      Tri tri;
      tri.v = {e.a, e.b, qi};
      tri.n = {kNone, kNone, e.outer};
      tri.c = {false, false, e.constrained};
      tri.inside = inside;
      const auto id = static_cast<Index>(t_.size());
      t_.push_back(tri);
      if (e.outer != kNone) {
        const int k = edge_index(t_[u(e.outer)], e.a, e.b);
        t_[u(e.outer)].n[u(k)] = id;
      }
      by_start[e.a] = id;
      by_end[e.b] = id;
      created.push_back(id);
    }
    for (const Index id : created) {
      Tri& tri = t_[u(id)];
      tri.n[0] = by_start.at(tri.v[1]);
      tri.n[1] = by_end.at(tri.v[0]);
      vt_[u(tri.v[0])] = id;
      vt_[u(tri.v[1])] = id;
    }
    if (vt_.size() <= u(qi)) {
      vt_.resize(u(qi) + 1, kNone);
    }
    vt_[u(qi)] = created.front();
    return created.front();
  }

  // Attempt to insert q as a refinement point; returns a new triangle or kNone if rejected.
  Index try_insert(const std::array<double, 2>& q, Index from, std::deque<Index>& queue) {
    bool crossed = false;
    const Index start = locate(q.data(), from, crossed);
    if (start == kNone || crossed || !t_[u(start)].inside) {
      return kNone;
    }
    std::vector<Index> tris;
    std::vector<CavityEdge> edges;
    if (!cavity(q.data(), start, tris, edges)) {
      return kNone;
    }
    for (const CavityEdge& e : edges) {
      if (!e.constrained) {
        continue;
      }
      // Inside the diametral lens of a boundary segment (the segment subtends more than 120
      // degrees at q; Shewchuk's lens, thinner than the diametral circle): the triangle q would
      // form on the segment is too obtuse to be repaired, so q is rejected.
      const auto& a = p_[u(e.a)];
      const auto& b = p_[u(e.b)];
      const double ax = a[0] - q[0];
      const double ay = a[1] - q[1];
      const double bx = b[0] - q[0];
      const double by = b[1] - q[1];
      if (ax * bx + ay * by < -0.5 * std::hypot(ax, ay) * std::hypot(bx, by)) {
        return kNone;
      }
    }
    const auto qi = static_cast<Index>(p_.size());
    p_.push_back(q);
    const auto first = static_cast<Index>(t_.size());
    const Index created = commit(qi, tris, edges, true);
    for (Index id = first; std::cmp_less(id, t_.size()); ++id) {
      queue.push_back(id);
    }
    return created;
  }

  [[nodiscard]] std::vector<Index> triangles_around(Index a) const {
    std::vector<Index> result;
    Index start = vt_[u(a)];
    if (start == kNone || !t_[u(start)].alive || vertex_index(t_[u(start)], a) < 0) {
      start = kNone;
      for (Index ti = 0; std::cmp_less(ti, t_.size()); ++ti) {
        if (t_[u(ti)].alive && vertex_index(t_[u(ti)], a) >= 0) {
          start = ti;
          break;
        }
      }
    }
    if (start == kNone) {
      return result;
    }
    result.push_back(start);
    for (std::size_t head = 0; head < result.size(); ++head) {
      const Tri& tri = t_[u(result[head])];
      for (const Index nb : tri.n) {
        if (nb == kNone || vertex_index(t_[u(nb)], a) < 0) {
          continue;
        }
        if (std::find(result.begin(), result.end(), nb) == result.end()) {
          result.push_back(nb);
        }
      }
    }
    return result;
  }

  [[nodiscard]] std::pair<Index, int> find_edge(Index a, Index b) const {
    for (const Index ti : triangles_around(a)) {
      const int k = edge_index(t_[u(ti)], a, b);
      if (k >= 0) {
        return {ti, k};
      }
    }
    return {kNone, -1};
  }

  void set_constrained(Index a, Index b) {
    const auto [ti, k] = find_edge(a, b);
    t_[u(ti)].c[u(k)] = true;
    const Index nb = t_[u(ti)].n[u(k)];
    if (nb != kNone) {
      t_[u(nb)].c[u(edge_index(t_[u(nb)], a, b))] = true;
    }
  }

  // Replace the diagonal shared by t1 (edge opposite t1.v[k1]) and its neighbour.
  void flip(Index t1, int k1) {
    const Index t2 = t_[u(t1)].n[u(k1)];
    const Tri a = t_[u(t1)];
    const Tri b = t_[u(t2)];
    const Index p0 = a.v[u(k1)];
    const Index p1 = a.v[u(next(k1))];
    const Index p2 = a.v[u(prev(k1))];
    const int k2 = edge_index(b, p1, p2);
    const Index q = b.v[u(k2)];
    // b is (q, p2, p1) counter-clockwise.
    const int bq = k2;
    const int bp2 = next(bq);
    const int bp1 = prev(bq);
    const Index n_a = a.n[u(prev(k1))];  // edge (p0, p1)
    const bool c_a = a.c[u(prev(k1))];
    const Index n_b = a.n[u(next(k1))];  // edge (p2, p0)
    const bool c_b = a.c[u(next(k1))];
    const Index n_c = b.n[u(bp2)];  // edge (p1, q)
    const bool c_c = b.c[u(bp2)];
    const Index n_d = b.n[u(bp1)];  // edge (q, p2)
    const bool c_d = b.c[u(bp1)];
    Tri& x = t_[u(t1)];
    x.v = {p0, p1, q};
    x.n = {n_c, t2, n_a};
    x.c = {c_c, false, c_a};
    Tri& y = t_[u(t2)];
    y.v = {p0, q, p2};
    y.n = {n_d, n_b, t1};
    y.c = {c_d, c_b, false};
    if (n_c != kNone) {
      t_[u(n_c)].n[u(edge_index(t_[u(n_c)], p1, q))] = t1;
    }
    if (n_b != kNone) {
      t_[u(n_b)].n[u(edge_index(t_[u(n_b)], p2, p0))] = t2;
    }
    vt_[u(p0)] = t1;
    vt_[u(p1)] = t1;
    vt_[u(q)] = t1;
    vt_[u(p2)] = t2;
  }

  void recover(Index a0, Index b0) {
    std::vector<std::pair<Index, Index>> stack = {{a0, b0}};
    while (!stack.empty()) {
      const auto [a, b] = stack.back();
      stack.pop_back();
      if (find_edge(a, b).first != kNone) {
        set_constrained(a, b);
        continue;
      }
      // March from a to b collecting the crossed edges (u right, w left of a -> b).
      std::vector<std::pair<Index, Index>> crossing;
      Index split_at = kNone;
      Index current = kNone;
      Index ur = kNone;
      Index wl = kNone;
      for (const Index ti : triangles_around(a)) {
        const Tri& tri = t_[u(ti)];
        const int ia = vertex_index(tri, a);
        const Index uu = tri.v[u(next(ia))];
        const Index ww = tri.v[u(prev(ia))];
        for (const Index x : {uu, ww}) {
          if (orient2d(point_ptr(a), point_ptr(b), point_ptr(x)) == 0 && between(a, b, x)) {
            split_at = x;
          }
        }
        if (split_at != kNone) {
          break;
        }
        if (orient2d(point_ptr(a), point_ptr(b), point_ptr(uu)) < 0 &&
            orient2d(point_ptr(a), point_ptr(b), point_ptr(ww)) > 0) {
          current = ti;
          ur = uu;
          wl = ww;
          break;
        }
      }
      while (split_at == kNone && current != kNone) {
        const int k = edge_index(t_[u(current)], ur, wl);
        if (t_[u(current)].c[u(k)]) {
          throw std::invalid_argument("triangulate: boundary segments intersect (segment " +
                                      std::to_string(a) + "-" + std::to_string(b) + ")");
        }
        crossing.emplace_back(ur, wl);
        const Index nb = t_[u(current)].n[u(k)];
        const Index x = t_[u(nb)].v[u(edge_index(t_[u(nb)], ur, wl))];
        if (x == b) {
          break;
        }
        const int ox = orient2d(point_ptr(a), point_ptr(b), point_ptr(x));
        if (ox == 0) {
          split_at = x;
          break;
        }
        if (ox < 0) {
          ur = x;
        } else {
          wl = x;
        }
        current = nb;
      }
      if (split_at != kNone) {
        stack.emplace_back(split_at, b);
        stack.emplace_back(a, split_at);
        continue;
      }
      if (current == kNone) {
        throw std::logic_error("triangulate: cannot start segment recovery");
      }
      std::deque<std::pair<Index, Index>> queue(crossing.begin(), crossing.end());
      std::size_t guard = 0;
      const std::size_t limit = 100 * (queue.size() + 1) * (queue.size() + 1) + 1000;
      while (!queue.empty()) {
        if (++guard > limit) {
          throw std::logic_error("triangulate: segment recovery does not terminate");
        }
        const auto [x, y] = queue.front();
        queue.pop_front();
        const auto [t1, k1] = find_edge(x, y);
        const Index t2 = t_[u(t1)].n[u(k1)];
        const Index p0 = t_[u(t1)].v[u(k1)];
        const Index q = t_[u(t2)].v[u(edge_index(t_[u(t2)], x, y))];
        const int sx = orient2d(point_ptr(p0), point_ptr(q), point_ptr(x));
        const int sy = orient2d(point_ptr(p0), point_ptr(q), point_ptr(y));
        if (sx * sy < 0) {
          flip(t1, k1);
          if (p0 != a && p0 != b && q != a && q != b &&
              orient2d(point_ptr(a), point_ptr(b), point_ptr(p0)) *
                      orient2d(point_ptr(a), point_ptr(b), point_ptr(q)) <
                  0) {
            queue.emplace_back(p0, q);
          }
        } else {
          queue.emplace_back(x, y);
        }
      }
      set_constrained(a, b);
    }
  }

  [[nodiscard]] bool between(Index a, Index b, Index x) const {
    const double dx = p_[u(b)][0] - p_[u(a)][0];
    const double dy = p_[u(b)][1] - p_[u(a)][1];
    const double t = (p_[u(x)][0] - p_[u(a)][0]) * dx + (p_[u(x)][1] - p_[u(a)][1]) * dy;
    return t > 0.0 && t < dx * dx + dy * dy;
  }

  void restore_delaunay() {
    bool changed = true;
    while (changed) {
      changed = false;
      for (Index ti = 0; std::cmp_less(ti, t_.size()); ++ti) {
        if (!t_[u(ti)].alive) {
          continue;
        }
        for (int k = 0; k < 3; ++k) {
          const Tri& tri = t_[u(ti)];
          const Index nb = tri.n[u(k)];
          if (tri.c[u(k)] || nb == kNone) {
            continue;
          }
          const Index q =
              t_[u(nb)].v[u(edge_index(t_[u(nb)], tri.v[u(next(k))], tri.v[u(prev(k))]))];
          if (incircle(point_ptr(tri.v[0]), point_ptr(tri.v[1]), point_ptr(tri.v[2]),
                       point_ptr(q)) > 0) {
            flip(ti, k);
            changed = true;
            break;
          }
        }
      }
    }
  }

  Index n_input_;
  std::vector<std::array<double, 2>> p_;
  std::vector<Tri> t_;
  std::vector<Index> vt_;
  std::vector<std::uint64_t> stamp_;
  std::uint64_t generation_ = 0;
  double scale_ = 1.0;
};

}  // namespace

TriangleMesh triangulate_domain(const TriangulationInput& input, const RefinementOptions& options) {
  for (const auto& s : input.segments) {
    for (const MeshIndex v : s) {
      if (v < 0 || std::cmp_greater_equal(v, input.points.size())) {
        throw std::invalid_argument("triangulate: segment refers to point " + std::to_string(v) +
                                    " out of range");
      }
    }
    if (s[0] == s[1]) {
      throw std::invalid_argument("triangulate: degenerate segment");
    }
  }
  if (!(options.min_angle_deg > 0.0 && options.min_angle_deg < 34.0)) {
    throw std::invalid_argument("triangulate: min angle must lie in (0, 34) degrees");
  }
  Mesher mesher(input);
  mesher.insert_boundary_points();
  mesher.recover_segments(input.segments);
  mesher.classify(input.segments);
  TriangleMesh out;
  if (options.refine) {
    mesher.refine(options, input.segments, out);
  }
  mesher.output(out);
  return out;
}

}  // namespace og
