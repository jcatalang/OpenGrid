#include "opengrid/core/unstructured/inflation3d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
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

Vec3 normalized(const Vec3& v) {
  const double n = norm(v);
  if (n == 0.0) {
    return v;
  }
  return v / n;
}

// Newell vector (twice the area vector) of a closed polygon.
Vec3 newell(const std::vector<Vec3>& p) {
  Vec3 n;
  for (std::size_t i = 0; i < p.size(); ++i) {
    n += cross(p[i], p[(i + 1) % p.size()]);
  }
  return n;
}

Vec3 min3(const Vec3& a, const Vec3& b) {
  return {std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)};
}

Vec3 max3(const Vec3& a, const Vec3& b) {
  return {std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)};
}

double axis(const Vec3& v, int k) {
  if (k == 0) {
    return v.x;
  }
  if (k == 1) {
    return v.y;
  }
  return v.z;
}

// A triangle for intersection queries: coordinates, vertex tags (point ids, or -(column + 1) for a
// column's pending node), the wall face it belongs to (or -1) and its face group (static faces,
// else -1). Oriented out of the region (fronts: towards their wall).
struct Tri {
  std::array<Vec3, 3> p;
  std::array<Index, 3> tag{};
  Index owner = -1;
  Index group = -1;
};

bool share_tag(const Tri& t, const std::array<Index, 3>& tags, std::size_t count) {
  for (std::size_t i = 0; i < count; ++i) {
    for (const Index a : t.tag) {
      if (a == tags[i]) {
        return true;
      }
    }
  }
  return false;
}

// Bounding-volume hierarchy over triangles (median split on the longest axis).
class TriTree {
 public:
  explicit TriTree(std::vector<Tri> tris) : tris_(std::move(tris)) {
    order_.resize(tris_.size());
    std::iota(order_.begin(), order_.end(), std::size_t{0});
    if (!tris_.empty()) {
      build(0, tris_.size());
    }
  }

  // Calls f(tri) for every triangle whose box overlaps [lo, hi].
  template <class F>
  void query(const Vec3& lo, const Vec3& hi, F&& f) const {
    if (nodes_.empty()) {
      return;
    }
    std::vector<std::size_t> stack{0};
    while (!stack.empty()) {
      const Node& node = nodes_[stack.back()];
      stack.pop_back();
      if (node.hi.x < lo.x || node.hi.y < lo.y || node.hi.z < lo.z || hi.x < node.lo.x ||
          hi.y < node.lo.y || hi.z < node.lo.z) {
        continue;
      }
      if (node.left == 0) {
        for (std::size_t i = node.begin; i < node.end; ++i) {
          f(tris_[order_[i]]);
        }
        continue;
      }
      stack.push_back(node.left);
      stack.push_back(node.right);
    }
  }

 private:
  struct Node {
    Vec3 lo;
    Vec3 hi;
    std::size_t left = 0;  // 0: leaf
    std::size_t right = 0;
    std::size_t begin = 0;
    std::size_t end = 0;
  };

  std::size_t build(std::size_t begin, std::size_t end) {
    const std::size_t id = nodes_.size();
    nodes_.emplace_back();
    Vec3 lo = tris_[order_[begin]].p[0];
    Vec3 hi = lo;
    for (std::size_t i = begin; i < end; ++i) {
      for (const Vec3& q : tris_[order_[i]].p) {
        lo = min3(lo, q);
        hi = max3(hi, q);
      }
    }
    nodes_[id].lo = lo;
    nodes_[id].hi = hi;
    nodes_[id].begin = begin;
    nodes_[id].end = end;
    if (end - begin <= 8) {
      return id;
    }
    int k = 0;
    const Vec3 extent = hi - lo;
    if (extent.y > extent.x && extent.y >= extent.z) {
      k = 1;
    } else if (extent.z > extent.x && extent.z > extent.y) {
      k = 2;
    }
    const std::size_t mid = begin + (end - begin) / 2;
    auto centre = [&](std::size_t t) {
      const Tri& tri = tris_[t];
      return axis(tri.p[0], k) + axis(tri.p[1], k) + axis(tri.p[2], k);
    };
    std::nth_element(order_.begin() + static_cast<std::ptrdiff_t>(begin),
                     order_.begin() + static_cast<std::ptrdiff_t>(mid),
                     order_.begin() + static_cast<std::ptrdiff_t>(end),
                     [&](std::size_t a, std::size_t b) {
                       const double ca = centre(a);
                       const double cb = centre(b);
                       if (ca != cb) {
                         return ca < cb;
                       }
                       return a < b;
                     });
    const std::size_t left = build(begin, mid);
    const std::size_t right = build(mid, end);
    nodes_[id].left = left;
    nodes_[id].right = right;
    return id;
  }

  std::vector<Tri> tris_;
  std::vector<std::size_t> order_;
  std::vector<Node> nodes_;
};

constexpr double kTol = 1e-9;

// Segment [a, b] meets the closed triangle (Moller-Trumbore with a relative tolerance, so
// touching counts); parallel segments never do. Returns the parameter along the segment.
bool segment_hits(const Vec3& a, const Vec3& b, const Tri& t, double* param) {
  const Vec3 d = b - a;
  const Vec3 e1 = t.p[1] - t.p[0];
  const Vec3 e2 = t.p[2] - t.p[0];
  const Vec3 pv = cross(d, e2);
  const double det = dot(e1, pv);
  if (std::abs(det) <= 1e-12 * norm(d) * norm(e1) * norm(e2)) {
    return false;
  }
  const double inv = 1.0 / det;
  const Vec3 s = a - t.p[0];
  const double bu = dot(s, pv) * inv;
  if (bu < -kTol || bu > 1.0 + kTol) {
    return false;
  }
  const Vec3 q = cross(s, e1);
  const double bv = dot(d, q) * inv;
  if (bv < -kTol || bu + bv > 1.0 + kTol) {
    return false;
  }
  const double tt = dot(e2, q) * inv;
  if (tt < -kTol || tt > 1.0 + kTol) {
    return false;
  }
  *param = tt;
  return true;
}

bool triangles_cross(const Tri& a, const Tri& b) {
  double t = 0.0;
  for (int i = 0; i < 3; ++i) {
    if (segment_hits(a.p[u(i)], a.p[u((i + 1) % 3)], b, &t) ||
        segment_hits(b.p[u(i)], b.p[u((i + 1) % 3)], a, &t)) {
      return true;
    }
  }
  return false;
}

enum class Kind { tet, pyramid, prism, hex };

struct Element {
  Kind kind = Kind::tet;
  std::array<Index, 8> v{};
  bool full = false;  // an undegenerate prism / hexahedron (scaled Jacobian threshold applies)
};

constexpr std::array<std::array<int, 4>, 6> kPrismCorners = {{
    {0, 1, 2, 3},
    {1, 2, 0, 4},
    {2, 0, 1, 5},
    {3, 5, 4, 0},
    {4, 3, 5, 1},
    {5, 4, 3, 2},
}};

constexpr std::array<std::array<int, 4>, 8> kHexCorners = {{
    {0, 1, 3, 4},
    {1, 2, 0, 5},
    {2, 3, 1, 6},
    {3, 0, 2, 7},
    {4, 7, 5, 0},
    {5, 4, 6, 1},
    {6, 5, 7, 2},
    {7, 6, 4, 3},
}};

double corner(const Vec3& o, const Vec3& a, const Vec3& b, const Vec3& c) {
  const Vec3 e1 = a - o;
  const Vec3 e2 = b - o;
  const Vec3 e3 = c - o;
  const double scale = norm(e1) * norm(e2) * norm(e3);
  if (scale == 0.0) {
    return 0.0;
  }
  return dot(e1, cross(e2, e3)) / scale;
}

struct Column {
  Index base = 0;
  Index cur = 0;
  bool active = true;
  int kind = 0;  // 0 free, 1 sliding, 2 fixed
  Index group = -1;
  Vec3 plane;
  double height = 0.0;
  int layers = 0;
  double origin_spacing = 0.0;
  Vec3 normal;
  double scale = 1.0;
  Vec3 next;
  double target = 0.0;
  std::string reason;
  std::vector<Index> faces;  // incident wall faces
  std::vector<Index> nbrs;   // neighbouring columns (front edges)
};

struct WallFace {
  Index face = 0;
  std::vector<Index> cols;
  bool split = false;         // a quad continuing as two triangles
  int diag = 0;               // 0: diagonal (0, 2); 1: (1, 3)
  std::vector<Index> across;  // per edge (i, i + 1): the face across it
  Vec3 inward;
};

class Inflator {
 public:
  Inflator(const Inflation3DInput& input, const Inflation3DOptions& options)
      : in_(input), opt_(options) {
    for (const auto& p : input.points) {
      points_.emplace_back(p);
    }
  }

  Inflation3DResult run() {
    setup();
    int max_layers = opt_.max_layers;
    if (opt_.total_thickness) {
      // Smallest N with first_height * (g^N - 1) / (g - 1) >= total_thickness (as in 2D).
      double sum = 0.0;
      int n = 0;
      double h = opt_.first_height;
      while (sum < *opt_.total_thickness && n < 100000) {
        sum += h;
        h *= opt_.growth_rate;
        ++n;
      }
      max_layers = std::min(max_layers, n);
    }
    double h = opt_.first_height;
    for (int k = 0; k < max_layers; ++k) {
      if (!any_active()) {
        break;
      }
      march(h);
      validate(h);
      commit(h);
      h *= opt_.growth_rate;
    }
    return finish();
  }

 private:
  void count(const std::string& reason) {
    stops_[reason] += 1;
  }

  void stop(Index c, const std::string& reason) {
    Column& col = cols_[u(c)];
    if (!col.active) {
      return;
    }
    col.active = false;
    col.next = points_[u(col.cur)];
    col.reason = reason;
    count(reason);
  }

  [[nodiscard]] bool any_active() const {
    return std::any_of(cols_.begin(), cols_.end(), [](const Column& c) { return c.active; });
  }

  void setup() {
    const std::size_t nf = in_.faces.size();
    if (in_.wall.size() != nf || in_.face_group.size() != nf ||
        in_.node_group.size() != points_.size()) {
      throw std::invalid_argument("inflate_3d: wall / face_group / node_group sizes do not match");
    }
    std::map<std::pair<Index, Index>, std::vector<Index>> edges;
    for (std::size_t f = 0; f < nf; ++f) {
      const auto& face = in_.faces[f];
      if (face.size() != 3 && face.size() != 4) {
        throw std::invalid_argument("inflate_3d: faces must have 3 or 4 nodes");
      }
      for (std::size_t i = 0; i < face.size(); ++i) {
        const Index a = face[i];
        const Index b = face[(i + 1) % face.size()];
        if (a < 0 || u(a) >= points_.size()) {
          throw std::invalid_argument("inflate_3d: face node index out of range");
        }
        edges[{std::min(a, b), std::max(a, b)}].push_back(sz(f));
      }
    }
    std::vector<Index> col_of(points_.size(), -1);
    for (std::size_t f = 0; f < nf; ++f) {
      if (!in_.wall[f]) {
        continue;
      }
      WallFace w;
      w.face = sz(f);
      const auto& face = in_.faces[f];
      for (std::size_t i = 0; i < face.size(); ++i) {
        const Index p = face[i];
        if (col_of[u(p)] < 0) {
          col_of[u(p)] = sz(cols_.size());
          Column c;
          c.base = p;
          c.cur = p;
          cols_.push_back(c);
        }
        w.cols.push_back(col_of[u(p)]);
        cols_[u(col_of[u(p)])].faces.push_back(sz(walls_.size()));
        const Index b = face[(i + 1) % face.size()];
        Index other = -1;
        for (const Index g : edges[{std::min(p, b), std::max(p, b)}]) {
          if (g != sz(f)) {
            other = g;
          }
        }
        w.across.push_back(other);
      }
      walls_.push_back(w);
    }
    if (walls_.empty()) {
      throw std::invalid_argument("inflate_3d: no wall face");
    }
    // Node kinds from the non-wall faces around each wall node.
    std::vector<std::vector<Index>> nonwall(points_.size());
    for (std::size_t f = 0; f < nf; ++f) {
      if (in_.wall[f]) {
        continue;
      }
      for (const Index p : in_.faces[f]) {
        if (col_of[u(p)] >= 0) {
          nonwall[u(p)].push_back(sz(f));
        }
      }
    }
    std::int64_t not_coplanar = 0;
    for (Column& c : cols_) {
      const auto& around = nonwall[u(c.base)];
      if (around.empty()) {
        continue;
      }
      const Index group = in_.node_group[u(c.base)];
      bool slide = group >= 0;
      Vec3 plane;
      for (const Index f : around) {
        if (in_.face_group[u(f)] != group) {
          slide = false;
        }
      }
      if (slide) {
        for (const Index f : around) {
          const Vec3 n = normalized(face_normal(f));
          if (norm(plane) == 0.0) {
            plane = n;
          } else if (norm(cross(plane, n)) > 1e-6) {
            slide = false;
          }
        }
        if (!slide) {
          ++not_coplanar;
        }
      }
      if (slide) {
        c.kind = 1;
        c.group = group;
        c.plane = plane;
      } else {
        c.kind = 2;
        c.active = false;
        c.reason = "fixed";
      }
    }
    if (not_coplanar > 0) {
      warnings_.push_back(std::to_string(not_coplanar) +
                          " sliding wall nodes are fixed: their plane's faces are not coplanar");
    }
    // Neighbours and wall spacing.
    std::set<std::pair<Index, Index>> seen;
    std::vector<double> length_sum(cols_.size(), 0.0);
    std::vector<int> length_count(cols_.size(), 0);
    for (const WallFace& w : walls_) {
      for (std::size_t i = 0; i < w.cols.size(); ++i) {
        const Index a = w.cols[i];
        const Index b = w.cols[(i + 1) % w.cols.size()];
        if (!seen.insert({std::min(a, b), std::max(a, b)}).second) {
          continue;
        }
        cols_[u(a)].nbrs.push_back(b);
        cols_[u(b)].nbrs.push_back(a);
        const double len = distance(points_[u(cols_[u(a)].base)], points_[u(cols_[u(b)].base)]);
        length_sum[u(a)] += len;
        length_sum[u(b)] += len;
        length_count[u(a)] += 1;
        length_count[u(b)] += 1;
      }
    }
    for (std::size_t c = 0; c < cols_.size(); ++c) {
      cols_[c].origin_spacing = length_sum[c] / std::max(length_count[c], 1);
      cols_[c].next = points_[u(cols_[c].base)];
    }
    // Static obstacles: the non-wall faces.
    std::vector<Tri> tris;
    for (std::size_t f = 0; f < nf; ++f) {
      if (in_.wall[f]) {
        continue;
      }
      const auto& face = in_.faces[f];
      for (std::size_t k = 1; k + 1 < face.size(); ++k) {
        Tri t;
        const std::array<Index, 3> ids = {face[0], face[k], face[k + 1]};
        for (std::size_t i = 0; i < 3; ++i) {
          t.p[i] = points_[u(ids[i])];
          t.tag[i] = ids[i];
        }
        t.group = in_.face_group[f];
        tris.push_back(t);
      }
    }
    static_ = std::make_unique<TriTree>(std::move(tris));
    for (const Column& c : cols_) {
      if (c.kind == 2) {
        count("fixed");
      }
    }
  }

  [[nodiscard]] Vec3 face_normal(Index f) const {
    std::vector<Vec3> p;
    for (const Index i : in_.faces[u(f)]) {
      p.push_back(points_[u(i)]);
    }
    return newell(p);
  }

  [[nodiscard]] Vec3 cur_pos(Index c) const {
    return points_[u(cols_[u(c)].cur)];
  }

  // Triangles of a wall face's current front (tags: current point ids).
  void front_tris(std::size_t w, std::vector<Tri>& out) const {
    const WallFace& f = walls_[w];
    for (const auto& tri : polygons(f)) {
      if (tri.size() == 3) {
        add_tri(f, w, tri[0], tri[1], tri[2], out);
      } else {
        add_tri(f, w, tri[0], tri[1], tri[2], out);
        add_tri(f, w, tri[0], tri[2], tri[3], out);
      }
    }
  }

  void add_tri(const WallFace& f, std::size_t w, std::size_t a, std::size_t b, std::size_t c,
               std::vector<Tri>& out) const {
    Tri t;
    const std::array<std::size_t, 3> local = {a, b, c};
    for (std::size_t i = 0; i < 3; ++i) {
      const Column& col = cols_[u(f.cols[local[i]])];
      t.p[i] = points_[u(col.cur)];
      t.tag[i] = col.cur;
    }
    t.owner = sz(w);
    out.push_back(t);
  }

  // Local node lists of a wall face's front polygons (one quad, one triangle or two triangles).
  static std::vector<std::vector<std::size_t>> polygons(const WallFace& f) {
    if (f.cols.size() == 3) {
      return {{0, 1, 2}};
    }
    if (!f.split) {
      return {{0, 1, 2, 3}};
    }
    if (f.diag == 0) {
      return {{0, 1, 2}, {0, 2, 3}};
    }
    return {{0, 1, 3}, {1, 2, 3}};
  }

  void march(double h) {
    // Inward normals of the current front faces.
    for (WallFace& w : walls_) {
      std::vector<Vec3> p;
      p.reserve(w.cols.size());
      for (const Index c : w.cols) {
        p.push_back(cur_pos(c));
      }
      w.inward = -normalized(newell(p));
    }
    const std::size_t n = cols_.size();
    std::vector<Vec3> normal(n);
    for (std::size_t c = 0; c < n; ++c) {
      const Column& col = cols_[c];
      if (!col.active) {
        continue;
      }
      Vec3 sum;
      for (const Index wi : col.faces) {
        const WallFace& w = walls_[u(wi)];
        const std::size_t m = w.cols.size();
        std::size_t i = 0;
        while (w.cols[i] != sz(c)) {
          ++i;
        }
        const Vec3 o = cur_pos(sz(c));
        const Vec3 a = cur_pos(w.cols[(i + 1) % m]) - o;
        const Vec3 b = cur_pos(w.cols[(i + m - 1) % m]) - o;
        const double angle = std::atan2(norm(cross(a, b)), dot(a, b));
        sum += w.inward * angle;
      }
      normal[c] = project(col, sum);
      if (visibility(col, normal[c]) < 0.05) {
        stop(sz(c), "visibility");
      }
    }
    for (int sweep = 0; sweep < opt_.normal_smoothing; ++sweep) {
      std::vector<Vec3> smoothed = normal;
      for (std::size_t c = 0; c < n; ++c) {
        const Column& col = cols_[c];
        if (!col.active) {
          continue;
        }
        Vec3 mean;
        int k = 0;
        for (const Index b : col.nbrs) {
          if (cols_[u(b)].active) {
            mean += normal[u(b)];
            ++k;
          }
        }
        if (k == 0) {
          continue;
        }
        const Vec3 candidate = project(col, normal[c] + mean * (0.5 / k));
        if (visibility(col, candidate) >= opt_.min_visibility) {
          smoothed[c] = candidate;
        }
      }
      normal = smoothed;
    }
    for (std::size_t c = 0; c < n; ++c) {
      Column& col = cols_[c];
      if (!col.active) {
        continue;
      }
      col.normal = normal[c];
      col.scale = 1.0 / std::max(visibility(col, col.normal), 1.0 / opt_.max_scale);
      double spacing = 0.0;
      for (const Index b : col.nbrs) {
        spacing += distance(cur_pos(sz(c)), cur_pos(b));
      }
      spacing /= static_cast<double>(std::max<std::size_t>(col.nbrs.size(), 1));
      col.target = std::max(spacing, col.origin_spacing);
      if (opt_.max_size) {
        col.target = std::min(col.target, *opt_.max_size);
      }
      if (opt_.transition && h >= opt_.transition_ratio * col.target) {
        stop(sz(c), "transition");
        continue;
      }
      col.next = cur_pos(sz(c)) + col.normal * (h * col.scale);
    }
    // Proximity: a facing front or boundary face ahead closer than proximity * size.
    std::vector<Tri> tris;
    for (std::size_t w = 0; w < walls_.size(); ++w) {
      front_tris(w, tris);
    }
    current_ = std::make_unique<TriTree>(std::move(tris));
    for (std::size_t c = 0; c < n; ++c) {
      Column& col = cols_[c];
      if (!col.active) {
        continue;
      }
      const Vec3 a = cur_pos(sz(c));
      const double step = h * col.scale;
      const double reach = step + opt_.proximity * col.target;
      const Vec3 b = a + col.normal * reach;
      const std::array<Index, 3> tags = {col.cur, 0, 0};
      double nearest = std::numeric_limits<double>::infinity();
      auto test = [&](const Tri& t) {
        if (share_tag(t, tags, 1) || (col.kind == 1 && t.group == col.group)) {
          return;
        }
        const Vec3 out = cross(t.p[1] - t.p[0], t.p[2] - t.p[0]);
        if (dot(out, col.normal) <= 0.0) {
          return;
        }
        double param = 0.0;
        if (segment_hits(a, b, t, &param)) {
          nearest = std::min(nearest, param * reach);
        }
      };
      static_->query(min3(a, b), max3(a, b), test);
      current_->query(min3(a, b), max3(a, b), test);
      if (nearest - step < opt_.proximity * col.target) {
        stop(sz(c), "proximity");
      }
    }
  }

  [[nodiscard]] static Vec3 project(const Column& col, const Vec3& v) {
    if (col.kind == 1) {
      return normalized(v - col.plane * dot(v, col.plane));
    }
    return normalized(v);
  }

  [[nodiscard]] double visibility(const Column& col, const Vec3& n) const {
    double worst = 1.0;
    for (const Index wi : col.faces) {
      worst = std::min(worst, dot(n, walls_[u(wi)].inward));
    }
    return worst;
  }

  // Elements of wall face w's layer cell for the current active flags; new column nodes get ids
  // top_id[c], centroids ids from next_extra (positions appended to extra).
  void build_cells(std::size_t w, const std::vector<Index>& top_id, Index* next_extra,
                   std::vector<Vec3>& extra, std::vector<Element>& out, bool* split,
                   int* diag) const {
    const WallFace& f = walls_[w];
    const std::size_t m = f.cols.size();
    std::vector<Index> bottom(m);
    std::vector<Index> top(m);
    std::vector<bool> act(m);
    int active = 0;
    for (std::size_t i = 0; i < m; ++i) {
      const Column& col = cols_[u(f.cols[i])];
      bottom[i] = col.cur;
      act[i] = col.active;
      top[i] = col.cur;
      if (col.active) {
        top[i] = top_id[u(f.cols[i])];
        ++active;
      }
    }
    *split = f.split;
    *diag = f.diag;
    if (active == 0) {
      return;
    }
    auto tri_cell = [&](std::size_t i0, std::size_t i1, std::size_t i2) {
      const std::array<std::size_t, 3> l = {i0, i1, i2};
      int on = 0;
      for (const std::size_t i : l) {
        on += static_cast<int>(act[i]);
      }
      if (on == 0) {
        return;
      }
      Element e;
      if (on == 3) {
        e.kind = Kind::prism;
        e.v = {bottom[i0], bottom[i2], bottom[i1], top[i0], top[i2], top[i1], 0, 0};
        e.full = true;
        out.push_back(e);
        return;
      }
      if (on == 2) {
        std::size_t s = 0;
        while (act[l[s]]) {
          ++s;
        }
        const std::size_t a = l[(s + 1) % 3];
        const std::size_t b = l[(s + 2) % 3];
        e.kind = Kind::pyramid;
        e.v = {bottom[a], bottom[b], top[b], top[a], bottom[l[s]], 0, 0, 0};
        out.push_back(e);
        return;
      }
      std::size_t t = 0;
      while (!act[l[t]]) {
        ++t;
      }
      e.kind = Kind::tet;
      e.v = {bottom[i0], bottom[i2], bottom[i1], top[l[t]], 0, 0, 0, 0};
      out.push_back(e);
    };
    if (m == 3) {
      tri_cell(0, 1, 2);
      return;
    }
    if (f.split) {
      if (f.diag == 0) {
        tri_cell(0, 1, 2);
        tri_cell(0, 2, 3);
      } else {
        tri_cell(0, 1, 3);
        tri_cell(1, 2, 3);
      }
      return;
    }
    Element e;
    if (active == 4) {
      e.kind = Kind::hex;
      e.v = {bottom[0], bottom[3], bottom[2], bottom[1], top[0], top[3], top[2], top[1]};
      e.full = true;
      out.push_back(e);
      return;
    }
    // From the first partial layer on, the face continues as two triangles: its front never
    // keeps a steep quad (a step between stopped and advancing columns), whose pyramid would
    // crowd its neighbours.
    *split = true;
    if (active == 1) {
      std::size_t t = 0;
      while (!act[t]) {
        ++t;
      }
      *diag = static_cast<int>(t % 2);
      e.kind = Kind::pyramid;
      e.v = {bottom[0], bottom[3], bottom[2], bottom[1], top[t], 0, 0, 0};
      out.push_back(e);
      return;
    }
    // One or two stopped columns: split from the centroid. The top is split along a diagonal
    // between two advancing columns where there is one, so one top triangle stays level (and
    // continues as prisms) instead of both sloping down to a stopped column.
    std::size_t s = 0;
    while (act[s]) {
      ++s;
    }
    *diag = static_cast<int>(s % 2);
    if (act[0] && act[2]) {
      *diag = 0;
    } else if (act[1] && act[3]) {
      *diag = 1;
    }
    std::vector<std::vector<Index>> faces;
    faces.push_back({bottom[0], bottom[1], bottom[2], bottom[3]});
    for (std::size_t i = 0; i < 4; ++i) {
      const std::size_t j = (i + 1) % 4;
      if (!act[i] && !act[j]) {
        continue;
      }
      if (!act[i]) {
        faces.push_back({bottom[j], bottom[i], top[j]});
      } else if (!act[j]) {
        faces.push_back({bottom[j], bottom[i], top[i]});
      } else {
        faces.push_back({bottom[j], bottom[i], top[i], top[j]});
      }
    }
    if (*diag == 0) {
      faces.push_back({top[0], top[2], top[1]});
      faces.push_back({top[0], top[3], top[2]});
    } else {
      faces.push_back({top[0], top[3], top[1]});
      faces.push_back({top[1], top[3], top[2]});
    }
    std::set<Index> distinct;
    for (const auto& face : faces) {
      distinct.insert(face.begin(), face.end());
    }
    Vec3 centre;
    for (const Index id : distinct) {
      centre += position(id, top_id, extra);
    }
    centre = centre / static_cast<double>(distinct.size());
    const Index cid = (*next_extra)++;
    extra.push_back(centre);
    for (const auto& face : faces) {
      Element sub;
      if (face.size() == 4) {
        sub.kind = Kind::pyramid;
        sub.v = {face[0], face[3], face[2], face[1], cid, 0, 0, 0};
      } else {
        sub.kind = Kind::tet;
        sub.v = {face[0], face[2], face[1], cid, 0, 0, 0, 0};
      }
      out.push_back(sub);
    }
  }

  // Position of a point id during validation: existing points, pending column nodes
  // (pending_[id - points_.size()]), then centroids (extra).
  [[nodiscard]] Vec3 position(Index id, const std::vector<Index>& /*top_id*/,
                              const std::vector<Vec3>& extra) const {
    if (u(id) < points_.size()) {
      return points_[u(id)];
    }
    const std::size_t k = u(id) - points_.size();
    if (k < pending_.size()) {
      return pending_[k];
    }
    return extra[k - pending_.size()];
  }

  [[nodiscard]] bool element_ok(const Element& e, const std::vector<Index>& top_id,
                                const std::vector<Vec3>& extra) const {
    auto p = [&](int i) { return position(e.v[u(i)], top_id, extra); };
    constexpr double kMin = 1e-9;
    if (e.kind == Kind::tet) {
      return corner(p(0), p(1), p(2), p(3)) > kMin;
    }
    if (e.kind == Kind::pyramid) {
      for (int k = 0; k < 4; ++k) {
        if (corner(p(k), p((k + 1) % 4), p((k + 3) % 4), p(4)) <= kMin) {
          return false;
        }
      }
      return true;
    }
    double limit = kMin;
    if (e.full) {
      limit = opt_.min_jacobian;
    }
    if (e.kind == Kind::prism) {
      for (const auto& c : kPrismCorners) {
        if (corner(p(c[0]), p(c[1]), p(c[2]), p(c[3])) <= limit) {
          return false;
        }
      }
      return true;
    }
    for (const auto& c : kHexCorners) {
      if (corner(p(c[0]), p(c[1]), p(c[2]), p(c[3])) <= limit) {
        return false;
      }
    }
    return true;
  }

  // Pending node ids of active columns (in column order) and their positions.
  std::vector<Index> assign_pending() {
    std::vector<Index> top_id(cols_.size(), -1);
    pending_.clear();
    for (std::size_t c = 0; c < cols_.size(); ++c) {
      if (cols_[c].active) {
        top_id[c] = sz(points_.size() + pending_.size());
        pending_.push_back(cols_[c].next);
      }
    }
    return top_id;
  }

  void slope() {
    bool stopped = true;
    while (stopped) {
      stopped = false;
      for (std::size_t c = 0; c < cols_.size(); ++c) {
        const Column& col = cols_[c];
        if (!col.active) {
          continue;
        }
        const Vec3 here = cur_pos(sz(c));
        for (const Index b : col.nbrs) {
          const Column& other = cols_[u(b)];
          if (other.active) {
            continue;
          }
          // Lateral gap: offset to the neighbour's top without its component along the normal.
          const Vec3 offset = points_[u(other.cur)] - here;
          const double along = dot(offset, col.normal);
          const double gap = norm(offset - col.normal * along);
          if (col.height + h_ - other.height > std::max(opt_.max_slope * gap, h_) * (1.0 + 1e-9)) {
            stop(sz(c), "slope");
            stopped = true;
            break;
          }
        }
      }
    }
  }

  void validate(double h) {
    h_ = h;
    bool changed = true;
    while (changed) {
      changed = false;
      slope();
      // Cells.
      std::vector<Index> top_id = assign_pending();
      std::vector<Vec3> extra;
      Index next_extra = sz(points_.size() + pending_.size());
      for (std::size_t w = 0; w < walls_.size(); ++w) {
        std::vector<Element> cells;
        bool split = false;
        int diag = 0;
        build_cells(w, top_id, &next_extra, extra, cells, &split, &diag);
        bool ok = true;
        for (const Element& e : cells) {
          ok = ok && element_ok(e, top_id, extra);
        }
        if (!ok) {
          // Stop the face's worst column(s): the most tilted marching normals (those cells are
          // re-checked in the next round).
          double worst = std::numeric_limits<double>::infinity();
          for (const Index c : walls_[w].cols) {
            const Column& col = cols_[u(c)];
            if (col.active) {
              worst = std::min(worst, visibility(col, col.normal));
            }
          }
          for (const Index c : walls_[w].cols) {
            const Column& col = cols_[u(c)];
            if (col.active && visibility(col, col.normal) <= worst + 1e-12) {
              stop(c, "quality");
              changed = true;
            }
          }
        }
      }
      if (changed) {
        continue;
      }
      // Collisions: new column edges and new front faces against the current fronts, the
      // non-wall faces and (column edges) the new front faces.
      std::vector<Tri> fresh;
      for (std::size_t w = 0; w < walls_.size(); ++w) {
        const WallFace& f = walls_[w];
        bool moving = false;
        for (const Index c : f.cols) {
          moving = moving || cols_[u(c)].active;
        }
        if (!moving) {
          continue;
        }
        for (const auto& poly : polygons(f)) {
          for (std::size_t k = 1; k + 1 < poly.size(); ++k) {
            Tri t;
            const std::array<std::size_t, 3> l = {poly[0], poly[k], poly[k + 1]};
            for (std::size_t i = 0; i < 3; ++i) {
              const Index c = f.cols[l[i]];
              const Column& col = cols_[u(c)];
              t.p[i] = col.next;
              t.tag[i] = col.cur;
              if (col.active) {
                t.tag[i] = -(c + 1);
              }
            }
            t.owner = sz(w);
            fresh.push_back(t);
          }
        }
      }
      for (const Tri& t : fresh) {
        Index slide_group = -1;
        for (const Index c : walls_[u(t.owner)].cols) {
          if (cols_[u(c)].kind == 1) {
            slide_group = cols_[u(c)].group;
          }
        }
        bool hit = false;
        Vec3 lo = min3(min3(t.p[0], t.p[1]), t.p[2]);
        Vec3 hi = max3(max3(t.p[0], t.p[1]), t.p[2]);
        auto test = [&](const Tri& o) {
          if (hit || o.owner == t.owner || share_tag(o, t.tag, 3) ||
              (slide_group >= 0 && o.group == slide_group)) {
            return;
          }
          if (triangles_cross(t, o)) {
            hit = true;
          }
        };
        static_->query(lo, hi, test);
        current_->query(lo, hi, test);
        if (hit) {
          for (const Index c : walls_[u(t.owner)].cols) {
            if (cols_[u(c)].active) {
              stop(c, "collision");
              changed = true;
            }
          }
        }
      }
      if (changed) {
        continue;
      }
      const TriTree fresh_tree(fresh);
      for (std::size_t c = 0; c < cols_.size(); ++c) {
        const Column& col = cols_[c];
        if (!col.active) {
          continue;
        }
        const Vec3 a = cur_pos(sz(c));
        const Vec3 b = col.next;
        const std::array<Index, 3> tags = {col.cur, -(sz(c) + 1), 0};
        bool hit = false;
        auto test = [&](const Tri& o) {
          if (hit || share_tag(o, tags, 2) || (col.kind == 1 && o.group == col.group)) {
            return;
          }
          double param = 0.0;
          if (segment_hits(a, b, o, &param)) {
            hit = true;
          }
        };
        static_->query(min3(a, b), max3(a, b), test);
        current_->query(min3(a, b), max3(a, b), test);
        fresh_tree.query(min3(a, b), max3(a, b), test);
        if (hit) {
          stop(sz(c), "collision");
          changed = true;
        }
      }
    }
  }

  void commit(double h) {
    std::vector<Index> top_id = assign_pending();
    if (pending_.empty()) {
      return;
    }
    std::vector<Vec3> extra;
    Index next_extra = sz(points_.size() + pending_.size());
    std::vector<Element> cells;
    std::vector<std::pair<bool, int>> modes(walls_.size());
    for (std::size_t w = 0; w < walls_.size(); ++w) {
      bool split = false;
      int diag = 0;
      build_cells(w, top_id, &next_extra, extra, cells, &split, &diag);
      modes[w] = {split, diag};
      // Lateral faces on sliding planes (edges between a wall face and a non-wall face).
      const WallFace& f = walls_[w];
      const std::size_t m = f.cols.size();
      for (std::size_t i = 0; i < m; ++i) {
        const Index other = f.across[i];
        if (other < 0 || in_.wall[u(other)]) {
          continue;
        }
        const Column& a = cols_[u(f.cols[i])];
        const Column& b = cols_[u(f.cols[(i + 1) % m])];
        if (!a.active && !b.active) {
          continue;
        }
        std::vector<Index> side;
        side.push_back(b.cur);
        side.push_back(a.cur);
        if (a.active) {
          side.push_back(top_id[u(f.cols[i])]);
        }
        if (b.active) {
          side.push_back(top_id[u(f.cols[(i + 1) % m])]);
        }
        sides_.push_back(side);
        side_face_.push_back(other);
      }
    }
    for (const Vec3& p : pending_) {
      points_.push_back(p);
    }
    for (const Vec3& p : extra) {
      points_.push_back(p);
    }
    for (const Element& e : cells) {
      if (e.kind == Kind::tet) {
        tets_.push_back({e.v[0], e.v[1], e.v[2], e.v[3]});
      } else if (e.kind == Kind::pyramid) {
        pyramids_.push_back({e.v[0], e.v[1], e.v[2], e.v[3], e.v[4]});
      } else if (e.kind == Kind::prism) {
        prisms_.push_back({e.v[0], e.v[1], e.v[2], e.v[3], e.v[4], e.v[5]});
      } else {
        hexes_.push_back(e.v);
      }
    }
    for (std::size_t w = 0; w < walls_.size(); ++w) {
      walls_[w].split = modes[w].first;
      walls_[w].diag = modes[w].second;
    }
    for (std::size_t c = 0; c < cols_.size(); ++c) {
      Column& col = cols_[c];
      if (!col.active) {
        continue;
      }
      col.cur = top_id[c];
      col.height += h;
      col.layers += 1;
      max_layers_ = std::max(max_layers_, col.layers);
    }
  }

  Inflation3DResult finish() {
    Inflation3DResult r;
    for (const Vec3& p : points_) {
      r.points.push_back(p.to_array());
    }
    r.tets = tets_;
    r.pyramids = pyramids_;
    r.prisms = prisms_;
    r.hexes = hexes_;
    for (const WallFace& f : walls_) {
      for (const auto& poly : polygons(f)) {
        std::vector<Index> ids;
        ids.reserve(poly.size());
        for (const std::size_t i : poly) {
          ids.push_back(cols_[u(f.cols[i])].cur);
        }
        r.front.push_back(ids);
        r.front_face.push_back(f.face);
      }
    }
    r.sides = sides_;
    r.side_face = side_face_;
    r.top.resize(in_.points.size());
    std::iota(r.top.begin(), r.top.end(), Index{0});
    r.height.assign(in_.points.size(), 0.0);
    r.layers.assign(in_.points.size(), 0);
    r.reason.assign(in_.points.size(), "");
    for (const Column& c : cols_) {
      r.reason[u(c.base)] = c.reason;
      if (c.active) {
        r.reason[u(c.base)] = "layers";
      }
      r.top[u(c.base)] = c.cur;
      r.height[u(c.base)] = c.height;
      r.layers[u(c.base)] = c.layers;
    }
    r.max_layers = max_layers_;
    r.stops = stops_;
    if (stops_.count("visibility") > 0) {
      r.warnings.push_back(std::to_string(stops_["visibility"]) +
                           " columns stopped: marching normal not visible from all their faces "
                           "(sharp ridges or corners)");
    }
    for (const std::string& w : warnings_) {
      r.warnings.push_back(w);
    }
    return r;
  }

  const Inflation3DInput& in_;
  const Inflation3DOptions& opt_;
  std::vector<Vec3> points_;
  std::vector<Column> cols_;
  std::vector<WallFace> walls_;
  std::unique_ptr<TriTree> static_;
  std::unique_ptr<TriTree> current_;
  std::vector<Vec3> pending_;
  double h_ = 0.0;
  std::vector<std::array<Index, 4>> tets_;
  std::vector<std::array<Index, 5>> pyramids_;
  std::vector<std::array<Index, 6>> prisms_;
  std::vector<std::array<Index, 8>> hexes_;
  std::vector<std::vector<Index>> sides_;
  std::vector<Index> side_face_;
  int max_layers_ = 0;
  std::map<std::string, std::int64_t> stops_;
  std::vector<std::string> warnings_;
};

}  // namespace

Inflation3DResult inflate_3d(const Inflation3DInput& input, const Inflation3DOptions& options) {
  if (!(options.first_height > 0.0) || !(options.growth_rate >= 1.0) || options.max_layers < 1 ||
      (options.total_thickness && !(*options.total_thickness > 0.0)) ||
      !(options.transition_ratio > 0.0) || (options.max_size && !(*options.max_size > 0.0)) ||
      !(options.proximity >= 0.0) || !(options.max_slope > 0.0) ||
      !(options.min_jacobian >= 0.0 && options.min_jacobian < 1.0) || !(options.max_scale >= 1.0) ||
      !(options.min_visibility > 0.0 && options.min_visibility < 1.0) ||
      options.normal_smoothing < 0) {
    throw std::invalid_argument("inflate_3d: invalid options");
  }
  Inflator inflator(input, options);
  return inflator.run();
}

}  // namespace og
