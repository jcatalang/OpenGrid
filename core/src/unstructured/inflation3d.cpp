#include "opengrid/core/unstructured/inflation3d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <numbers>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
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

constexpr double kDegree = std::numbers::pi / 180.0;

double angle_between(const Vec3& a, const Vec3& b) {
  return std::atan2(norm(cross(a, b)), dot(a, b));
}

// Spherical interpolation between unit vectors (normalised linear interpolation when they are
// (nearly) opposite).
Vec3 slerp(const Vec3& a, const Vec3& b, double t) {
  const double w = angle_between(a, b);
  if (w < 1e-12) {
    return a;
  }
  if (std::sin(w) < 1e-6) {
    return normalized(a * (1.0 - t) + b * t);
  }
  return a * (std::sin((1.0 - t) * w) / std::sin(w)) + b * (std::sin(t * w) / std::sin(w));
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
  bool fan = false;          // a fan column at a ridge (several columns share its base node)
  Vec3 start_normal;         // fan columns: the normal of the first layer; others: their sector's
};

struct WallFace {
  Index face = 0;
  std::vector<Index> cols;
  bool split = false;         // a quad continuing as two triangles
  int diag = 0;               // 0: diagonal (0, 2); 1: (1, 3)
  std::vector<Index> across;  // per edge (i, i + 1): the non-wall face across it, or -1
  Vec3 inward;
  bool virt = false;  // a fan / cap face (zero area at the wall)
  bool cap = false;   // a cap triangle at a ridge corner: its cells need only be valid
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
    // Columns sharing a node (fans) stop together in their first layer: a part of a fan left on
    // the wall would pinch the front against the ridge edge.
    if (col.layers == 0) {
      for (const Index d : node_cols_[u(col.base)]) {
        stop(d, reason);
      }
    }
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
    // Directed wall edges -> the wall face traversing them.
    std::map<std::pair<Index, Index>, Index> directed;
    std::vector<bool> wall_node(points_.size(), false);
    for (std::size_t f = 0; f < nf; ++f) {
      if (!in_.wall[f]) {
        continue;
      }
      const auto& face = in_.faces[f];
      for (std::size_t i = 0; i < face.size(); ++i) {
        directed[{face[i], face[(i + 1) % face.size()]}] = sz(f);
        wall_node[u(face[i])] = true;
      }
    }
    if (directed.empty()) {
      throw std::invalid_argument("inflate_3d: no wall face");
    }
    // Node kinds from the non-wall faces around each wall node.
    std::vector<std::vector<Index>> nonwall(points_.size());
    for (std::size_t f = 0; f < nf; ++f) {
      if (in_.wall[f]) {
        continue;
      }
      for (const Index p : in_.faces[f]) {
        if (wall_node[u(p)]) {
          nonwall[u(p)].push_back(sz(f));
        }
      }
    }
    std::vector<int> kind(points_.size(), 0);
    std::vector<Vec3> plane(points_.size());
    std::int64_t not_coplanar = 0;
    for (std::size_t p = 0; p < points_.size(); ++p) {
      const auto& around = nonwall[p];
      if (around.empty()) {
        continue;
      }
      const Index group = in_.node_group[p];
      bool slide = group >= 0;
      for (const Index f : around) {
        if (in_.face_group[u(f)] != group) {
          slide = false;
        }
      }
      if (slide) {
        for (const Index f : around) {
          const Vec3 n = normalized(face_normal(f));
          if (norm(plane[p]) == 0.0) {
            plane[p] = n;
          } else if (norm(cross(plane[p], n)) > 1e-6) {
            slide = false;
          }
        }
        if (!slide) {
          ++not_coplanar;
        }
      }
      kind[p] = 2;
      if (slide) {
        kind[p] = 1;
      }
    }
    if (not_coplanar > 0) {
      warnings_.push_back(std::to_string(not_coplanar) +
                          " sliding wall nodes are fixed: their plane's faces are not coplanar");
    }
    // Fluid-side unit normals and centroids of the wall faces.
    std::vector<Vec3> fluid(nf);
    std::vector<Vec3> centroid(nf);
    for (std::size_t f = 0; f < nf; ++f) {
      if (!in_.wall[f]) {
        continue;
      }
      fluid[f] = -normalized(face_normal(sz(f)));
      for (const Index p : in_.faces[f]) {
        centroid[f] += points_[u(p)];
      }
      centroid[f] = centroid[f] / static_cast<double>(in_.faces[f].size());
    }
    // Ridges: convex wall edges whose fluid-side normals turn by more than fan_threshold.
    const double threshold = opt_.fan_threshold_deg * kDegree;
    std::set<std::pair<Index, Index>> ridge;  // (a, b) as traversed by the first face
    for (const auto& [e, f1] : directed) {
      const auto it = directed.find({e.second, e.first});
      if (it == directed.end() || e.first > e.second) {
        continue;
      }
      const Index f2 = it->second;
      const double turn = angle_between(fluid[u(f1)], fluid[u(f2)]);
      const bool convex = dot(centroid[u(f2)] - centroid[u(f1)], fluid[u(f1)]) < 0.0;
      if (convex && turn > threshold && opt_.fan_threshold_deg < 180.0) {
        ridge.insert(e);
      }
    }
    auto is_ridge = [&](Index a, Index b) {
      return ridge.count({a, b}) > 0 || ridge.count({b, a}) > 0;
    };
    // Sectors: the wall faces around a node, connected across its non-ridge edges.
    std::map<std::pair<Index, Index>, Index> sector_root;  // (node, face) -> root face
    std::vector<std::vector<Index>> node_faces(points_.size());
    for (std::size_t f = 0; f < nf; ++f) {
      if (in_.wall[f]) {
        for (const Index p : in_.faces[f]) {
          node_faces[u(p)].push_back(sz(f));
        }
      }
    }
    for (std::size_t p = 0; p < points_.size(); ++p) {
      // Union-find over the faces of node p (indices into node_faces[p]).
      const auto& around = node_faces[p];
      std::vector<std::size_t> up(around.size());
      std::iota(up.begin(), up.end(), std::size_t{0});
      auto find = [&](std::size_t i) {
        while (up[i] != i) {
          up[i] = up[up[i]];
          i = up[i];
        }
        return i;
      };
      for (std::size_t i = 0; i < around.size(); ++i) {
        const auto& face = in_.faces[u(around[i])];
        std::size_t k = 0;
        while (face[k] != sz(p)) {
          ++k;
        }
        const Index x = face[(k + 1) % face.size()];
        const auto it = directed.find({x, sz(p)});
        if (it == directed.end() || is_ridge(sz(p), x)) {
          continue;
        }
        const auto j = static_cast<std::size_t>(
            std::find(around.begin(), around.end(), it->second) - around.begin());
        const std::size_t ri = find(i);
        const std::size_t rj = find(j);
        // The smaller face index is the root, so sectors are numbered deterministically.
        if (around[ri] < around[rj]) {
          up[rj] = ri;
        } else if (ri != rj) {
          up[ri] = rj;
        }
      }
      for (std::size_t i = 0; i < around.size(); ++i) {
        sector_root[{sz(p), around[i]}] = around[find(i)];
      }
    }
    // Sector normals (angle-weighted fluid-side normals of the sector's faces at the node).
    std::map<std::pair<Index, Index>, Vec3> sector_normal;
    for (const auto& [key, root] : sector_root) {
      const auto& face = in_.faces[u(key.second)];
      const std::size_t m = face.size();
      std::size_t k = 0;
      while (face[k] != key.first) {
        ++k;
      }
      const Vec3 o = points_[u(key.first)];
      const Vec3 a = points_[u(face[(k + 1) % m])] - o;
      const Vec3 b = points_[u(face[(k + m - 1) % m])] - o;
      sector_normal[{key.first, root}] += fluid[u(key.second)] * angle_between(a, b);
    }
    auto column_normal = [&](Index p, const Vec3& n) {
      if (kind[u(p)] == 1) {
        return normalized(n - plane[u(p)] * dot(n, plane[u(p)]));
      }
      return normalized(n);
    };
    // Columns: one per sector of every wall node (one per node without ridges).
    std::map<std::pair<Index, Index>, Index> col_of;  // (node, sector root) -> column
    auto new_column = [&](Index p) {
      Column c;
      c.base = p;
      c.cur = p;
      c.kind = kind[u(p)];
      if (c.kind == 1) {
        c.group = in_.node_group[u(p)];
        c.plane = plane[u(p)];
      } else if (c.kind == 2) {
        c.active = false;
        c.reason = "fixed";
      }
      cols_.push_back(c);
      return sz(cols_.size() - 1);
    };
    std::vector<std::vector<Index>> face_cols;  // per real wall face (walls_ order)
    std::vector<Index> real_faces;
    for (std::size_t f = 0; f < nf; ++f) {
      if (!in_.wall[f]) {
        continue;
      }
      std::vector<Index> list;
      for (const Index p : in_.faces[f]) {
        const Index root = sector_root[{p, sz(f)}];
        auto it = col_of.find({p, root});
        if (it == col_of.end()) {
          it = col_of.emplace(std::make_pair(p, root), new_column(p)).first;
          cols_[u(it->second)].start_normal = column_normal(p, sector_normal[{p, root}]);
        }
        list.push_back(it->second);
      }
      face_cols.push_back(list);
      real_faces.push_back(sz(f));
    }
    // Fans: along every ridge edge, columns at both ends with normals turning from one side's
    // sector to the other's. A fan at node p is keyed by (p, the two sectors in increasing
    // order), so a ridge line shares its fans between consecutive edges whichever way they are
    // traversed; a chain of ridge edges sharing fans uses one count of fan columns (the largest
    // turn needs at most fan_step degrees per step). A ridge ending at a node (one sector) has
    // its own fan there.
    struct FanKey {
      Index p;
      Index first;
      Index second;
      bool operator<(const FanKey& o) const {
        return std::tie(p, first, second) < std::tie(o.p, o.first, o.second);
      }
    };
    std::vector<std::pair<Index, Index>> ridges(ridge.begin(), ridge.end());
    std::vector<std::array<FanKey, 2>> ridge_keys;
    std::map<FanKey, std::vector<std::size_t>> key_ridges;
    for (std::size_t r = 0; r < ridges.size(); ++r) {
      const auto [a, b] = ridges[r];
      const Index f1 = directed[{a, b}];
      const Index f2 = directed[{b, a}];
      std::array<FanKey, 2> keys{};
      const std::array<Index, 2> ends = {a, b};
      for (std::size_t k = 0; k < 2; ++k) {
        const Index x = sector_root[{ends[k], f1}];
        const Index y = sector_root[{ends[k], f2}];
        keys[k] = {ends[k], std::min(x, y), std::max(x, y)};
        key_ridges[keys[k]].push_back(r);
      }
      ridge_keys.push_back(keys);
    }
    std::vector<std::size_t> chain(ridges.size());
    std::iota(chain.begin(), chain.end(), std::size_t{0});
    auto find_chain = [&](std::size_t i) {
      while (chain[i] != i) {
        chain[i] = chain[chain[i]];
        i = chain[i];
      }
      return i;
    };
    for (const auto& [key, list] : key_ridges) {
      for (const std::size_t r : list) {
        const std::size_t x = find_chain(r);
        const std::size_t y = find_chain(list.front());
        chain[std::max(x, y)] = std::min(x, y);
      }
    }
    // Turn of a fan: between the two sectors' normals, or (one sector: a ridge ending at the
    // node) between the normals of the two faces of the ridge edge.
    auto fan_ends = [&](const FanKey& key, std::size_t r) {
      if (key.first != key.second) {
        return std::make_pair(sector_normal[{key.p, key.first}],
                              sector_normal[{key.p, key.second}]);
      }
      const auto [a, b] = ridges[r];
      return std::make_pair(fluid[u(directed[{a, b}])], fluid[u(directed[{b, a}])]);
    };

    const double step = std::max(opt_.fan_step_deg, 1.0) * kDegree;
    std::map<std::size_t, int> chain_fans;
    for (std::size_t r = 0; r < ridges.size(); ++r) {
      for (const FanKey& key : ridge_keys[r]) {
        const auto [n0, n1] = fan_ends(key, r);
        const int k = std::max(1, static_cast<int>(std::ceil(angle_between(n0, n1) / step)) - 1);
        int& fans = chain_fans[find_chain(r)];
        fans = std::max(fans, k);
      }
    }
    std::map<FanKey, std::vector<Index>> fan_cols;  // the fan's columns, first to second sector
    for (std::size_t r = 0; r < ridges.size(); ++r) {
      for (const FanKey& key : ridge_keys[r]) {
        if (fan_cols.count(key) > 0) {
          continue;
        }
        const int k = chain_fans[find_chain(r)];
        const auto [n0, n1] = fan_ends(key, r);
        std::vector<Index> list = {col_of[{key.p, key.first}]};
        if (key.first == key.second) {
          // A ridge ending at the node (one sector): its single column cannot serve both sides
          // of the ridge, so it does not grow and the fan narrows to it (the layers taper
          // towards the ridge end).
          Column& end = cols_[u(list.front())];
          if (end.active) {
            end.active = false;
            end.reason = "ridge end";
            count("ridge end");
          }
          fan_cols[key] = std::vector<Index>(static_cast<std::size_t>(k) + 2, list.front());
          continue;
        }
        for (int i = 1; i <= k; ++i) {
          const Index c = new_column(key.p);
          cols_[u(c)].fan = true;
          cols_[u(c)].start_normal = column_normal(
              key.p, slerp(normalized(n0), normalized(n1), static_cast<double>(i) / (k + 1)));
          list.push_back(c);
        }
        list.push_back(col_of[{key.p, key.second}]);
        fan_cols[key] = list;
      }
    }
    // The fan's columns from the side of the ridge edge's first face to its second face's.
    auto fan_from_first = [&](const FanKey& key, std::size_t r) {
      std::vector<Index> list = fan_cols.at(key);
      const auto [a, b] = ridges[r];
      if (key.first != key.second && sector_root[{key.p, directed[{a, b}]}] != key.first) {
        std::reverse(list.begin(), list.end());
      }
      return list;
    };
    // Wall faces: the real ones (as before), then the fan quads along every ridge edge.
    std::vector<std::vector<Index>> virt_cols;
    std::vector<Index> virt_face;
    for (std::size_t r = 0; r < ridges.size(); ++r) {
      const std::vector<Index> sa = fan_from_first(ridge_keys[r][0], r);
      const std::vector<Index> sb = fan_from_first(ridge_keys[r][1], r);
      for (std::size_t i = 0; i + 1 < sa.size(); ++i) {
        // Without repeated columns (a fan narrowing to a ridge end becomes triangles).
        std::vector<Index> list;
        for (const Index c : {sb[i], sa[i], sa[i + 1], sb[i + 1]}) {
          if (list.empty() || list.back() != c) {
            list.push_back(c);
          }
        }
        if (list.size() > 1 && list.front() == list.back()) {
          list.pop_back();
        }
        if (list.size() >= 3) {
          virt_cols.push_back(list);
          virt_face.push_back(directed[{ridges[r].first, ridges[r].second}]);
        }
      }
    }
    const std::size_t n_fan_quads = virt_cols.size();
    // Caps: at interior nodes (no non-wall face) the fan quads leave open edges between columns
    // of the node; they form cycles, closed by triangles around a sector column (a ridge
    // ending at the node) or around a new centre column (a corner of several ridges).
    std::map<Index, std::vector<std::pair<Index, Index>>> open;  // node -> directed column edges
    {
      std::set<std::pair<Index, Index>> used;
      for (const auto& list : virt_cols) {
        for (std::size_t i = 0; i < list.size(); ++i) {
          const Index x = list[i];
          const Index y = list[(i + 1) % list.size()];
          if (cols_[u(x)].base == cols_[u(y)].base) {
            used.insert({x, y});
          }
        }
      }
      for (const auto& e : used) {
        if (used.count({e.second, e.first}) == 0) {
          open[cols_[u(e.first)].base].push_back(e);
        }
      }
    }
    for (const auto& [p, list] : open) {
      if (!nonwall[u(p)].empty()) {
        continue;  // a boundary node: the open edges lie on its non-wall faces
      }
      std::map<Index, Index> next;
      for (const auto& [x, y] : list) {
        if (!next.emplace(x, y).second) {
          throw std::logic_error("inflate_3d: fan columns at a node do not form a cycle");
        }
      }
      std::set<Index> done;
      for (const auto& [x0, y0] : next) {
        if (done.count(x0) > 0) {
          continue;
        }
        std::vector<Index> cycle;
        Index x = x0;
        while (done.count(x) == 0) {
          done.insert(x);
          cycle.push_back(x);
          const auto it = next.find(x);
          if (it == next.end()) {
            throw std::logic_error("inflate_3d: fan columns at a node do not form a cycle");
          }
          x = it->second;
        }
        // The cap traverses the cycle backwards.
        std::reverse(cycle.begin(), cycle.end());
        std::vector<std::size_t> sector;
        for (std::size_t i = 0; i < cycle.size(); ++i) {
          if (!cols_[u(cycle[i])].fan) {
            sector.push_back(i);
          }
        }
        const Index face = node_faces[u(p)].front();
        if (sector.size() == 1) {
          std::rotate(cycle.begin(), cycle.begin() + static_cast<std::ptrdiff_t>(sector.front()),
                      cycle.end());
          for (std::size_t i = 1; i + 1 < cycle.size(); ++i) {
            virt_cols.push_back({cycle[0], cycle[i], cycle[i + 1]});
            virt_face.push_back(face);
          }
          continue;
        }
        Vec3 mean;
        for (const Index c : cycle) {
          mean += cols_[u(c)].start_normal;
        }
        const Index centre = new_column(p);
        cols_[u(centre)].fan = true;
        cols_[u(centre)].start_normal = column_normal(p, mean);
        for (std::size_t i = 0; i < cycle.size(); ++i) {
          virt_cols.push_back({centre, cycle[i], cycle[(i + 1) % cycle.size()]});
          virt_face.push_back(face);
        }
      }
    }
    // The wall faces in column space, the face across every edge and the columns' faces.
    std::map<std::pair<Index, Index>, std::size_t> col_edge;  // directed column edge -> face
    auto add_face = [&](Index face, const std::vector<Index>& list, bool virt) {
      WallFace w;
      w.face = face;
      w.cols = list;
      w.virt = virt;
      for (std::size_t i = 0; i < list.size(); ++i) {
        col_edge[{list[i], list[(i + 1) % list.size()]}] = walls_.size();
        cols_[u(list[i])].faces.push_back(sz(walls_.size()));
      }
      walls_.push_back(w);
    };
    for (std::size_t i = 0; i < real_faces.size(); ++i) {
      add_face(real_faces[i], face_cols[i], false);
    }
    for (std::size_t i = 0; i < virt_cols.size(); ++i) {
      add_face(virt_face[i], virt_cols[i], true);
      walls_.back().cap = i >= n_fan_quads;
    }
    for (WallFace& w : walls_) {
      const std::size_t m = w.cols.size();
      for (std::size_t i = 0; i < m; ++i) {
        const Index x = w.cols[i];
        const Index y = w.cols[(i + 1) % m];
        Index other = -1;
        if (col_edge.count({y, x}) == 0) {
          const Index a = cols_[u(x)].base;
          const Index b = cols_[u(y)].base;
          if (a != b) {
            for (const Index g : edges[{std::min(a, b), std::max(a, b)}]) {
              if (!in_.wall[u(g)]) {
                other = g;
              }
            }
          } else if (!nonwall[u(a)].empty()) {
            other = nonwall[u(a)].front();
          }
        }
        w.across.push_back(other);
      }
    }
    // Neighbours and wall spacing (over neighbours on other nodes).
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
        if (cols_[u(a)].base == cols_[u(b)].base) {
          continue;
        }
        const double len = distance(points_[u(cols_[u(a)].base)], points_[u(cols_[u(b)].base)]);
        length_sum[u(a)] += len;
        length_sum[u(b)] += len;
        length_count[u(a)] += 1;
        length_count[u(b)] += 1;
      }
    }
    // Columns without neighbours on other nodes (corner centres) take their node's spacing.
    std::map<Index, std::pair<double, int>> node_length;
    for (std::size_t c = 0; c < cols_.size(); ++c) {
      auto& [sum, count] = node_length[cols_[c].base];
      sum += length_sum[c];
      count += length_count[c];
    }
    for (std::size_t c = 0; c < cols_.size(); ++c) {
      if (length_count[c] == 0) {
        const auto& [sum, count] = node_length[cols_[c].base];
        length_sum[c] = sum;
        length_count[c] = count;
      }
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
    fans_ = 0;
    node_cols_.assign(points_.size(), {});
    for (std::size_t c = 0; c < cols_.size(); ++c) {
      fans_ += static_cast<int>(cols_[c].fan);
      node_cols_[u(cols_[c].base)].push_back(sz(c));
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

  // Triangles of a wall face's current front (tags: current point ids; polygons of zero area at
  // fans are skipped).
  void front_tris(std::size_t w, std::vector<Tri>& out) const {
    const WallFace& f = walls_[w];
    for (const auto& poly : polygons(f)) {
      std::vector<Index> ids;
      ids.reserve(poly.size());
      for (const std::size_t i : poly) {
        ids.push_back(cols_[u(f.cols[i])].cur);
      }
      ids = distinct_cycle(ids);
      for (std::size_t k = 1; k + 1 < ids.size(); ++k) {
        Tri t;
        const std::array<Index, 3> tri = {ids[0], ids[k], ids[k + 1]};
        for (std::size_t i = 0; i < 3; ++i) {
          t.p[i] = points_[u(tri[i])];
          t.tag[i] = tri[i];
        }
        t.owner = sz(w);
        out.push_back(t);
      }
    }
  }

  // A polygon without repeated nodes (consecutive repeats, as at fans, removed); empty when
  // fewer than three distinct nodes remain.
  static std::vector<Index> distinct_cycle(const std::vector<Index>& ids) {
    std::vector<Index> g;
    for (const Index v : ids) {
      if (g.empty() || g.back() != v) {
        g.push_back(v);
      }
    }
    while (g.size() > 1 && g.front() == g.back()) {
      g.pop_back();
    }
    if (g.size() < 3 || std::set<Index>(g.begin(), g.end()).size() != g.size()) {
      return {};
    }
    return g;
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
      // Fan and cap faces do not count: a ridge's sector columns follow their own side's faces
      // (with the fan's faces they would turn into the fan, layer after layer).
      Vec3 sum;
      for (const Index wi : col.faces) {
        const WallFace& w = walls_[u(wi)];
        if (w.virt) {
          continue;
        }
        const std::size_t m = w.cols.size();
        std::size_t i = 0;
        while (w.cols[i] != sz(c)) {
          ++i;
        }
        const Vec3 o = cur_pos(sz(c));
        const Vec3 a = cur_pos(w.cols[(i + 1) % m]) - o;
        const Vec3 b = cur_pos(w.cols[(i + m - 1) % m]) - o;
        sum += w.inward * angle_between(a, b);
      }
      // Fan columns keep their direction (radial about the ridge): recomputed from their faces,
      // neighbouring fan columns can converge where the layers taper, folding the front.
      if (col.fan) {
        sum = col.start_normal;
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
        if (!col.active || col.fan) {
          continue;
        }
        Vec3 mean;
        int k = 0;
        for (const Index b : col.nbrs) {
          // Columns on the same node (a fan) are not averaged: their normals differ on purpose.
          if (cols_[u(b)].active && cols_[u(b)].base != col.base) {
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
      if (col.fan) {
        col.scale = 1.0;
      }
      // Spacing: mean distance to the neighbours on other nodes (a corner's centre column has
      // none: those of all the node's columns).
      double spacing = 0.0;
      int others = 0;
      for (const Index b : col.nbrs) {
        if (cols_[u(b)].base != col.base) {
          spacing += distance(cur_pos(sz(c)), cur_pos(b));
          ++others;
        }
      }
      if (others == 0) {
        for (const Index d : node_cols_[u(col.base)]) {
          for (const Index b : cols_[u(d)].nbrs) {
            if (cols_[u(b)].base != col.base) {
              spacing += distance(cur_pos(d), cur_pos(b));
              ++others;
            }
          }
        }
      }
      spacing /= static_cast<double>(std::max(others, 1));
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

  // Smallest n . inward over the column's real faces (fan and cap faces do not count).
  [[nodiscard]] double visibility(const Column& col, const Vec3& n) const {
    double worst = 1.0;
    for (const Index wi : col.faces) {
      const WallFace& w = walls_[u(wi)];
      if (!w.virt) {
        worst = std::min(worst, dot(n, w.inward));
      }
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
    // Polygons with columns sharing a node (fans and caps at ridges): cells from their faces.
    auto repeats = [&](const std::vector<std::size_t>& poly) {
      for (std::size_t i = 0; i < poly.size(); ++i) {
        for (std::size_t j = i + 1; j < poly.size(); ++j) {
          if (bottom[poly[i]] == bottom[poly[j]]) {
            return true;
          }
        }
      }
      return false;
    };
    const auto polys = polygons(f);
    bool degenerate = false;
    for (const auto& poly : polys) {
      degenerate = degenerate || repeats(poly);
    }
    if (degenerate) {
      for (const auto& poly : polys) {
        if (repeats(poly)) {
          collapsed_cell(poly, bottom, top, top_id, next_extra, extra, out);
        } else {
          tri_cell(poly[0], poly[1], poly[2]);
        }
      }
      return;
    }
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

  // Cells on corner caps are only required to be valid: the cap's triangles join columns of very
  // different directions, so their prisms cannot meet the scaled-Jacobian threshold of the
  // layers.
  void relax_caps(std::size_t w, std::vector<Element>& cells) const {
    if (walls_[w].cap) {
      for (Element& e : cells) {
        e.full = false;
      }
    }
  }

  // The cell between a polygon of a wall face whose bottom repeats nodes (columns sharing a
  // node) and its top: its faces (bottom, lateral, top) without repeated nodes, faces of zero
  // area dropped and coincident opposite faces cancelled. A tetrahedron, pyramid or prism is
  // recognised as such (a fan quad's first layer is a prism lying on its side, a cap
  // triangle's a tetrahedron); anything else is split from its centroid.
  void collapsed_cell(const std::vector<std::size_t>& poly, const std::vector<Index>& bottom,
                      const std::vector<Index>& top, const std::vector<Index>& top_id,
                      Index* next_extra, std::vector<Vec3>& extra,
                      std::vector<Element>& out) const {
    std::vector<std::vector<Index>> faces;
    auto add = [&](const std::vector<Index>& raw) {
      std::vector<Index> g;
      for (const Index v : raw) {
        if (g.empty() || g.back() != v) {
          g.push_back(v);
        }
      }
      while (g.size() > 1 && g.front() == g.back()) {
        g.pop_back();
      }
      const std::set<Index> distinct(g.begin(), g.end());
      if (distinct.size() >= 3 && distinct.size() == g.size()) {
        faces.push_back(g);
      }
    };
    const std::size_t m = poly.size();
    std::vector<Index> b(m);
    std::vector<Index> t(m);
    for (std::size_t i = 0; i < m; ++i) {
      b[i] = bottom[poly[i]];
      t[i] = top[poly[i]];
    }
    add(b);
    for (std::size_t i = 0; i < m; ++i) {
      const std::size_t j = (i + 1) % m;
      add({b[j], b[i], t[i], t[j]});
    }
    add(std::vector<Index>(t.rbegin(), t.rend()));
    // Cancel pairs of coincident faces of opposite orientation (zero-volume parts).
    auto opposite = [](const std::vector<Index>& x, const std::vector<Index>& y) {
      if (x.size() != y.size()) {
        return false;
      }
      const std::size_t n = x.size();
      for (std::size_t s = 0; s < n; ++s) {
        bool same = true;
        for (std::size_t k = 0; k < n && same; ++k) {
          same = x[k] == y[(s + n - k) % n];
        }
        if (same) {
          return true;
        }
      }
      return false;
    };
    std::vector<bool> gone(faces.size(), false);
    for (std::size_t i = 0; i < faces.size(); ++i) {
      for (std::size_t j = i + 1; j < faces.size() && !gone[i]; ++j) {
        if (!gone[j] && opposite(faces[i], faces[j])) {
          gone[i] = true;
          gone[j] = true;
        }
      }
    }
    std::vector<std::vector<Index>> kept;
    std::set<Index> vertices;
    std::vector<std::size_t> tris;
    std::vector<std::size_t> quads;
    for (std::size_t i = 0; i < faces.size(); ++i) {
      if (gone[i]) {
        continue;
      }
      if (faces[i].size() == 3) {
        tris.push_back(kept.size());
      } else {
        quads.push_back(kept.size());
      }
      vertices.insert(faces[i].begin(), faces[i].end());
      kept.push_back(faces[i]);
    }
    if (kept.empty()) {
      return;
    }
    auto other_than = [&](const std::vector<Index>& face) {
      for (const Index v : vertices) {
        if (std::find(face.begin(), face.end(), v) == face.end()) {
          return v;
        }
      }
      return Index{-1};
    };
    Element e;
    if (vertices.size() == 4 && tris.size() == 4 && quads.empty()) {
      const auto& f = kept[tris[0]];
      e.kind = Kind::tet;
      e.v = {f[0], f[2], f[1], other_than(f), 0, 0, 0, 0};
      out.push_back(e);
      return;
    }
    if (vertices.size() == 5 && tris.size() == 4 && quads.size() == 1) {
      const auto& f = kept[quads[0]];
      e.kind = Kind::pyramid;
      e.v = {f[0], f[3], f[2], f[1], other_than(f), 0, 0, 0};
      out.push_back(e);
      return;
    }
    if (vertices.size() == 6 && tris.size() == 2 && quads.size() == 3) {
      const auto& t1 = kept[tris[0]];
      const auto& t2 = kept[tris[1]];
      // Partner of x in the other triangle: the node next to x in a quad.
      auto partner = [&](Index x) {
        for (const std::size_t q : quads) {
          const auto& f = kept[q];
          for (std::size_t k = 0; k < 4; ++k) {
            if (f[k] != x) {
              continue;
            }
            for (const Index y : {f[(k + 1) % 4], f[(k + 3) % 4]}) {
              if (std::find(t2.begin(), t2.end(), y) != t2.end()) {
                return y;
              }
            }
          }
        }
        return Index{-1};
      };
      e.kind = Kind::prism;
      e.v = {t1[0], t1[2], t1[1], partner(t1[0]), partner(t1[2]), partner(t1[1]), 0, 0};
      e.full = true;
      out.push_back(e);
      return;
    }
    Vec3 centre;
    for (const Index id : vertices) {
      centre += position(id, top_id, extra);
    }
    centre = centre / static_cast<double>(vertices.size());
    const Index cid = (*next_extra)++;
    extra.push_back(centre);
    for (const auto& face : kept) {
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
        relax_caps(w, cells);
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
          // Tags: the current node of a stopped column, -(column + 1) for a new node.
          std::vector<Index> tags;
          std::map<Index, Vec3> where;
          for (const std::size_t i : poly) {
            const Index c = f.cols[i];
            const Column& col = cols_[u(c)];
            Index tag = col.cur;
            if (col.active) {
              tag = -(c + 1);
            }
            tags.push_back(tag);
            where[tag] = col.next;
          }
          tags = distinct_cycle(tags);
          for (std::size_t k = 1; k + 1 < tags.size(); ++k) {
            Tri t;
            const std::array<Index, 3> l = {tags[0], tags[k], tags[k + 1]};
            for (std::size_t i = 0; i < 3; ++i) {
              t.p[i] = where[l[i]];
              t.tag[i] = l[i];
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
        std::vector<Index> side = {b.cur, a.cur, a.cur, b.cur};
        if (a.active) {
          side[2] = top_id[u(f.cols[i])];
        }
        if (b.active) {
          side[3] = top_id[u(f.cols[(i + 1) % m])];
        }
        side = distinct_cycle(side);
        if (side.empty()) {
          continue;
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
        ids = distinct_cycle(ids);
        if (ids.empty()) {
          continue;
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
    // Nodes with several columns (fans) report their column with the most layers.
    std::vector<bool> seen(in_.points.size(), false);
    for (const Column& c : cols_) {
      const std::size_t p = u(c.base);
      if (seen[p] && c.layers <= r.layers[p]) {
        continue;
      }
      seen[p] = true;
      r.reason[p] = c.reason;
      if (c.active) {
        r.reason[p] = "layers";
      }
      r.top[p] = c.cur;
      r.height[p] = c.height;
      r.layers[p] = c.layers;
    }
    r.max_layers = max_layers_;
    r.fan_columns = fans_;
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
  int fans_ = 0;                               // fan columns
  std::vector<std::vector<Index>> node_cols_;  // per input point: its columns
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
