#include "opengrid/core/quality/volume_quality.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <numbers>
#include <set>
#include <stdexcept>
#include <utility>

#include "opengrid/core/geom/vec3.hpp"
#include "opengrid/core/util/parallel.hpp"

namespace og {

namespace {

constexpr double kTiny = 1e-300;

struct Type {
  std::size_t size = 4;                     // nodes
  std::vector<std::vector<int>> faces;      // outward
  std::vector<std::array<int, 4>> corners;  // node, three neighbours
  double factor = 1.0;                      // ideal corner -> 1
  double ideal_lo = 90.0;                   // ideal smallest / largest dihedral angle
  double ideal_hi = 90.0;
  std::vector<std::pair<int, int>> edge_faces;  // the two faces at every edge
  std::vector<std::pair<int, int>> edges;       // unique edges
};

Type make(std::size_t size, std::vector<std::vector<int>> faces,
          std::vector<std::array<int, 4>> corners, double factor, double lo, double hi) {
  Type t;
  t.size = size;
  t.faces = std::move(faces);
  t.corners = std::move(corners);
  t.factor = factor;
  t.ideal_lo = lo;
  t.ideal_hi = hi;
  std::map<std::pair<int, int>, int> owner;
  for (std::size_t f = 0; f < t.faces.size(); ++f) {
    const auto& face = t.faces[f];
    for (std::size_t k = 0; k < face.size(); ++k) {
      owner[{face[k], face[(k + 1) % face.size()]}] = static_cast<int>(f);
    }
  }
  std::set<std::pair<int, int>> pairs;
  std::set<std::pair<int, int>> edges;
  for (const auto& [edge, f] : owner) {
    const int g = owner.at({edge.second, edge.first});
    pairs.insert({std::min(f, g), std::max(f, g)});
    edges.insert({std::min(edge.first, edge.second), std::max(edge.first, edge.second)});
  }
  t.edge_faces.assign(pairs.begin(), pairs.end());
  t.edges.assign(edges.begin(), edges.end());
  return t;
}

const Type& type_of(int n) {
  const double tet_angle = std::acos(1.0 / 3.0) * 180.0 / std::numbers::pi;
  static const Type tet = make(4, {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {0, 3, 2}},
                               {{0, 1, 2, 3}, {1, 2, 0, 3}, {2, 0, 1, 3}, {3, 0, 2, 1}},
                               std::numbers::sqrt2, tet_angle, tet_angle);
  static const Type pyramid =
      make(5, {{0, 3, 2, 1}, {0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}},
           {{0, 1, 3, 4}, {1, 2, 0, 4}, {2, 3, 1, 4}, {3, 0, 2, 4}}, std::numbers::sqrt2,
           std::atan(std::numbers::sqrt2) * 180.0 / std::numbers::pi, 180.0 - tet_angle);
  static const Type prism =
      make(6, {{0, 2, 1}, {3, 4, 5}, {0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}},
           {{0, 1, 2, 3}, {1, 2, 0, 4}, {2, 0, 1, 5}, {3, 5, 4, 0}, {4, 3, 5, 1}, {5, 4, 3, 2}},
           2.0 / std::numbers::sqrt3, 60.0, 90.0);
  static const Type hex =
      make(8, {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}},
           {{0, 1, 3, 4},
            {1, 2, 0, 5},
            {2, 3, 1, 6},
            {3, 0, 2, 7},
            {4, 7, 5, 0},
            {5, 4, 6, 1},
            {6, 5, 7, 2},
            {7, 6, 4, 3}},
           1.0, 90.0, 90.0);
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
  throw std::invalid_argument("volume_quality: a cell has 4, 5, 6 or 8 nodes");
}

double degrees(double radians) {
  return radians * 180.0 / std::numbers::pi;
}

// Area vector of a triangle or a bilinear quad.
Vec3 area_vector(const std::array<Vec3, 4>& p, std::size_t n) {
  if (n == 3) {
    return 0.5 * cross(p[1] - p[0], p[2] - p[0]);
  }
  return 0.5 * cross(p[2] - p[0], p[3] - p[1]);
}

double sign(double value) {
  if (value > 0.0) {
    return 1.0;
  }
  if (value < 0.0) {
    return -1.0;
  }
  return 0.0;
}

}  // namespace

VolumeQuality volume_quality(const double* points, std::size_t n_points,
                             const std::vector<CellArray>& cells) {
  // Global cell numbering and the first face of every cell (faces are numbered cell by cell).
  std::size_t n_cells = 0;
  std::size_t n_faces = 0;
  for (const CellArray& array : cells) {
    n_cells += array.count;
    n_faces += array.count * type_of(array.size).faces.size();
  }
  std::vector<const std::int64_t*> node_of(n_cells);
  std::vector<const Type*> type(n_cells);
  std::vector<std::size_t> first_face(n_cells + 1, 0);
  {
    std::size_t c = 0;
    for (const CellArray& array : cells) {
      const Type& t = type_of(array.size);
      for (std::size_t k = 0; k < array.count; ++k, ++c) {
        node_of[c] = array.nodes + k * static_cast<std::size_t>(array.size);
        type[c] = &t;
        first_face[c + 1] = first_face[c] + t.faces.size();
      }
    }
  }
  auto point = [&](std::int64_t i) {
    const double* p = points + 3 * static_cast<std::size_t>(i);
    return Vec3{p[0], p[1], p[2]};
  };
  auto face_points = [&](std::size_t c, const std::vector<int>& face) {
    std::array<Vec3, 4> p{};
    for (std::size_t k = 0; k < face.size(); ++k) {
      p[k] = point(node_of[c][face[k]]);
    }
    return p;
  };

  VolumeQuality q;
  for (auto* v : {&q.scaled_jacobian, &q.aspect_ratio, &q.min_angle, &q.max_angle, &q.skewness,
                  &q.volume, &q.non_orthogonality}) {
    v->assign(n_cells, 0.0);
  }
  q.volume_ratio.assign(n_cells, 1.0);
  std::vector<Vec3> centroid(n_cells);

  // Per-cell metrics.
  parallel_for(0, static_cast<std::int64_t>(n_cells), [&](std::int64_t index) {
    const auto c = static_cast<std::size_t>(index);
    const Type& t = *type[c];
    const std::int64_t* nodes = node_of[c];
    const std::size_t size = t.size;
    std::array<Vec3, 8> p{};
    Vec3 centre;
    for (std::size_t k = 0; k < size; ++k) {
      p[k] = point(nodes[k]);
      centre += p[k];
    }
    centroid[c] = centre / static_cast<double>(size);
    // Volume and face area vectors.
    std::array<Vec3, 6> normal{};
    double volume = 0.0;
    for (std::size_t f = 0; f < t.faces.size(); ++f) {
      const auto& face = t.faces[f];
      std::array<Vec3, 4> fp{};
      Vec3 mean;
      for (std::size_t k = 0; k < face.size(); ++k) {
        fp[k] = p[static_cast<std::size_t>(face[k])];
        mean += fp[k];
      }
      mean = mean / static_cast<double>(face.size());
      normal[f] = area_vector(fp, face.size());
      volume += dot(mean, normal[f]);
    }
    volume /= 3.0;
    q.volume[c] = volume;
    // Scaled Jacobian.
    double worst = 1.0;
    double worst_abs = 1.0;
    for (const auto& corner : t.corners) {
      const Vec3 o = p[static_cast<std::size_t>(corner[0])];
      const Vec3 e1 = p[static_cast<std::size_t>(corner[1])] - o;
      const Vec3 e2 = p[static_cast<std::size_t>(corner[2])] - o;
      const Vec3 e3 = p[static_cast<std::size_t>(corner[3])] - o;
      const double lengths = norm(e1) * norm(e2) * norm(e3);
      const double value =
          std::min(dot(e1, cross(e2, e3)) / std::max(lengths, kTiny) * t.factor, 1.0);
      worst = std::min(worst, value);
      worst_abs = std::min(worst_abs, std::fabs(value));
    }
    if (size == 4) {
      q.scaled_jacobian[c] = worst_abs * sign(volume);
    } else {
      q.scaled_jacobian[c] = worst;
    }
    // Aspect ratio.
    if (size == 4) {
      constexpr std::array<std::array<int, 3>, 4> kTetFaces = {
          {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}}};
      double area_sum = 0.0;
      for (const auto& f : kTetFaces) {
        area_sum +=
            0.5 *
            norm(cross(p[static_cast<std::size_t>(f[1])] - p[static_cast<std::size_t>(f[0])],
                       p[static_cast<std::size_t>(f[2])] - p[static_cast<std::size_t>(f[0])]));
      }
      const double inradius = 3.0 * std::fabs(volume) / std::max(area_sum, kTiny);
      const Vec3 a = p[1] - p[0];
      const Vec3 b = p[2] - p[0];
      const Vec3 d = p[3] - p[0];
      const Vec3 num = dot(d, d) * cross(a, b) + dot(a, a) * cross(b, d) + dot(b, b) * cross(d, a);
      const double circum = norm(num) / std::max(12.0 * std::fabs(volume), kTiny);
      q.aspect_ratio[c] = circum / std::max(3.0 * inradius, kTiny);
    } else {
      double longest = 0.0;
      double shortest = 0.0;
      bool first = true;
      for (const auto& [a, b] : t.edges) {
        const double length =
            distance(p[static_cast<std::size_t>(a)], p[static_cast<std::size_t>(b)]);
        if (first) {
          longest = length;
          shortest = length;
          first = false;
        }
        longest = std::max(longest, length);
        shortest = std::min(shortest, length);
      }
      q.aspect_ratio[c] = longest / std::max(shortest, kTiny);
    }
    // Dihedral angles and skewness.
    double lo = 0.0;
    double hi = 0.0;
    bool first = true;
    for (const auto& [fa, fb] : t.edge_faces) {
      const Vec3& na = normal[static_cast<std::size_t>(fa)];
      const Vec3& nb = normal[static_cast<std::size_t>(fb)];
      const double cosine = dot(na, nb) / std::max(norm(na) * norm(nb), kTiny);
      const double dihedral = 180.0 - degrees(std::acos(std::clamp(cosine, -1.0, 1.0)));
      if (first) {
        lo = dihedral;
        hi = dihedral;
        first = false;
      }
      lo = std::min(lo, dihedral);
      hi = std::max(hi, dihedral);
    }
    q.min_angle[c] = lo;
    q.max_angle[c] = hi;
    q.skewness[c] =
        std::max({(hi - t.ideal_hi) / (180.0 - t.ideal_hi), (t.ideal_lo - lo) / t.ideal_lo, 0.0});
  });

  // Face matching: bucket faces by their smallest node (counting sort), then compare sorted node
  // lists within a bucket.
  std::vector<std::size_t> face_cell(n_faces);
  std::vector<std::array<std::int64_t, 4>> key(n_faces);
  std::vector<std::uint8_t> key_size(n_faces);
  parallel_for(0, static_cast<std::int64_t>(n_cells), [&](std::int64_t index) {
    const auto c = static_cast<std::size_t>(index);
    const Type& t = *type[c];
    for (std::size_t f = 0; f < t.faces.size(); ++f) {
      const std::size_t id = first_face[c] + f;
      const auto& face = t.faces[f];
      std::array<std::int64_t, 4> k{-1, -1, -1, -1};
      for (std::size_t m = 0; m < face.size(); ++m) {
        k[m] = node_of[c][face[m]];
      }
      std::sort(k.begin(), k.begin() + static_cast<std::ptrdiff_t>(face.size()));
      key[id] = k;
      key_size[id] = static_cast<std::uint8_t>(face.size());
      face_cell[id] = c;
    }
  });
  std::vector<std::size_t> bucket_start(n_points + 1, 0);
  for (std::size_t id = 0; id < n_faces; ++id) {
    bucket_start[static_cast<std::size_t>(key[id][0]) + 1] += 1;
  }
  for (std::size_t i = 0; i < n_points; ++i) {
    bucket_start[i + 1] += bucket_start[i];
  }
  std::vector<std::size_t> bucket(n_faces);
  {
    std::vector<std::size_t> fill(bucket_start.begin(), bucket_start.end() - 1);
    for (std::size_t id = 0; id < n_faces; ++id) {
      bucket[fill[static_cast<std::size_t>(key[id][0])]++] = id;
    }
  }
  constexpr auto kNone = static_cast<std::size_t>(-1);
  std::vector<std::size_t> other(n_faces, kNone);
  parallel_for(0, static_cast<std::int64_t>(n_points), [&](std::int64_t index) {
    const auto i = static_cast<std::size_t>(index);
    for (std::size_t a = bucket_start[i]; a < bucket_start[i + 1]; ++a) {
      const std::size_t fa = bucket[a];
      if (other[fa] != kNone) {
        continue;
      }
      for (std::size_t b = a + 1; b < bucket_start[i + 1]; ++b) {
        const std::size_t fb = bucket[b];
        if (other[fb] == kNone && key_size[fa] == key_size[fb] && key[fa] == key[fb]) {
          other[fa] = face_cell[fb];
          other[fb] = face_cell[fa];
          break;
        }
      }
    }
  });

  // Face-based metrics, reduced per cell (its faces are contiguous).
  parallel_for(0, static_cast<std::int64_t>(n_cells), [&](std::int64_t index) {
    const auto c = static_cast<std::size_t>(index);
    const Type& t = *type[c];
    double non_orth = 0.0;
    double ratio = 1.0;
    for (std::size_t f = 0; f < t.faces.size(); ++f) {
      const std::size_t id = first_face[c] + f;
      const auto& face = t.faces[f];
      const std::array<Vec3, 4> fp = face_points(c, face);
      const Vec3 area = area_vector(fp, face.size());
      Vec3 link;
      if (other[id] != kNone) {
        link = centroid[other[id]] - centroid[c];
        const double va = std::fabs(q.volume[c]);
        const double vb = std::fabs(q.volume[other[id]]);
        const double r = std::max(va / vb, vb / va);
        if (std::isfinite(r)) {
          ratio = std::max(ratio, r);
        }
      } else {
        Vec3 mean;
        for (std::size_t k = 0; k < face.size(); ++k) {
          mean += fp[k];
        }
        link = mean / static_cast<double>(face.size()) - centroid[c];
      }
      const double cosine = std::fabs(dot(area, link)) / std::max(norm(area) * norm(link), kTiny);
      non_orth = std::max(non_orth, degrees(std::acos(std::clamp(cosine, 0.0, 1.0))));
    }
    q.non_orthogonality[c] = non_orth;
    q.volume_ratio[c] = ratio;
  });
  return q;
}

}  // namespace og
