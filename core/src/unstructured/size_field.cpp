#include "opengrid/core/unstructured/size_field.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace og {
namespace {

double segment_distance(double x, double y, const std::array<double, 2>& a,
                        const std::array<double, 2>& b) {
  const double dx = b[0] - a[0];
  const double dy = b[1] - a[1];
  const double len2 = dx * dx + dy * dy;
  double t = 0.0;
  if (len2 > 0.0) {
    t = std::clamp(((x - a[0]) * dx + (y - a[1]) * dy) / len2, 0.0, 1.0);
  }
  return std::hypot(x - (a[0] + t * dx), y - (a[1] + t * dy));
}

constexpr std::size_t kLeafSize = 8;

double box_distance(double x, double y, const std::array<double, 4>& box) {
  const double dx = std::max({box[0] - x, 0.0, x - box[2]});
  const double dy = std::max({box[1] - y, 0.0, y - box[3]});
  return std::hypot(dx, dy);
}

}  // namespace

SizeField::SizeField(std::vector<std::array<double, 2>> points, std::vector<double> sizes,
                     const SizeFieldOptions& options)
    : points_(std::move(points)),
      sizes_(std::move(sizes)),
      options_(options),
      gradient_(options.growth_rate - 1.0) {
  if (points_.empty() || points_.size() != sizes_.size()) {
    throw std::invalid_argument("size field: need one size per boundary point");
  }
  if (!(options.growth_rate > 1.0)) {
    throw std::invalid_argument("size field: growth rate must be greater than 1");
  }
  for (const double s : sizes_) {
    if (!(s > 0.0)) {
      throw std::invalid_argument("size field: sizes must be positive");
    }
  }
  for (const SizeSource& source : options.sources) {
    if (!(source.size > 0.0) || source.radius < 0.0) {
      throw std::invalid_argument("size field: sources need a positive size and radius >= 0");
    }
  }
  if (options.max_size && !(*options.max_size > 0.0)) {
    throw std::invalid_argument("size field: max size must be positive");
  }
  // The tree build reorders points_ and sizes_ together.
  nodes_.reserve(2 * points_.size() / kLeafSize + 2);
  build(0, points_.size(), 0);
}

double SizeField::operator()(double x, double y) const {
  double best = std::numeric_limits<double>::infinity();
  if (options_.max_size) {
    best = *options_.max_size;
  }
  for (const SizeSource& source : options_.sources) {
    const double d = segment_distance(x, y, source.p0, source.p1);
    best = std::min(best, source.size + gradient_ * std::max(0.0, d - source.radius));
  }
  // Depth-first branch and bound, nearer child first.
  std::array<std::size_t, 128> stack{};
  std::size_t top = 0;
  stack[top++] = 0;
  while (top > 0) {
    const Node& node = nodes_[stack[--top]];
    if (node.min_size + gradient_ * box_distance(x, y, node.box) >= best) {
      continue;
    }
    if (node.left == 0) {
      for (std::size_t k = node.begin; k < node.end; ++k) {
        const double d = std::hypot(x - points_[k][0], y - points_[k][1]);
        best = std::min(best, sizes_[k] + gradient_ * d);
      }
      continue;
    }
    const double dl = box_distance(x, y, nodes_[node.left].box);
    const double dr = box_distance(x, y, nodes_[node.right].box);
    if (dl <= dr) {
      stack[top++] = node.right;
      stack[top++] = node.left;
    } else {
      stack[top++] = node.left;
      stack[top++] = node.right;
    }
  }
  return best;
}

std::size_t SizeField::build(std::size_t begin, std::size_t end, int depth) {
  const std::size_t id = nodes_.size();
  nodes_.emplace_back();
  Node node;
  node.begin = begin;
  node.end = end;
  node.box = {points_[begin][0], points_[begin][1], points_[begin][0], points_[begin][1]};
  node.min_size = sizes_[begin];
  for (std::size_t k = begin; k < end; ++k) {
    node.box[0] = std::min(node.box[0], points_[k][0]);
    node.box[1] = std::min(node.box[1], points_[k][1]);
    node.box[2] = std::max(node.box[2], points_[k][0]);
    node.box[3] = std::max(node.box[3], points_[k][1]);
    node.min_size = std::min(node.min_size, sizes_[k]);
  }
  if (end - begin > kLeafSize && depth < 60) {
    // Split at the median of the longer box side (points and sizes move together).
    std::size_t axis = 0;
    if (node.box[3] - node.box[1] > node.box[2] - node.box[0]) {
      axis = 1;
    }
    const std::size_t mid = begin + (end - begin) / 2;
    std::vector<std::size_t> idx(end - begin);
    for (std::size_t k = 0; k < idx.size(); ++k) {
      idx[k] = begin + k;
    }
    std::nth_element(idx.begin(), idx.begin() + static_cast<std::ptrdiff_t>(mid - begin), idx.end(),
                     [&](std::size_t a, std::size_t b) {
                       if (points_[a][axis] != points_[b][axis]) {
                         return points_[a][axis] < points_[b][axis];
                       }
                       return a < b;
                     });
    std::vector<std::array<double, 2>> p(idx.size());
    std::vector<double> h(idx.size());
    for (std::size_t k = 0; k < idx.size(); ++k) {
      p[k] = points_[idx[k]];
      h[k] = sizes_[idx[k]];
    }
    std::copy(p.begin(), p.end(), points_.begin() + static_cast<std::ptrdiff_t>(begin));
    std::copy(h.begin(), h.end(), sizes_.begin() + static_cast<std::ptrdiff_t>(begin));
    node.left = build(begin, mid, depth + 1);
    node.right = build(mid, end, depth + 1);
  }
  nodes_[id] = node;
  return id;
}

}  // namespace og
