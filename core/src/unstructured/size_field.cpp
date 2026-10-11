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

}  // namespace

SizeField::SizeField(std::vector<std::array<double, 2>> points, std::vector<double> sizes,
                     const SizeFieldOptions& options)
    : points_(std::move(points)),
      sizes_(std::move(sizes)),
      options_(options),
      gradient_(options.growth_rate - 1.0),
      min_size_(std::numeric_limits<double>::infinity()) {
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
    min_size_ = std::min(min_size_, s);
  }
  for (const SizeSource& source : options.sources) {
    if (!(source.size > 0.0) || source.radius < 0.0) {
      throw std::invalid_argument("size field: sources need a positive size and radius >= 0");
    }
  }
  if (options.max_size && !(*options.max_size > 0.0)) {
    throw std::invalid_argument("size field: max size must be positive");
  }
  double xmin = points_[0][0];
  double xmax = xmin;
  double ymin = points_[0][1];
  double ymax = ymin;
  double mean = 0.0;
  for (std::size_t i = 0; i < points_.size(); ++i) {
    xmin = std::min(xmin, points_[i][0]);
    xmax = std::max(xmax, points_[i][0]);
    ymin = std::min(ymin, points_[i][1]);
    ymax = std::max(ymax, points_[i][1]);
    mean += sizes_[i];
  }
  mean /= static_cast<double>(points_.size());
  // About a few points per bucket.
  cell_ = std::max(4.0 * mean, 1e-300);
  origin_x_ = xmin;
  origin_y_ = ymin;
  nx_ = std::max(1L, static_cast<long>(std::floor((xmax - xmin) / cell_)) + 1);
  ny_ = std::max(1L, static_cast<long>(std::floor((ymax - ymin) / cell_)) + 1);
  if (nx_ * ny_ > 4'000'000) {
    const double factor = std::sqrt(static_cast<double>(nx_ * ny_) / 4'000'000.0);
    cell_ *= factor;
    nx_ = std::max(1L, static_cast<long>(std::floor((xmax - xmin) / cell_)) + 1);
    ny_ = std::max(1L, static_cast<long>(std::floor((ymax - ymin) / cell_)) + 1);
  }
  buckets_.assign(static_cast<std::size_t>(nx_ * ny_), {});
  for (std::size_t i = 0; i < points_.size(); ++i) {
    const long cx = std::min(nx_ - 1, static_cast<long>((points_[i][0] - origin_x_) / cell_));
    const long cy = std::min(ny_ - 1, static_cast<long>((points_[i][1] - origin_y_) / cell_));
    buckets_[static_cast<std::size_t>(cy * nx_ + cx)].push_back(i);
  }
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
  const long cx = std::clamp(static_cast<long>(std::floor((x - origin_x_) / cell_)), -1L, nx_);
  const long cy = std::clamp(static_cast<long>(std::floor((y - origin_y_) / cell_)), -1L, ny_);
  const long max_ring = std::max(nx_, ny_) + 2;
  for (long ring = 0; ring <= max_ring; ++ring) {
    // Points in ring r are at least (r - 1) cells away.
    const double bound =
        min_size_ + gradient_ * static_cast<double>(std::max(0L, ring - 1)) * cell_;
    if (bound >= best) {
      break;
    }
    for (long j = cy - ring; j <= cy + ring; ++j) {
      for (long i = cx - ring; i <= cx + ring; ++i) {
        if (std::max(std::labs(i - cx), std::labs(j - cy)) != ring) {
          continue;
        }
        if (i < 0 || j < 0 || i >= nx_ || j >= ny_) {
          continue;
        }
        for (const std::size_t k : buckets_[static_cast<std::size_t>(j * nx_ + i)]) {
          const double d = std::hypot(x - points_[k][0], y - points_[k][1]);
          best = std::min(best, sizes_[k] + gradient_ * d);
        }
      }
    }
  }
  return best;
}

}  // namespace og
