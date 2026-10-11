// Mesh size field for 2D unstructured meshing.
//
// h(x) = min( max_size,
//             min_j ( h_j + g |x - x_j| )                       over boundary nodes j,
//             min_s ( size_s + g max(0, dist(x, source_s) - radius_s) ) over sources )
// with g = growth_rate - 1: the size grows linearly away from the boundary (whose local size h_j
// is the mean length of the segments at node j) and from point/line sources. The boundary term
// is evaluated with a uniform bucket grid and an expanding ring search with an exact lower
// bound, so it equals the brute-force minimum.
#pragma once

#include <array>
#include <optional>
#include <vector>

namespace og {

struct SizeSource {
  std::array<double, 2> p0{};
  std::array<double, 2> p1{};  // equal to p0 for a point source
  double size = 1.0;
  double radius = 0.0;
};

struct SizeFieldOptions {
  double growth_rate = 1.2;  // > 1
  std::optional<double> max_size;
  std::vector<SizeSource> sources;
};

class SizeField {
 public:
  // Throws std::invalid_argument for growth_rate <= 1, non-positive sizes or empty input.
  SizeField(std::vector<std::array<double, 2>> points, std::vector<double> sizes,
            const SizeFieldOptions& options);

  [[nodiscard]] double operator()(double x, double y) const;

 private:
  std::vector<std::array<double, 2>> points_;
  std::vector<double> sizes_;
  SizeFieldOptions options_;
  double gradient_;
  double min_size_;
  double origin_x_ = 0.0;
  double origin_y_ = 0.0;
  double cell_ = 1.0;
  long nx_ = 1;
  long ny_ = 1;
  std::vector<std::vector<std::size_t>> buckets_;
};

}  // namespace og
