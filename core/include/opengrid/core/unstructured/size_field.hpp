// Mesh size field for 2D unstructured meshing.
//
// h(x) = min( max_size,
//             min_j ( h_j + g |x - x_j| )                       over boundary nodes j,
//             min_s ( size_s + g max(0, dist(x, source_s) - radius_s) ) over sources )
// with g = growth_rate - 1: the size grows linearly away from the boundary (whose local size h_j
// is the mean length of the segments at node j) and from point/line sources. The boundary term
// is evaluated by branch and bound on a k-d tree whose nodes store their bounding box and the
// smallest h_j below them: a subtree is skipped when min h + g dist(x, box) cannot improve the
// current minimum, so the result equals the brute-force minimum exactly (the minimum does not
// depend on the visiting order) at O(log n) typical cost per query.
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

  struct Node {
    std::array<double, 4> box{};  // xmin, ymin, xmax, ymax
    double min_size = 0.0;
    std::size_t begin = 0;  // range in points_ (reordered) for leaves
    std::size_t end = 0;
    std::size_t left = 0;  // child node indices; 0 for a leaf (the root is never a child)
    std::size_t right = 0;
  };
  std::size_t build(std::size_t begin, std::size_t end, int depth);
  std::vector<Node> nodes_;
};

}  // namespace og
