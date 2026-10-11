#include "opengrid/core/quality/quality_tri.hpp"

#include "opengrid/core/quality/quality_cells.hpp"

namespace og {

Quality2D quality_triangles(const std::vector<std::array<double, 2>>& points,
                            const std::vector<std::array<std::int64_t, 3>>& triangles) {
  return quality_cells(points, triangles, {});
}

}  // namespace og
