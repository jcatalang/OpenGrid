// Sharp corners of a closed polygonal loop, shared by the all-quad coarse boundary
// (domain_mesher.cpp) and the boundary-layer front (inflation.cpp) so both agree on which nodes
// must stay corner nodes of the coarse boundary.
#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <vector>

namespace og {

// Absolute turning angle (radians) of the loop at each of its nodes.
inline std::vector<double> loop_turning(const std::vector<std::array<double, 2>>& points,
                                        const std::vector<std::int64_t>& loop) {
  const std::size_t m = loop.size();
  std::vector<double> turn(m, 0.0);
  for (std::size_t k = 0; k < m; ++k) {
    const auto& a = points[static_cast<std::size_t>(loop[(k + m - 1) % m])];
    const auto& b = points[static_cast<std::size_t>(loop[k])];
    const auto& c = points[static_cast<std::size_t>(loop[(k + 1) % m])];
    const double ux = b[0] - a[0];
    const double uy = b[1] - a[1];
    const double vx = c[0] - b[0];
    const double vy = c[1] - b[1];
    turn[k] = std::fabs(std::atan2(ux * vy - uy * vx, ux * vx + uy * vy));
  }
  return turn;
}

// Positions of the loop's corners: nodes that turn by more than 35 degrees.
inline std::vector<std::size_t> loop_corners(const std::vector<std::array<double, 2>>& points,
                                             const std::vector<std::int64_t>& loop) {
  constexpr double kCorner = 35.0 * std::numbers::pi / 180.0;
  const std::vector<double> turn = loop_turning(points, loop);
  std::vector<std::size_t> corners;
  for (std::size_t k = 0; k < turn.size(); ++k) {
    if (turn[k] > kCorner) {
      corners.push_back(k);
    }
  }
  return corners;
}

}  // namespace og
