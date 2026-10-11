#include "opengrid/core/geom/polyline.hpp"

#include <stdexcept>
#include <string>
#include <utility>

namespace og {

std::vector<double> chord_length_knots(const std::vector<Vec3>& points, const char* what) {
  if (points.size() < 2) {
    throw std::invalid_argument(std::string(what) + ": at least 2 points are required");
  }
  std::vector<double> knots(points.size(), 0.0);
  for (std::size_t i = 1; i < points.size(); ++i) {
    const double chord = distance(points[i], points[i - 1]);
    if (chord == 0.0) {
      throw std::invalid_argument(std::string(what) + ": points[" + std::to_string(i) +
                                  "] coincides with points[" + std::to_string(i - 1) + "]");
    }
    knots[i] = knots[i - 1] + chord;
  }
  const double total = knots.back();
  for (double& k : knots) {
    k /= total;
  }
  knots.back() = 1.0;
  return knots;
}

Polyline::Polyline(std::vector<Vec3> points)
    : points_(std::move(points)), knots_(chord_length_knots(points_, "polyline")) {}

Vec3 Polyline::eval(double t) const {
  const std::size_t i = find_segment(knots_, t);
  const double u = (t - knots_[i]) / (knots_[i + 1] - knots_[i]);
  return points_[i] + u * (points_[i + 1] - points_[i]);
}

Vec3 Polyline::d1(double t) const {
  const std::size_t i = find_segment(knots_, t);
  return (points_[i + 1] - points_[i]) / (knots_[i + 1] - knots_[i]);
}

Vec3 Polyline::d2(double /*t*/) const {
  return {};
}

std::vector<double> Polyline::breakpoints() const {
  return knots_;
}

}  // namespace og
