// Minimal 3-component double-precision vector used throughout the geometry kernel.
#pragma once

#include <array>
#include <cmath>

namespace og {

struct Vec3 {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;

  constexpr Vec3() = default;
  constexpr Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
  explicit constexpr Vec3(const std::array<double, 3>& a) : x(a[0]), y(a[1]), z(a[2]) {}

  [[nodiscard]] constexpr std::array<double, 3> to_array() const {
    return {x, y, z};
  }

  constexpr Vec3& operator+=(const Vec3& o) {
    x += o.x;
    y += o.y;
    z += o.z;
    return *this;
  }
  constexpr Vec3& operator-=(const Vec3& o) {
    x -= o.x;
    y -= o.y;
    z -= o.z;
    return *this;
  }
  constexpr Vec3& operator*=(double s) {
    x *= s;
    y *= s;
    z *= s;
    return *this;
  }
};

constexpr Vec3 operator+(Vec3 a, const Vec3& b) {
  a += b;
  return a;
}
constexpr Vec3 operator-(Vec3 a, const Vec3& b) {
  a -= b;
  return a;
}
constexpr Vec3 operator-(const Vec3& a) {
  return {-a.x, -a.y, -a.z};
}
constexpr Vec3 operator*(Vec3 a, double s) {
  a *= s;
  return a;
}
constexpr Vec3 operator*(double s, Vec3 a) {
  a *= s;
  return a;
}
constexpr Vec3 operator/(const Vec3& a, double s) {
  return {a.x / s, a.y / s, a.z / s};
}

constexpr double dot(const Vec3& a, const Vec3& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}
constexpr Vec3 cross(const Vec3& a, const Vec3& b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm(const Vec3& a) {
  return std::sqrt(dot(a, a));
}
inline double distance(const Vec3& a, const Vec3& b) {
  return norm(a - b);
}

}  // namespace og
