// nanobind module `opengrid._core`: Python view of the C++ core library.
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdlib>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "opengrid/core/dist/distribution.hpp"
#include "opengrid/core/geom/arc.hpp"
#include "opengrid/core/geom/cubic_spline.hpp"
#include "opengrid/core/geom/curve_chain.hpp"
#include "opengrid/core/geom/line.hpp"
#include "opengrid/core/geom/nurbs.hpp"
#include "opengrid/core/geom/polyline.hpp"
#include "opengrid/core/topo/connector_nodes.hpp"
#include "opengrid/core/version.hpp"

namespace nb = nanobind;

namespace {

using Array3 = std::array<double, 3>;
using PointArray = nb::ndarray<nb::numpy, double, nb::shape<-1, 3>>;
using ParamArray = nb::ndarray<const double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

std::vector<og::Vec3> to_vec3(const std::vector<Array3>& points) {
  std::vector<og::Vec3> out;
  out.reserve(points.size());
  for (const Array3& p : points) {
    out.emplace_back(p);
  }
  return out;
}

std::optional<og::Vec3> to_optional_vec3(const std::optional<Array3>& a) {
  if (a) {
    return og::Vec3(*a);
  }
  return std::nullopt;
}

// (n, 3) NumPy array owning a copy of the given points.
PointArray to_numpy(const std::vector<og::Vec3>& points) {
  auto* data = new double[points.size() * 3];
  for (std::size_t i = 0; i < points.size(); ++i) {
    data[3 * i] = points[i].x;
    data[3 * i + 1] = points[i].y;
    data[3 * i + 2] = points[i].z;
  }
  const nb::capsule owner(data, [](void* ptr) noexcept { delete[] static_cast<double*>(ptr); });
  return PointArray(data, {points.size(), 3}, owner);
}

using ScalarArray = nb::ndarray<nb::numpy, double, nb::ndim<1>>;

ScalarArray to_numpy(const std::vector<double>& values) {
  auto* data = new double[values.size()];
  std::copy(values.begin(), values.end(), data);
  const nb::capsule owner(data, [](void* ptr) noexcept { delete[] static_cast<double*>(ptr); });
  return ScalarArray(data, {values.size()}, owner);
}

og::DistributionType distribution_type(const std::string& name) {
  if (name == "uniform") {
    return og::DistributionType::Uniform;
  }
  if (name == "geometric") {
    return og::DistributionType::Geometric;
  }
  if (name == "tanh") {
    return og::DistributionType::Tanh;
  }
  if (name == "sine") {
    return og::DistributionType::Sine;
  }
  throw std::invalid_argument("unknown distribution type '" + name + "'");
}

og::DistributionSpec make_spec(const std::string& type, std::optional<double> begin_spacing,
                               std::optional<double> end_spacing,
                               std::optional<double> growth_ratio) {
  og::DistributionSpec spec;
  spec.type = distribution_type(type);
  spec.begin_spacing = begin_spacing;
  spec.end_spacing = end_spacing;
  spec.growth_ratio = growth_ratio;
  return spec;
}

// (start, end, n_points, type, begin_spacing, end_spacing, growth_ratio)
using SubTuple = std::tuple<double, double, int, std::string, std::optional<double>,
                            std::optional<double>, std::optional<double>>;
using SegmentTuple = std::tuple<std::shared_ptr<og::Curve>, double, double>;

std::vector<og::Vec3> eval_params(const og::Curve& curve, const double* ts, std::size_t n) {
  std::vector<og::Vec3> points;
  points.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    points.push_back(curve.eval(ts[i]));
  }
  return points;
}

}  // namespace

NB_MODULE(_core, m) {  // NOLINT(bugprone-signed-bitwise): inside the nanobind macro
  m.doc() = "OpenGrid C++ core";
  // Qt/VTK keep Python callbacks (and through them curve objects) alive until after the
  // interpreter has torn this module down, which nanobind reports as leaks at exit. These are
  // shutdown-order artefacts, not reference-counting errors; the API-only test suite runs
  // leak-free with the warnings enabled (OG_LEAK_WARNINGS=1 re-enables them).
  nb::set_leak_warnings(std::getenv("OG_LEAK_WARNINGS") != nullptr);
  m.def("version", &og::version, "Version string of the C++ core library.");
  m.def("add", &og::add, nb::arg("a"), nb::arg("b"), "Add two 64-bit integers (build smoke test).");

  nb::class_<og::Curve>(m, "Curve", "Parametric curve x(t), t in [0, 1].")
      .def(
          "eval", [](const og::Curve& c, double t) { return c.eval(t).to_array(); }, nb::arg("t"))
      .def(
          "d1", [](const og::Curve& c, double t) { return c.d1(t).to_array(); }, nb::arg("t"),
          "First derivative dx/dt.")
      .def(
          "d2", [](const og::Curve& c, double t) { return c.d2(t).to_array(); }, nb::arg("t"),
          "Second derivative d2x/dt2.")
      .def(
          "eval_many",
          [](const og::Curve& c, const ParamArray& ts) {
            return to_numpy(eval_params(c, ts.data(), ts.shape(0)));
          },
          nb::arg("t"), "Evaluate at an array of parameters; returns an (n, 3) array.")
      .def("breakpoints", &og::Curve::breakpoints)
      .def("length", &og::Curve::length)
      .def("length_between", &og::Curve::length_between, nb::arg("t0"), nb::arg("t1"))
      .def("t_at_length", &og::Curve::t_at_length, nb::arg("s"))
      .def(
          "project",
          [](const og::Curve& c, const Array3& p) {
            const og::CurveProjection proj = c.project(og::Vec3(p));
            return std::make_tuple(proj.t, proj.point.to_array(), proj.distance);
          },
          nb::arg("point"), "Closest point: returns (t, point, distance).")
      .def("tessellate", &og::Curve::tessellate, nb::arg("rel_tol") = 1e-4,
           "Parameters of an adaptive polyline approximation.")
      .def(
          "tessellate_points",
          [](const og::Curve& c, double rel_tol) {
            const std::vector<double> ts = c.tessellate(rel_tol);
            return to_numpy(eval_params(c, ts.data(), ts.size()));
          },
          nb::arg("rel_tol") = 1e-4, "Points of an adaptive polyline approximation, (n, 3).");

  nb::class_<og::Line, og::Curve>(m, "Line").def(
      "__init__",
      [](og::Line* self, const Array3& p0, const Array3& p1) {
        new (self) og::Line(og::Vec3(p0), og::Vec3(p1));
      },
      nb::arg("p0"), nb::arg("p1"));

  nb::class_<og::Arc, og::Curve>(m, "Arc").def(
      "__init__",
      [](og::Arc* self, const Array3& center, double radius, const Array3& normal,
         const Array3& x_axis, double start_angle, double end_angle) {
        new (self) og::Arc(og::Vec3(center), radius, og::Vec3(normal), og::Vec3(x_axis),
                           start_angle, end_angle);
      },
      nb::arg("center"), nb::arg("radius"), nb::arg("normal"), nb::arg("x_axis"),
      nb::arg("start_angle"), nb::arg("end_angle"));

  nb::class_<og::Polyline, og::Curve>(m, "Polyline")
      .def(
          "__init__",
          [](og::Polyline* self, const std::vector<Array3>& points) {
            new (self) og::Polyline(to_vec3(points));
          },
          nb::arg("points"));

  nb::enum_<og::SplineMethod>(m, "SplineMethod")
      .value("NATURAL", og::SplineMethod::Natural)
      .value("CLAMPED", og::SplineMethod::Clamped)
      .value("CATMULL_ROM", og::SplineMethod::CatmullRom);

  nb::class_<og::CubicSpline, og::Curve>(m, "CubicSpline")
      .def(
          "__init__",
          [](og::CubicSpline* self, const std::vector<Array3>& points, og::SplineMethod method,
             const std::optional<Array3>& start_tangent, const std::optional<Array3>& end_tangent) {
            new (self) og::CubicSpline(to_vec3(points), method, to_optional_vec3(start_tangent),
                                       to_optional_vec3(end_tangent));
          },
          nb::arg("points"), nb::arg("method"), nb::arg("start_tangent") = nb::none(),
          nb::arg("end_tangent") = nb::none())
      .def_prop_ro("knots", &og::CubicSpline::knots);

  m.def(
      "distribute",
      [](int n_points, const std::string& type, std::optional<double> begin_spacing,
         std::optional<double> end_spacing, std::optional<double> growth_ratio) {
        const og::Distribution d =
            og::distribute(n_points, make_spec(type, begin_spacing, end_spacing, growth_ratio));
        return std::make_tuple(to_numpy(d.s), d.warnings);
      },
      nb::arg("n_points"), nb::arg("type") = "uniform", nb::arg("begin_spacing") = nb::none(),
      nb::arg("end_spacing") = nb::none(), nb::arg("growth_ratio") = nb::none(),
      "Distribution on [0, 1] (spacings as fractions); returns (s, warnings).");
  m.def("max_stretching_ratio", &og::max_stretching_ratio, nb::arg("positions"));

  nb::class_<og::CurveChain>(m, "CurveChain", "Chain of curve pieces (curve, t0, t1).")
      .def(
          "__init__",
          [](og::CurveChain* self, const std::vector<SegmentTuple>& segments, double tolerance) {
            std::vector<og::ChainSegment> pieces;
            pieces.reserve(segments.size());
            for (const auto& [curve, t0, t1] : segments) {
              pieces.push_back({curve, t0, t1});
            }
            new (self) og::CurveChain(std::move(pieces), tolerance);
          },
          nb::arg("segments"), nb::arg("tolerance"))
      .def("length", &og::CurveChain::length)
      .def("start", [](const og::CurveChain& c) { return c.start().to_array(); })
      .def("end", [](const og::CurveChain& c) { return c.end().to_array(); })
      .def(
          "point_at_length",
          [](const og::CurveChain& c, double s) { return c.point_at_length(s).to_array(); },
          nb::arg("s"))
      .def(
          "point_at_param",
          [](const og::CurveChain& c, double u) { return c.point_at_param(u).to_array(); },
          nb::arg("u"));

  m.def(
      "connector_nodes",
      [](const og::CurveChain& chain, const std::string& mode, const std::vector<SubTuple>& subs) {
        og::SpacingMode spacing_mode = og::SpacingMode::ArcLength;
        if (mode == "parametric") {
          spacing_mode = og::SpacingMode::Parametric;
        } else if (mode != "arc_length") {
          throw std::invalid_argument("unknown spacing mode '" + mode + "'");
        }
        std::vector<og::SubConnectorSpec> specs;
        for (const auto& [start, end, n_points, type, begin, finish, ratio] : subs) {
          og::SubConnectorSpec spec;
          spec.start = start;
          spec.end = end;
          spec.n_points = n_points;
          spec.distribution = make_spec(type, begin, finish, ratio);
          specs.push_back(spec);
        }
        const og::ConnectorNodes nodes = og::connector_nodes(chain, spacing_mode, specs);
        return std::make_tuple(to_numpy(nodes.points), to_numpy(nodes.arc_positions),
                               nodes.warnings, nodes.max_stretching_ratio);
      },
      nb::arg("chain"), nb::arg("mode"), nb::arg("subs"),
      "Connector nodes: returns (points (n, 3), arc_positions (n,), warnings, max ratio).");

  nb::class_<og::Nurbs, og::Curve>(m, "Nurbs")
      .def(
          "__init__",
          [](og::Nurbs* self, int degree, const std::vector<double>& knots,
             const std::vector<Array3>& control_points, const std::vector<double>& weights) {
            new (self) og::Nurbs(degree, knots, to_vec3(control_points), weights);
          },
          nb::arg("degree"), nb::arg("knots"), nb::arg("control_points"), nb::arg("weights"));
}
