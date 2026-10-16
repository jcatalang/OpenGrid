// nanobind module `opengrid._core`: Python view of the C++ core library.
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdlib>
#include <map>
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
#include "opengrid/core/gpu/gpu.hpp"
#include "opengrid/core/highorder/elements.hpp"
#include "opengrid/core/highorder/nodes.hpp"
#include "opengrid/core/io/cgns_structured.hpp"
#include "opengrid/core/quality/quality2d.hpp"
#include "opengrid/core/quality/quality_cells.hpp"
#include "opengrid/core/quality/quality_tri.hpp"
#include "opengrid/core/structured/elliptic.hpp"
#include "opengrid/core/structured/multiblock_elliptic.hpp"
#include "opengrid/core/structured/tfi.hpp"
#include "opengrid/core/topo/connector_nodes.hpp"
#include "opengrid/core/unstructured/advancing_front.hpp"
#include "opengrid/core/unstructured/domain_mesher.hpp"
#include "opengrid/core/unstructured/triangulate.hpp"
#include "opengrid/core/util/parallel.hpp"
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
using GridIn = nb::ndarray<const double, nb::ndim<3>, nb::c_contig, nb::device::cpu>;
using GridOut = nb::ndarray<nb::numpy, double, nb::ndim<3>>;
using CellsOut = nb::ndarray<nb::numpy, double, nb::ndim<2>>;
using PointsIn = nb::ndarray<const double, nb::shape<-1, 3>, nb::c_contig, nb::device::cpu>;

std::vector<og::Vec3> points_from(const PointsIn& a) {
  std::vector<og::Vec3> out;
  out.reserve(a.shape(0));
  for (std::size_t k = 0; k < a.shape(0); ++k) {
    out.emplace_back(a(k, 0), a(k, 1), a(k, 2));
  }
  return out;
}

og::Grid2D grid_from(const GridIn& a) {
  if (a.shape(2) != 3) {
    throw std::invalid_argument("grid must have shape (nj, ni, 3)");
  }
  og::Grid2D g(static_cast<int>(a.shape(1)), static_cast<int>(a.shape(0)));
  const double* data = a.data();
  for (std::size_t k = 0; k < g.x.size(); ++k) {
    g.x[k] = {data[3 * k], data[3 * k + 1], data[3 * k + 2]};
  }
  return g;
}

GridOut grid_to(const og::Grid2D& g) {
  auto* data = new double[g.x.size() * 3];
  for (std::size_t k = 0; k < g.x.size(); ++k) {
    data[3 * k] = g.x[k].x;
    data[3 * k + 1] = g.x[k].y;
    data[3 * k + 2] = g.x[k].z;
  }
  const nb::capsule owner(data, [](void* ptr) noexcept { delete[] static_cast<double*>(ptr); });
  return GridOut(data, {static_cast<std::size_t>(g.nj), static_cast<std::size_t>(g.ni), 3}, owner);
}

CellsOut cells_to(const std::vector<double>& values, int ci, int cj) {
  auto* data = new double[values.size()];
  std::copy(values.begin(), values.end(), data);
  const nb::capsule owner(data, [](void* ptr) noexcept { delete[] static_cast<double*>(ptr); });
  return CellsOut(data, {static_cast<std::size_t>(cj), static_cast<std::size_t>(ci)}, owner);
}

// (spacing or None, decay, relaxation, max_source)
using WallTuple = std::tuple<std::optional<double>, double, double, double>;

og::WallControl wall_from(const std::optional<WallTuple>& wall) {
  og::WallControl control;
  if (wall) {
    control.enabled = true;
    control.spacing = std::get<0>(*wall);
    control.decay = std::get<1>(*wall);
    control.relaxation = std::get<2>(*wall);
    control.max_source = std::get<3>(*wall);
  }
  return control;
}

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
  if (name == "explicit") {
    return og::DistributionType::Explicit;
  }
  throw std::invalid_argument("unknown distribution type '" + name + "'");
}

og::DistributionSpec make_spec(const std::string& type, std::optional<double> begin_spacing,
                               std::optional<double> end_spacing,
                               std::optional<double> growth_ratio,
                               const std::optional<std::vector<double>>& positions) {
  og::DistributionSpec spec;
  spec.type = distribution_type(type);
  spec.begin_spacing = begin_spacing;
  spec.end_spacing = end_spacing;
  spec.growth_ratio = growth_ratio;
  if (positions) {
    spec.positions = *positions;
  }
  return spec;
}

// (start, end, n_points, type, begin_spacing, end_spacing, growth_ratio, positions)
using SubTuple =
    std::tuple<double, double, int, std::string, std::optional<double>, std::optional<double>,
               std::optional<double>, std::optional<std::vector<double>>>;
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

void bind_3d(nb::module_& m);

NB_MODULE(_core, m) {  // NOLINT(bugprone-signed-bitwise): inside the nanobind macro
  m.doc() = "OpenGrid C++ core";
  // Qt/VTK keep Python callbacks (and through them curve objects) alive until after the
  // interpreter has torn this module down, which nanobind reports as leaks at exit. These are
  // shutdown-order artefacts, not reference-counting errors; the API-only test suite runs
  // leak-free with the warnings enabled (OG_LEAK_WARNINGS=1 re-enables them).
  nb::set_leak_warnings(std::getenv("OG_LEAK_WARNINGS") != nullptr);
  m.def("version", &og::version, "Version string of the C++ core library.");
  m.def("add", &og::add, nb::arg("a"), nb::arg("b"), "Add two 64-bit integers (build smoke test).");
  m.def("thread_count", &og::thread_count, "Threads used by parallel loops (1 without OpenMP).");
  m.def("set_thread_count", &og::set_thread_count, nb::arg("n"),
        "Set the number of threads of parallel loops (results do not depend on it).");
  m.def("openmp_enabled", &og::openmp_enabled, "Whether the core was built with OpenMP.");
  m.def(
      "gpu_info",
      [] {
        const og::GpuInfo info = og::gpu_info();
        std::map<std::string, std::string> out;
        out["built"] = std::to_string(static_cast<int>(info.built));
        out["available"] = std::to_string(static_cast<int>(info.available));
        out["platform"] = info.platform;
        out["device"] = info.device;
        out["version"] = info.version;
        out["compute_units"] = std::to_string(info.compute_units);
        out["float64"] = std::to_string(static_cast<int>(info.float64));
        out["reason"] = info.reason;
        return out;
      },
      "GPU compute (OpenCL): built, available (a device passed the self-test), platform, "
      "device, version, compute_units, float64, reason (all strings).");
  m.def(
      "set_gpu_mode",
      [](const std::string& mode) {
        if (mode == "off") {
          og::set_gpu_mode(og::GpuMode::Off);
        } else if (mode == "auto") {
          og::set_gpu_mode(og::GpuMode::Auto);
        } else if (mode == "on") {
          og::set_gpu_mode(og::GpuMode::On);
        } else {
          throw std::invalid_argument("GPU mode must be 'off', 'auto' or 'on'");
        }
      },
      nb::arg("mode"),
      "GPU use: 'off', 'auto' (devices that passed the self-test, grids with at least "
      "gpu_min_nodes nodes) or 'on'.");
  m.def(
      "gpu_mode",
      [] {
        const og::GpuMode mode = og::gpu_mode();
        if (mode == og::GpuMode::Off) {
          return std::string("off");
        }
        if (mode == og::GpuMode::On) {
          return std::string("on");
        }
        return std::string("auto");
      },
      "The GPU mode: 'off', 'auto' or 'on'.");
  m.def(
      "nodes_1d",
      [](const std::string& set, int count) {
        return to_numpy(og::nodes_1d(og::node_set_from(set), count));
      },
      nb::arg("set"), nb::arg("count"),
      "Points of a 1D node set on [-1, 1]: 'equispaced', 'gll' or 'gauss'.");
  m.def(
      "weights_1d",
      [](const std::string& set, int count) {
        return to_numpy(og::weights_1d(og::node_set_from(set), count));
      },
      nb::arg("set"), nb::arg("count"), "Quadrature weights of 'gll' or 'gauss' points.");
  m.def(
      "high_order_elements",
      [](const nb::ndarray<const double, nb::c_contig, nb::device::cpu>& grid, int order,
         const std::string& nodes) {
        int ni = 0;
        int nj = 0;
        int nk = 1;
        if (grid.ndim() == 3 && grid.shape(2) == 3) {
          nj = static_cast<int>(grid.shape(0));
          ni = static_cast<int>(grid.shape(1));
        } else if (grid.ndim() == 4 && grid.shape(3) == 3) {
          nk = static_cast<int>(grid.shape(0));
          nj = static_cast<int>(grid.shape(1));
          ni = static_cast<int>(grid.shape(2));
        } else {
          throw std::invalid_argument("grid must be (nj, ni, 3) or (nk, nj, ni, 3)");
        }
        const std::size_t count = static_cast<std::size_t>(ni) * static_cast<std::size_t>(nj) *
                                  static_cast<std::size_t>(nk);
        std::vector<og::Vec3> points(count);
        const double* data = grid.data();
        for (std::size_t g = 0; g < count; ++g) {
          points[g] = {data[3 * g], data[3 * g + 1], data[3 * g + 2]};
        }
        og::HighOrderElements e;
        {
          const nb::gil_scoped_release release;
          e = og::high_order_from_grid(points, ni, nj, nk, order, og::node_set_from(nodes));
        }
        const auto n1 = static_cast<std::size_t>(order) + 1;
        std::size_t per = n1 * n1;
        if (e.dim == 3) {
          per *= n1;
        }
        auto* xyz = new double[e.xyz.size()];
        std::copy(e.xyz.begin(), e.xyz.end(), xyz);
        const nb::capsule owner(xyz, [](void* p) noexcept { delete[] static_cast<double*>(p); });
        nb::ndarray<nb::numpy, double, nb::ndim<3>> nodes_out(
            xyz, {static_cast<std::size_t>(e.count), per, 3}, owner);
        return std::make_tuple(nodes_out, to_numpy(e.min_jacobian), to_numpy(e.max_jacobian), e.ne);
      },
      nb::arg("grid"), nb::arg("order"), nb::arg("nodes") = "equispaced",
      "Curved order-p elements of a structured grid ((nj, ni, 3) or (nk, nj, ni, 3); every "
      "interval count a multiple of p): (nodes (elements, (p + 1)^d, 3) in tensor order with i "
      "fastest, min Jacobian, max Jacobian (per element, at the (p + 2)^d GLL points), "
      "elements per direction (i, j, k)). nodes: 'equispaced', 'gll' or 'gauss'.");
  m.attr("gpu_min_nodes") = og::kGpuMinNodes;
  m.attr("gpu_tolerance") = og::kGpuTolerance;

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
         std::optional<double> end_spacing, std::optional<double> growth_ratio,
         const std::optional<std::vector<double>>& positions) {
        const og::Distribution d = og::distribute(
            n_points, make_spec(type, begin_spacing, end_spacing, growth_ratio, positions));
        return std::make_tuple(to_numpy(d.s), d.warnings);
      },
      nb::arg("n_points"), nb::arg("type") = "uniform", nb::arg("begin_spacing") = nb::none(),
      nb::arg("end_spacing") = nb::none(), nb::arg("growth_ratio") = nb::none(),
      nb::arg("positions") = nb::none(),
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
          nb::arg("u"))
      .def("locate_length", &og::CurveChain::locate_length, nb::arg("s"),
           "(piece index, curve parameter) of the point at arc length s.");

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
        for (const auto& [start, end, n_points, type, begin, finish, ratio, positions] : subs) {
          og::SubConnectorSpec spec;
          spec.start = start;
          spec.end = end;
          spec.n_points = n_points;
          spec.distribution = make_spec(type, begin, finish, ratio, positions);
          specs.push_back(spec);
        }
        const og::ConnectorNodes nodes = og::connector_nodes(chain, spacing_mode, specs);
        return std::make_tuple(to_numpy(nodes.points), to_numpy(nodes.arc_positions),
                               nodes.warnings, nodes.max_stretching_ratio);
      },
      nb::arg("chain"), nb::arg("mode"), nb::arg("subs"),
      "Connector nodes: returns (points (n, 3), arc_positions (n,), warnings, max ratio).");

  m.def(
      "tfi",
      [](const PointsIn& bottom, const PointsIn& right, const PointsIn& top, const PointsIn& left,
         const std::string& method, bool hermite_bottom, bool hermite_top, double strength) {
        og::Boundary2D b{points_from(bottom), points_from(right), points_from(top),
                         points_from(left)};
        if (method == "linear") {
          return grid_to(og::tfi_linear(b));
        }
        if (method == "hermite") {
          og::HermiteOptions options;
          options.bottom = hermite_bottom;
          options.top = hermite_top;
          options.strength = strength;
          return grid_to(og::tfi_hermite(b, options));
        }
        throw std::invalid_argument("unknown TFI method '" + method + "'");
      },
      nb::arg("bottom"), nb::arg("right"), nb::arg("top"), nb::arg("left"),
      nb::arg("method") = "linear", nb::arg("hermite_bottom") = true,
      nb::arg("hermite_top") = false, nb::arg("strength") = 1.0,
      "Transfinite interpolation; boundaries ordered by increasing index; returns (nj, ni, 3).");

  nb::class_<og::SmoothingProgress>(m, "SmoothingProgress",
                                    "Iteration count, residual history and cancel flag.")
      .def(nb::init<>())
      .def("cancel", &og::SmoothingProgress::cancel)
      .def_prop_ro("cancelled", &og::SmoothingProgress::cancelled)
      .def_prop_ro("iteration", &og::SmoothingProgress::iteration)
      .def("residuals", [](const og::SmoothingProgress& p) { return to_numpy(p.residuals()); });

  m.def(
      "smooth_elliptic",
      [](const GridIn& grid, const std::string& control,
         const std::vector<std::tuple<bool, double, double, double>>& lines,
         const std::vector<std::tuple<double, double, double, double>>& points,
         const std::optional<WallTuple>& bottom_wall, const std::optional<WallTuple>& top_wall,
         const std::string& solver, double omega, int max_iterations, double tolerance,
         og::SmoothingProgress* progress, int acceleration) {
        og::EllipticOptions options;
        if (control == "laplace") {
          options.control = og::ControlType::Laplace;
        } else if (control == "thomas_middlecoff") {
          options.control = og::ControlType::ThomasMiddlecoff;
        } else if (control == "ttm") {
          options.control = og::ControlType::Ttm;
        } else {
          throw std::invalid_argument("unknown control '" + control + "'");
        }
        for (const auto& [xi_line, position, amplitude, decay] : lines) {
          options.lines.push_back({xi_line, position, amplitude, decay});
        }
        for (const auto& [xi, eta, amplitude, decay] : points) {
          options.points.push_back({xi, eta, amplitude, decay});
        }
        options.bottom = wall_from(bottom_wall);
        options.top = wall_from(top_wall);
        if (solver == "point_sor") {
          options.solver = og::EllipticSolver::PointSor;
        } else if (solver == "line_sor") {
          options.solver = og::EllipticSolver::LineSor;
        } else {
          throw std::invalid_argument("unknown solver '" + solver + "'");
        }
        options.omega = omega;
        options.max_iterations = max_iterations;
        options.tolerance = tolerance;
        options.acceleration = acceleration;
        const og::Grid2D initial = grid_from(grid);
        og::EllipticResult result;
        {
          const nb::gil_scoped_release release;
          result = og::smooth_elliptic(initial, options, progress);
        }
        return std::make_tuple(grid_to(result.grid), to_numpy(result.residuals), result.iterations,
                               result.converged, result.cancelled, result.warnings);
      },
      nb::arg("grid"), nb::arg("control") = "laplace", nb::arg("lines") = nb::list(),
      nb::arg("points") = nb::list(), nb::arg("bottom_wall") = nb::none(),
      nb::arg("top_wall") = nb::none(), nb::arg("solver") = "line_sor", nb::arg("omega") = 1.0,
      nb::arg("max_iterations") = 1000, nb::arg("tolerance") = 1e-6,
      nb::arg("progress") = nb::none(), nb::arg("acceleration") = 0,
      "Elliptic smoothing (GIL released); returns (grid, residuals, iterations, converged, "
      "cancelled, warnings).");

  m.def(
      "smooth_multiblock",
      [](const PointsIn& points,
         const nb::ndarray<const std::int64_t, nb::shape<-1, 4>, nb::c_contig, nb::device::cpu>&
             quads,
         const nb::ndarray<const std::uint8_t, nb::ndim<1>, nb::c_contig, nb::device::cpu>& fixed,
         const std::vector<std::shared_ptr<og::CurveChain>>& chains,
         const std::vector<std::tuple<std::int64_t, std::size_t>>& slides,
         const std::vector<std::tuple<std::int64_t, double, double, double, double, int, double>>&
             walls,
         double omega, int max_iterations, double tolerance, int acceleration,
         og::SmoothingProgress* progress, double slide_limit, double absolute_tolerance) {
        const std::size_t n = points.shape(0);
        if (fixed.shape(0) != n) {
          throw std::invalid_argument("fixed must have one entry per point");
        }
        og::MultiblockInput input;
        input.xy.resize(2 * n);
        for (std::size_t k = 0; k < n; ++k) {
          input.xy[2 * k] = points(k, 0);
          input.xy[2 * k + 1] = points(k, 1);
        }
        if (n > 0) {
          input.z = points(0, 2);
        }
        input.quads.resize(quads.shape(0));
        for (std::size_t q = 0; q < quads.shape(0); ++q) {
          for (std::size_t c = 0; c < 4; ++c) {
            input.quads[q][c] = quads(q, c);
          }
        }
        input.fixed.assign(fixed.data(), fixed.data() + n);
        for (const auto& chain : chains) {
          input.chains.push_back(chain);
        }
        for (const auto& [node, chain] : slides) {
          input.slides.push_back({node, chain});
        }
        for (const auto& [node, spacing, decay, relaxation, max_source, layers, growth] : walls) {
          input.walls.push_back({node, spacing, decay, relaxation, max_source, layers, growth});
        }
        og::MultiblockOptions options;
        options.omega = omega;
        options.max_iterations = max_iterations;
        options.tolerance = tolerance;
        options.acceleration = acceleration;
        options.slide_limit = slide_limit;
        options.absolute_tolerance = absolute_tolerance;
        og::MultiblockResult result;
        {
          const nb::gil_scoped_release release;
          result = og::smooth_multiblock(input, options, progress);
        }
        std::vector<og::Vec3> out(n);
        for (std::size_t k = 0; k < n; ++k) {
          out[k] = {result.xy[2 * k], result.xy[2 * k + 1], input.z};
        }
        auto* classes = new std::uint8_t[n];
        std::copy(result.node_class.begin(), result.node_class.end(), classes);
        const nb::capsule owner(
            classes, [](void* ptr) noexcept { delete[] static_cast<std::uint8_t*>(ptr); });
        nb::ndarray<nb::numpy, std::uint8_t, nb::ndim<1>> node_class(classes, {n}, owner);
        return std::make_tuple(to_numpy(out), node_class, to_numpy(result.residuals),
                               result.iterations, result.converged, result.cancelled, result.lines,
                               result.colours, result.warnings);
      },
      nb::arg("points"), nb::arg("quads"), nb::arg("fixed"), nb::arg("chains"), nb::arg("slides"),
      nb::arg("walls"), nb::arg("omega") = 1.0, nb::arg("max_iterations") = 2000,
      nb::arg("tolerance") = 1e-6, nb::arg("acceleration") = 20, nb::arg("progress") = nb::none(),
      nb::arg("slide_limit") = 0.5, nb::arg("absolute_tolerance") = 0.0,
      "Elliptic smoothing of a planar multi-block grid as one grid (GIL released). quads are "
      "counter-clockwise node ids; slides are (node, chain index; a node moves at most "
      "slide_limit of its initial boundary intervals); walls are (node, spacing "
      "(<= 0: initial), decay, relaxation, max_source, layers, growth). Returns (points, "
      "node_class (0 fixed, 1 regular, 2 irregular, 3 sliding, 4 wall layer), residuals, "
      "iterations, converged, "
      "cancelled, lines, colours, warnings).");

  m.def(
      "quality_2d",
      [](const GridIn& grid) {
        const og::Grid2D g = grid_from(grid);
        const og::Quality2D q = og::quality_2d(g);
        const int ci = g.ni - 1;
        const int cj = g.nj - 1;
        std::map<std::string, CellsOut> out;
        out.emplace("scaled_jacobian", cells_to(q.scaled_jacobian, ci, cj));
        out.emplace("skewness", cells_to(q.skewness, ci, cj));
        out.emplace("aspect_ratio", cells_to(q.aspect_ratio, ci, cj));
        out.emplace("min_angle", cells_to(q.min_angle, ci, cj));
        out.emplace("max_angle", cells_to(q.max_angle, ci, cj));
        out.emplace("non_orthogonality", cells_to(q.non_orthogonality, ci, cj));
        out.emplace("area_ratio", cells_to(q.area_ratio, ci, cj));
        out.emplace("area", cells_to(q.area, ci, cj));
        return out;
      },
      nb::arg("grid"), "Per-cell quality metrics, each of shape (nj - 1, ni - 1).");

  m.def("cgns_available", &og::cgns_available, "True if built with the CGNS library.");
  m.def(
      "write_cgns_2d",
      [](const std::string& path, const std::string& base_name,
         const std::vector<std::tuple<std::string, int, int, std::vector<double>,
                                      std::vector<double>, std::string>>& zones,
         const std::vector<std::tuple<
             std::string, std::vector<double>, std::vector<double>,
             nb::ndarray<const std::int64_t, nb::shape<-1, 3>, nb::c_contig, nb::device::cpu>,
             nb::ndarray<const std::int64_t, nb::shape<-1, 4>, nb::c_contig, nb::device::cpu>,
             std::string>>& unstructured_zones,
         const std::vector<std::tuple<std::string, std::string>>& families,
         const std::vector<std::tuple<int, std::string, std::string, std::array<int, 4>>>& bcs,
         const std::vector<std::tuple<int, std::string, std::string,
                                      std::vector<std::array<std::int64_t, 2>>>>& edge_bcs,
         const std::vector<std::tuple<int, std::string, int, std::array<int, 4>, std::array<int, 4>,
                                      std::array<int, 2>, std::optional<std::array<double, 2>>>>&
             connections,
         const std::vector<
             std::tuple<int, std::string, int, std::vector<std::int64_t>, std::vector<std::int64_t>,
                        std::optional<std::array<double, 2>>>>& vertex_connections) {
        std::vector<og::CgnsZone> z;
        z.reserve(zones.size());
        for (const auto& [name, ni, nj, x, y, family] : zones) {
          z.push_back({name, ni, nj, x, y, family});
        }
        std::vector<og::CgnsUnstructuredZone> uz;
        uz.reserve(unstructured_zones.size());
        for (const auto& [name, x, y, triangles, quads, family] : unstructured_zones) {
          og::CgnsUnstructuredZone zone{name, x, y, {}, {}, family};
          for (std::size_t k = 0; k < triangles.shape(0); ++k) {
            zone.triangles.push_back({triangles(k, 0), triangles(k, 1), triangles(k, 2)});
          }
          for (std::size_t k = 0; k < quads.shape(0); ++k) {
            zone.quads.push_back({quads(k, 0), quads(k, 1), quads(k, 2), quads(k, 3)});
          }
          uz.push_back(std::move(zone));
        }
        std::vector<og::CgnsFamily> f;
        f.reserve(families.size());
        for (const auto& [name, bc_type] : families) {
          f.push_back({name, bc_type});
        }
        std::vector<og::CgnsBoundary> b;
        b.reserve(bcs.size());
        for (const auto& [zone, name, family, range] : bcs) {
          b.push_back({zone, name, family, range});
        }
        std::vector<og::CgnsEdgeBoundary> eb;
        eb.reserve(edge_bcs.size());
        for (const auto& [zone, name, family, edges] : edge_bcs) {
          eb.push_back({zone, name, family, edges});
        }
        std::vector<og::CgnsConnection> c;
        for (const auto& [zone, name, donor, range, donor_range, transform, translation] :
             connections) {
          og::CgnsConnection connection{zone, name, donor, range, donor_range, transform};
          if (translation) {
            connection.periodic = true;
            connection.translation = *translation;
          }
          c.push_back(connection);
        }
        std::vector<og::CgnsVertexConnection> vc;
        for (const auto& [zone, name, donor, points, donor_points, translation] :
             vertex_connections) {
          og::CgnsVertexConnection connection{zone, name, donor, points, donor_points};
          if (translation) {
            connection.periodic = true;
            connection.translation = *translation;
          }
          vc.push_back(connection);
        }
        const nb::gil_scoped_release release;
        og::write_cgns_2d(path, base_name, z, uz, f, b, eb, c, vc);
      },
      nb::arg("path"), nb::arg("base_name"), nb::arg("zones"), nb::arg("unstructured_zones"),
      nb::arg("families"), nb::arg("boundaries"), nb::arg("edge_boundaries"),
      nb::arg("connections"), nb::arg("vertex_connections"),
      "Write a 2D CGNS file with structured and unstructured zones. unstructured_zones: (name, x, "
      "y, triangles (m, 3), quads (k, 4), family), 0-based vertex ids; edge_boundaries: "
      "(unstructured zone index, name, family, [(v0, v1), ...]); vertex_connections: (zone, name, "
      "donor zone, points, donor points, periodic translation or None) with global zone indices "
      "(structured first) and 1-based points ((i, j) pairs for structured zones). Other "
      "arguments as write_cgns_structured_2d.");

  m.def(
      "write_cgns_structured_2d",
      [](const std::string& path, const std::string& base_name,
         const std::vector<std::tuple<std::string, int, int, std::vector<double>,
                                      std::vector<double>, std::string>>& zones,
         const std::vector<std::tuple<std::string, std::string>>& families,
         const std::vector<std::tuple<int, std::string, std::string, std::array<int, 4>>>& bcs,
         const std::vector<std::tuple<int, std::string, int, std::array<int, 4>, std::array<int, 4>,
                                      std::array<int, 2>, std::optional<std::array<double, 2>>>>&
             connections) {
        std::vector<og::CgnsZone> z;
        z.reserve(zones.size());
        for (const auto& [name, ni, nj, x, y, family] : zones) {
          z.push_back({name, ni, nj, x, y, family});
        }
        std::vector<og::CgnsFamily> f;
        f.reserve(families.size());
        for (const auto& [name, bc_type] : families) {
          f.push_back({name, bc_type});
        }
        std::vector<og::CgnsBoundary> b;
        b.reserve(bcs.size());
        for (const auto& [zone, name, family, range] : bcs) {
          b.push_back({zone, name, family, range});
        }
        std::vector<og::CgnsConnection> c;
        c.reserve(connections.size());
        for (const auto& [zone, name, donor, range, donor_range, transform, translation] :
             connections) {
          og::CgnsConnection connection{zone, name, donor, range, donor_range, transform};
          if (translation) {
            connection.periodic = true;
            connection.translation = *translation;
          }
          c.push_back(connection);
        }
        const nb::gil_scoped_release release;
        og::write_cgns_structured_2d(path, base_name, z, f, b, c);
      },
      nb::arg("path"), nb::arg("base_name"), nb::arg("zones"), nb::arg("families"),
      nb::arg("boundaries"), nb::arg("connections"),
      "Write a structured 2D CGNS file. zones: (name, ni, nj, x, y, family); families: (name, "
      "bc_type or ''); boundaries: (zone index, name, family, (imin, jmin, imax, jmax)); "
      "connections: (zone, name, donor zone, range, donor range, transform, periodic translation "
      "or None). 1-based ranges.");

  m.def(
      "triangulate",
      [](const nb::ndarray<const double, nb::shape<-1, 2>, nb::c_contig, nb::device::cpu>& points,
         const nb::ndarray<const std::int64_t, nb::shape<-1, 2>, nb::c_contig, nb::device::cpu>&
             segments,
         bool refine, double min_angle, double growth_rate, std::optional<double> max_size,
         const std::vector<std::tuple<double, double, double, double, double, double>>& sources,
         std::int64_t max_points, const std::string& algorithm) {
        if (algorithm != "delaunay" && algorithm != "advancing_front") {
          throw std::invalid_argument("algorithm must be 'delaunay' or 'advancing_front'");
        }
        og::TriangulationInput input;
        for (std::size_t k = 0; k < points.shape(0); ++k) {
          input.points.push_back({points(k, 0), points(k, 1)});
        }
        for (std::size_t k = 0; k < segments.shape(0); ++k) {
          input.segments.push_back({segments(k, 0), segments(k, 1)});
        }
        og::RefinementOptions options;
        options.refine = refine;
        options.min_angle_deg = min_angle;
        options.size.growth_rate = growth_rate;
        options.size.max_size = max_size;
        for (const auto& [x0, y0, x1, y1, size, radius] : sources) {
          options.size.sources.push_back({{x0, y0}, {x1, y1}, size, radius});
        }
        options.max_points = max_points;
        og::TriangleMesh mesh;
        {
          const nb::gil_scoped_release release;
          if (algorithm == "advancing_front") {
            mesh = og::advancing_front(input, options);
          } else {
            mesh = og::triangulate_domain(input, options);
          }
        }
        auto* xy = new double[mesh.points.size() * 2];
        for (std::size_t k = 0; k < mesh.points.size(); ++k) {
          xy[2 * k] = mesh.points[k][0];
          xy[2 * k + 1] = mesh.points[k][1];
        }
        auto* tri = new std::int64_t[mesh.triangles.size() * 3];
        for (std::size_t k = 0; k < mesh.triangles.size(); ++k) {
          for (std::size_t j = 0; j < 3; ++j) {
            tri[3 * k + j] = mesh.triangles[k][j];
          }
        }
        const nb::capsule own_xy(xy,
                                 [](void* ptr) noexcept { delete[] static_cast<double*>(ptr); });
        const nb::capsule own_tri(
            tri, [](void* ptr) noexcept { delete[] static_cast<std::int64_t*>(ptr); });
        return std::make_tuple(
            nb::ndarray<nb::numpy, double, nb::ndim<2>>(xy, {mesh.points.size(), 2}, own_xy),
            nb::ndarray<nb::numpy, std::int64_t, nb::ndim<2>>(tri, {mesh.triangles.size(), 3},
                                                              own_tri),
            mesh.rejected, mesh.warnings);
      },
      nb::arg("points"), nb::arg("segments"), nb::arg("refine") = true, nb::arg("min_angle") = 25.0,
      nb::arg("growth_rate") = 1.2, nb::arg("max_size") = nb::none(),
      nb::arg("sources") = nb::list(), nb::arg("max_points") = 20'000'000,
      nb::arg("algorithm") = "delaunay",
      "Triangulation of a domain bounded by oriented segments (domain on the left): "
      "algorithm='delaunay' (constrained Delaunay + refinement) or 'advancing_front' (refine and "
      "min_angle are ignored). Segments oriented with the domain on the "
      "left. sources: (x0, y0, x1, y1, size, radius). Returns (points (n, 2), triangles (m, 3), "
      "rejected, warnings).");

  m.def(
      "quality_triangles",
      [](const nb::ndarray<const double, nb::shape<-1, 2>, nb::c_contig, nb::device::cpu>& points,
         const nb::ndarray<const std::int64_t, nb::shape<-1, 3>, nb::c_contig, nb::device::cpu>&
             triangles) {
        std::vector<std::array<double, 2>> p;
        p.reserve(points.shape(0));
        for (std::size_t k = 0; k < points.shape(0); ++k) {
          p.push_back({points(k, 0), points(k, 1)});
        }
        std::vector<std::array<std::int64_t, 3>> t;
        t.reserve(triangles.shape(0));
        for (std::size_t k = 0; k < triangles.shape(0); ++k) {
          t.push_back({triangles(k, 0), triangles(k, 1), triangles(k, 2)});
        }
        const og::Quality2D q = og::quality_triangles(p, t);
        const auto n = static_cast<int>(t.size());
        std::map<std::string, CellsOut> out;
        out.emplace("scaled_jacobian", cells_to(q.scaled_jacobian, n, 1));
        out.emplace("skewness", cells_to(q.skewness, n, 1));
        out.emplace("aspect_ratio", cells_to(q.aspect_ratio, n, 1));
        out.emplace("min_angle", cells_to(q.min_angle, n, 1));
        out.emplace("max_angle", cells_to(q.max_angle, n, 1));
        out.emplace("non_orthogonality", cells_to(q.non_orthogonality, n, 1));
        out.emplace("area_ratio", cells_to(q.area_ratio, n, 1));
        out.emplace("area", cells_to(q.area, n, 1));
        return out;
      },
      nb::arg("points"), nb::arg("triangles"),
      "Per-triangle quality metrics, each of shape (1, m) (same keys as quality_2d).");

  m.def(
      "quality_cells",
      [](const nb::ndarray<const double, nb::shape<-1, 2>, nb::c_contig, nb::device::cpu>& points,
         const nb::ndarray<const std::int64_t, nb::shape<-1, 3>, nb::c_contig, nb::device::cpu>&
             triangles,
         const nb::ndarray<const std::int64_t, nb::shape<-1, 4>, nb::c_contig, nb::device::cpu>&
             quads) {
        std::vector<std::array<double, 2>> p;
        p.reserve(points.shape(0));
        for (std::size_t k = 0; k < points.shape(0); ++k) {
          p.push_back({points(k, 0), points(k, 1)});
        }
        std::vector<std::array<std::int64_t, 3>> t;
        t.reserve(triangles.shape(0));
        for (std::size_t k = 0; k < triangles.shape(0); ++k) {
          t.push_back({triangles(k, 0), triangles(k, 1), triangles(k, 2)});
        }
        std::vector<std::array<std::int64_t, 4>> qd;
        qd.reserve(quads.shape(0));
        for (std::size_t k = 0; k < quads.shape(0); ++k) {
          qd.push_back({quads(k, 0), quads(k, 1), quads(k, 2), quads(k, 3)});
        }
        const og::Quality2D q = og::quality_cells(p, t, qd);
        const auto n = static_cast<int>(t.size() + qd.size());
        std::map<std::string, CellsOut> out;
        out.emplace("scaled_jacobian", cells_to(q.scaled_jacobian, n, 1));
        out.emplace("skewness", cells_to(q.skewness, n, 1));
        out.emplace("aspect_ratio", cells_to(q.aspect_ratio, n, 1));
        out.emplace("min_angle", cells_to(q.min_angle, n, 1));
        out.emplace("max_angle", cells_to(q.max_angle, n, 1));
        out.emplace("non_orthogonality", cells_to(q.non_orthogonality, n, 1));
        out.emplace("area_ratio", cells_to(q.area_ratio, n, 1));
        out.emplace("area", cells_to(q.area, n, 1));
        return out;
      },
      nb::arg("points"), nb::arg("triangles"), nb::arg("quads"),
      "Per-cell quality metrics of a triangle/quad mesh (triangles first, then quads), each of "
      "shape (1, m) (same keys as quality_2d).");

  m.def(
      "mesh_domain",
      [](const nb::ndarray<const double, nb::shape<-1, 2>, nb::c_contig, nb::device::cpu>& points,
         const nb::ndarray<const std::int64_t, nb::shape<-1, 2>, nb::c_contig, nb::device::cpu>&
             segments,
         const std::string& algorithm, const std::string& element_type, bool refine,
         double min_angle, double growth_rate, std::optional<double> max_size,
         const std::vector<std::tuple<double, double, double, double, double, double>>& sources,
         std::int64_t max_points, const std::string& smoothing, int smoothing_iterations,
         bool edge_swap, const std::vector<bool>& wall, double first_height,
         double inflation_growth, int layers, std::optional<double> total_thickness,
         bool transition) {
        og::DomainMeshOptions options;
        options.wall = wall;
        options.inflation.first_height = first_height;
        options.inflation.growth_rate = inflation_growth;
        options.inflation.max_layers = layers;
        options.inflation.total_thickness = total_thickness;
        options.inflation.transition = transition;
        if (algorithm == "delaunay") {
          options.algorithm = og::TriangleAlgorithm::delaunay;
        } else if (algorithm == "advancing_front") {
          options.algorithm = og::TriangleAlgorithm::advancing_front;
        } else {
          throw std::invalid_argument("algorithm must be 'delaunay' or 'advancing_front'");
        }
        if (element_type == "triangles") {
          options.element_type = og::ElementType::triangles;
        } else if (element_type == "quad_dominant") {
          options.element_type = og::ElementType::quad_dominant;
        } else if (element_type == "quads") {
          options.element_type = og::ElementType::quads;
        } else {
          throw std::invalid_argument(
              "element_type must be 'triangles', 'quad_dominant' or 'quads'");
        }
        if (smoothing == "none") {
          options.smoothing.method = og::SmoothingMethod::none;
        } else if (smoothing == "laplacian") {
          options.smoothing.method = og::SmoothingMethod::laplacian;
        } else if (smoothing == "length_weighted") {
          options.smoothing.method = og::SmoothingMethod::length_weighted;
        } else if (smoothing == "angle_based") {
          options.smoothing.method = og::SmoothingMethod::angle_based;
        } else {
          throw std::invalid_argument(
              "smoothing must be 'none', 'laplacian', 'length_weighted' or 'angle_based'");
        }
        options.smoothing.iterations = smoothing_iterations;
        options.smoothing.edge_swap = edge_swap;
        og::TriangulationInput input;
        for (std::size_t k = 0; k < points.shape(0); ++k) {
          input.points.push_back({points(k, 0), points(k, 1)});
        }
        for (std::size_t k = 0; k < segments.shape(0); ++k) {
          input.segments.push_back({segments(k, 0), segments(k, 1)});
        }
        options.refinement.refine = refine;
        options.refinement.min_angle_deg = min_angle;
        options.refinement.size.growth_rate = growth_rate;
        options.refinement.size.max_size = max_size;
        for (const auto& [x0, y0, x1, y1, size, radius] : sources) {
          options.refinement.size.sources.push_back({{x0, y0}, {x1, y1}, size, radius});
        }
        options.refinement.max_points = max_points;
        og::MixedMesh mesh;
        {
          const nb::gil_scoped_release release;
          mesh = og::mesh_domain(input, options);
        }
        auto* xy = new double[mesh.points.size() * 2];
        for (std::size_t k = 0; k < mesh.points.size(); ++k) {
          xy[2 * k] = mesh.points[k][0];
          xy[2 * k + 1] = mesh.points[k][1];
        }
        auto* tri = new std::int64_t[mesh.triangles.size() * 3];
        for (std::size_t k = 0; k < mesh.triangles.size(); ++k) {
          for (std::size_t j = 0; j < 3; ++j) {
            tri[3 * k + j] = mesh.triangles[k][j];
          }
        }
        auto* quad = new std::int64_t[mesh.quads.size() * 4];
        for (std::size_t k = 0; k < mesh.quads.size(); ++k) {
          for (std::size_t j = 0; j < 4; ++j) {
            quad[4 * k + j] = mesh.quads[k][j];
          }
        }
        const nb::capsule own_xy(xy,
                                 [](void* ptr) noexcept { delete[] static_cast<double*>(ptr); });
        const nb::capsule own_tri(
            tri, [](void* ptr) noexcept { delete[] static_cast<std::int64_t*>(ptr); });
        const nb::capsule own_quad(
            quad, [](void* ptr) noexcept { delete[] static_cast<std::int64_t*>(ptr); });
        return std::make_tuple(
            nb::ndarray<nb::numpy, double, nb::ndim<2>>(xy, {mesh.points.size(), 2}, own_xy),
            nb::ndarray<nb::numpy, std::int64_t, nb::ndim<2>>(tri, {mesh.triangles.size(), 3},
                                                              own_tri),
            nb::ndarray<nb::numpy, std::int64_t, nb::ndim<2>>(quad, {mesh.quads.size(), 4},
                                                              own_quad),
            mesh.rejected, mesh.warnings, mesh.layers);
      },
      nb::arg("points"), nb::arg("segments"), nb::arg("algorithm") = "delaunay",
      nb::arg("element_type") = "triangles", nb::arg("refine") = true, nb::arg("min_angle") = 25.0,
      nb::arg("growth_rate") = 1.2, nb::arg("max_size") = nb::none(),
      nb::arg("sources") = nb::list(), nb::arg("max_points") = 20'000'000,
      nb::arg("smoothing") = "none", nb::arg("smoothing_iterations") = 5,
      nb::arg("edge_swap") = true, nb::arg("wall") = std::vector<bool>{},
      nb::arg("first_height") = 1e-3, nb::arg("inflation_growth") = 1.2, nb::arg("layers") = 30,
      nb::arg("total_thickness") = nb::none(), nb::arg("transition") = true,
      "Mesh a 2D domain bounded by oriented segments (domain on the left): triangulation "
      "(algorithm), element type ('triangles', 'quad_dominant', 'quads'), smoothing ('none', "
      "'laplacian', 'length_weighted', 'angle_based') with optional edge swapping, untangling. "
      "wall: one flag per segment; flagged segments get inflation layers (first_height, "
      "inflation_growth, layers, total_thickness, transition). Returns (points (n, 2), triangles "
      "(m, 3), quads (k, 4), rejected, warnings, layers).");

  nb::class_<og::Nurbs, og::Curve>(m, "Nurbs")
      .def(
          "__init__",
          [](og::Nurbs* self, int degree, const std::vector<double>& knots,
             const std::vector<Array3>& control_points, const std::vector<double>& weights) {
            new (self) og::Nurbs(degree, knots, to_vec3(control_points), weights);
          },
          nb::arg("degree"), nb::arg("knots"), nb::arg("control_points"), nb::arg("weights"));

  bind_3d(m);
}
