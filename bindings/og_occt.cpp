// Optional module opengrid._occt: STEP / IGES import through OpenCASCADE (backends/occt).
// The model crosses the boundary as plain Python data (dicts, lists, NumPy arrays); entity
// creation happens in opengrid.io.cad.
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include "opengrid/backends/occt/cad_import.hpp"

namespace nb = nanobind;

namespace {

using Points = nb::ndarray<nb::numpy, double, nb::shape<-1, 3>>;
using Values = nb::ndarray<nb::numpy, double, nb::shape<-1>>;

Points points_array(const std::vector<std::array<double, 3>>& points) {
  auto* data = new std::vector<double>();
  data->reserve(points.size() * 3);
  for (const auto& p : points) {
    data->insert(data->end(), p.begin(), p.end());
  }
  nb::capsule owner(data, [](void* p) noexcept { delete static_cast<std::vector<double>*>(p); });
  return Points(data->data(), {points.size(), 3}, owner);
}

Values values_array(const std::vector<double>& values) {
  auto* data = new std::vector<double>(values);
  nb::capsule owner(data, [](void* p) noexcept { delete static_cast<std::vector<double>*>(p); });
  return Values(data->data(), {values.size()}, owner);
}

nb::dict curve_dict(const og::cad::BSplineCurveData& c) {
  nb::dict d;
  d["degree"] = c.degree;
  d["knots"] = values_array(c.knots);
  d["poles"] = points_array(c.poles);
  d["weights"] = values_array(c.weights);
  return d;
}

nb::dict surface_dict(const og::cad::BSplineSurfaceData& s) {
  nb::dict d;
  d["degree_u"] = s.degree_u;
  d["degree_v"] = s.degree_v;
  d["knots_u"] = values_array(s.knots_u);
  d["knots_v"] = values_array(s.knots_v);
  d["nu"] = s.nu;
  d["nv"] = s.nv;
  d["poles"] = points_array(s.poles);
  d["weights"] = values_array(s.weights);
  return d;
}

nb::dict read_cad(const std::string& path, double tolerance, double scale, bool heal, bool sew,
                  bool split_closed) {
  og::cad::CadOptions options;
  options.tolerance = tolerance;
  options.scale = scale;
  options.heal = heal;
  options.sew = sew;
  options.split_closed = split_closed;
  og::cad::CadModel model;
  {
    nb::gil_scoped_release release;
    model = og::cad::read_cad(path, options);
  }
  nb::list edges;
  for (const auto& e : model.edges) {
    nb::dict d;
    d["curve"] = curve_dict(e.curve);
    d["v0"] = e.v0;
    d["v1"] = e.v1;
    d["length"] = e.length;
    d["tolerance"] = e.tolerance;
    d["faces"] = e.faces;
    d["type"] = e.type;
    edges.append(d);
  }
  nb::list faces;
  for (const auto& f : model.faces) {
    nb::dict d;
    if (f.has_surface) {
      d["surface"] = surface_dict(f.surface);
    } else {
      d["surface"] = nb::none();
    }
    d["type"] = f.type;
    d["degenerate"] = f.degenerate;
    nb::list wires;
    for (const auto& w : f.wires) {
      nb::dict wd;
      wd["edges"] = w.edges;
      nb::list reversed;
      for (const bool r : w.reversed) {
        reversed.append(r);
      }
      wd["reversed"] = reversed;
      wd["outer"] = w.outer;
      wires.append(wd);
    }
    d["wires"] = wires;
    faces.append(d);
  }
  nb::dict result;
  result["vertices"] = points_array(model.vertices);
  result["edges"] = edges;
  result["faces"] = faces;
  result["solids"] = model.solids;
  result["shells"] = model.shells;
  result["max_tolerance"] = model.max_tolerance;
  result["log"] = model.log;
  return result;
}

}  // namespace

NB_MODULE(_occt, m) {  // NOLINT(bugprone-signed-bitwise): inside the nanobind macro
  m.doc() = "OpenGrid CAD import (OpenCASCADE backend)";
  m.def("read_cad", &read_cad, nb::arg("path"), nb::arg("tolerance") = 1e-6, nb::arg("scale") = 1.0,
        nb::arg("heal") = true, nb::arg("sew") = true, nb::arg("split_closed") = true,
        "Read a STEP / IGES file: healed topology with B-spline curves and surfaces "
        "(see backends/occt/include/opengrid/backends/occt/cad_import.hpp).");
}
