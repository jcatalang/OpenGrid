// Bindings of the 3D structured core: surfaces, surface TFI, block TFI, extrusion, 3D quality
// and the 3D CGNS writer. Grids cross the boundary as NumPy arrays of shape (nk, nj, ni, 3)
// (blocks) and (nj, ni, 3) (faces / surface grids).
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <array>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "opengrid/core/geom/curve.hpp"
#include "opengrid/core/geom/surface.hpp"
#include "opengrid/core/io/cgns_structured.hpp"
#include "opengrid/core/quality/quality3d.hpp"
#include "opengrid/core/structured/elliptic3d.hpp"
#include "opengrid/core/structured/extrude.hpp"
#include "opengrid/core/structured/hyperbolic.hpp"
#include "opengrid/core/structured/surface_elliptic.hpp"
#include "opengrid/core/structured/surface_tfi.hpp"
#include "opengrid/core/structured/tfi3d.hpp"
#include "opengrid/core/unstructured/inflation3d.hpp"
#include "opengrid/core/unstructured/optimize3d.hpp"

namespace nb = nanobind;

namespace {

using Face = nb::ndarray<const double, nb::ndim<3>, nb::c_contig, nb::device::cpu>;
using Block = nb::ndarray<const double, nb::ndim<4>, nb::c_contig, nb::device::cpu>;
using FaceOut = nb::ndarray<nb::numpy, double, nb::ndim<3>>;
using BlockOut = nb::ndarray<nb::numpy, double, nb::ndim<4>>;
using CellsOut = nb::ndarray<nb::numpy, double, nb::ndim<3>>;

og::Vec3 vec(const std::array<double, 3>& a) {
  return {a[0], a[1], a[2]};
}

og::Grid2D face_from(const Face& a) {
  if (a.shape(2) != 3) {
    throw std::invalid_argument("face grids must have shape (n2, n1, 3)");
  }
  og::Grid2D g(static_cast<int>(a.shape(1)), static_cast<int>(a.shape(0)));
  const double* d = a.data();
  for (std::size_t k = 0; k < g.x.size(); ++k) {
    g.x[k] = {d[3 * k], d[3 * k + 1], d[3 * k + 2]};
  }
  return g;
}

FaceOut face_to(const og::Grid2D& g) {
  auto* data = new double[g.x.size() * 3];
  for (std::size_t k = 0; k < g.x.size(); ++k) {
    data[3 * k] = g.x[k].x;
    data[3 * k + 1] = g.x[k].y;
    data[3 * k + 2] = g.x[k].z;
  }
  const nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<double*>(p); });
  return FaceOut(data, {static_cast<std::size_t>(g.nj), static_cast<std::size_t>(g.ni), 3}, owner);
}

og::Grid3D block_from(const Block& a) {
  if (a.shape(3) != 3) {
    throw std::invalid_argument("block grids must have shape (nk, nj, ni, 3)");
  }
  og::Grid3D g(static_cast<int>(a.shape(2)), static_cast<int>(a.shape(1)),
               static_cast<int>(a.shape(0)));
  const double* d = a.data();
  for (std::size_t k = 0; k < g.x.size(); ++k) {
    g.x[k] = {d[3 * k], d[3 * k + 1], d[3 * k + 2]};
  }
  return g;
}

BlockOut block_to(const og::Grid3D& g) {
  auto* data = new double[g.x.size() * 3];
  for (std::size_t k = 0; k < g.x.size(); ++k) {
    data[3 * k] = g.x[k].x;
    data[3 * k + 1] = g.x[k].y;
    data[3 * k + 2] = g.x[k].z;
  }
  const nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<double*>(p); });
  return BlockOut(data,
                  {static_cast<std::size_t>(g.nk), static_cast<std::size_t>(g.nj),
                   static_cast<std::size_t>(g.ni), 3},
                  owner);
}

CellsOut cells3_to(const std::vector<double>& values, int ci, int cj, int ck) {
  auto* data = new double[values.size()];
  std::copy(values.begin(), values.end(), data);
  const nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<double*>(p); });
  return CellsOut(
      data,
      {static_cast<std::size_t>(ck), static_cast<std::size_t>(cj), static_cast<std::size_t>(ci)},
      owner);
}

std::vector<og::Vec3> side(
    const nb::ndarray<const double, nb::shape<-1, 3>, nb::c_contig, nb::device::cpu>& a) {
  std::vector<og::Vec3> out;
  out.reserve(a.shape(0));
  for (std::size_t k = 0; k < a.shape(0); ++k) {
    out.emplace_back(a(k, 0), a(k, 1), a(k, 2));
  }
  return out;
}

template <std::size_t N>
nb::ndarray<nb::numpy, std::int64_t, nb::ndim<2>> index_array(
    const std::vector<std::array<std::int64_t, N>>& cells) {
  auto* data = new std::int64_t[cells.size() * N];
  for (std::size_t c = 0; c < cells.size(); ++c) {
    for (std::size_t k = 0; k < N; ++k) {
      data[c * N + k] = cells[c][k];
    }
  }
  const nb::capsule owner(data,
                          [](void* ptr) noexcept { delete[] static_cast<std::int64_t*>(ptr); });
  return {data, {cells.size(), N}, owner};
}

template <std::size_t N>
using Cells = nb::ndarray<const std::int64_t, nb::shape<-1, N>, nb::c_contig, nb::device::cpu>;

template <std::size_t N>
std::vector<std::array<std::int64_t, N>> index_input(const Cells<N>& a) {
  std::vector<std::array<std::int64_t, N>> out(a.shape(0));
  const std::int64_t* d = a.data();
  for (std::size_t c = 0; c < out.size(); ++c) {
    for (std::size_t k = 0; k < N; ++k) {
      out[c][k] = d[c * N + k];
    }
  }
  return out;
}

}  // namespace

void bind_3d(nb::module_& m) {  // NOLINT(misc-use-internal-linkage): called from og_core.cpp
  using Points = nb::ndarray<const double, nb::shape<-1, 3>, nb::c_contig, nb::device::cpu>;

  nb::class_<og::Surface>(m, "Surface", "Parametric surface x(u, v), (u, v) in [0, 1]^2.")
      .def("eval",
           [](const og::Surface& s, double u, double v) {
             const og::Vec3 p = s.eval(u, v);
             return std::array<double, 3>{p.x, p.y, p.z};
           })
      .def("normal",
           [](const og::Surface& s, double u, double v) {
             const og::Vec3 p = s.normal(u, v);
             return std::array<double, 3>{p.x, p.y, p.z};
           })
      .def(
          "project",
          [](const og::Surface& s, const std::array<double, 3>& p) {
            const og::SurfaceProjection q = s.project(vec(p));
            return std::make_tuple(q.u, q.v, std::array<double, 3>{q.point.x, q.point.y, q.point.z},
                                   q.distance);
          },
          "Closest point: (u, v, point, distance).")
      .def(
          "grid",
          [](const og::Surface& s, int nu, int nv) {
            og::Grid2D g(nu, nv);
            for (int j = 0; j < nv; ++j) {
              for (int i = 0; i < nu; ++i) {
                g.at(i, j) =
                    s.eval(static_cast<double>(i) / (nu - 1), static_cast<double>(j) / (nv - 1));
              }
            }
            return face_to(g);
          },
          nb::arg("nu"), nb::arg("nv"), "Points on a uniform parameter grid, shape (nv, nu, 3).")
      .def_prop_ro("size", &og::Surface::size);

  nb::class_<og::PlaneSurface, og::Surface>(m, "PlaneSurface")
      .def("__init__", [](og::PlaneSurface* self, const std::array<double, 3>& origin,
                          const std::array<double, 3>& a, const std::array<double, 3>& b) {
        new (self) og::PlaneSurface(vec(origin), vec(a), vec(b));
      });
  nb::class_<og::RuledSurface, og::Surface>(m, "RuledSurface")
      .def("__init__", [](og::RuledSurface* self, std::shared_ptr<og::Curve> c0,
                          std::shared_ptr<og::Curve> c1) {
        new (self) og::RuledSurface(std::move(c0), std::move(c1));
      });
  nb::class_<og::ExtrudedSurface, og::Surface>(m, "ExtrudedSurface")
      .def("__init__", [](og::ExtrudedSurface* self, std::shared_ptr<og::Curve> curve,
                          const std::array<double, 3>& direction) {
        new (self) og::ExtrudedSurface(std::move(curve), vec(direction));
      });
  nb::class_<og::RevolvedSurface, og::Surface>(m, "RevolvedSurface")
      .def("__init__", [](og::RevolvedSurface* self, std::shared_ptr<og::Curve> curve,
                          const std::array<double, 3>& point, const std::array<double, 3>& axis,
                          double angle0, double angle1) {
        new (self) og::RevolvedSurface(std::move(curve), vec(point), vec(axis), angle0, angle1);
      });
  nb::class_<og::NurbsSurface, og::Surface>(m, "NurbsSurface")
      .def("__init__",
           [](og::NurbsSurface* self, int pu, int pv, std::vector<double> ku,
              std::vector<double> kv, int nu, int nv,
              const std::vector<std::array<double, 3>>& control, std::vector<double> weights) {
             std::vector<og::Vec3> c;
             c.reserve(control.size());
             for (const auto& p : control) {
               c.push_back(vec(p));
             }
             new (self) og::NurbsSurface(pu, pv, std::move(ku), std::move(kv), nu, nv, std::move(c),
                                         std::move(weights));
           });

  m.def(
      "tfi_on_surface",
      [](const og::Surface& surface, const Points& bottom, const Points& right, const Points& top,
         const Points& left) {
        og::Boundary2D b;
        b.bottom = side(bottom);
        b.right = side(right);
        b.top = side(top);
        b.left = side(left);
        og::SurfaceTfi r;
        {
          const nb::gil_scoped_release release;
          r = og::tfi_on_surface(surface, b);
        }
        return std::make_tuple(face_to(r.grid), r.max_boundary_distance);
      },
      nb::arg("surface"), nb::arg("bottom"), nb::arg("right"), nb::arg("top"), nb::arg("left"),
      "Parametric TFI on a surface; returns (grid (nj, ni, 3), largest boundary-node distance).");

  m.def(
      "tfi_block",
      [](const Face& imin, const Face& imax, const Face& jmin, const Face& jmax, const Face& kmin,
         const Face& kmax) {
        og::BlockFaces f;
        f.imin = face_from(imin);
        f.imax = face_from(imax);
        f.jmin = face_from(jmin);
        f.jmax = face_from(jmax);
        f.kmin = face_from(kmin);
        f.kmax = face_from(kmax);
        og::Grid3D g;
        {
          const nb::gil_scoped_release release;
          g = og::tfi_block(f);
        }
        return block_to(g);
      },
      nb::arg("imin"), nb::arg("imax"), nb::arg("jmin"), nb::arg("jmax"), nb::arg("kmin"),
      nb::arg("kmax"),
      "Block TFI from six faces: imin/imax (nk, nj, 3), jmin/jmax (nk, ni, 3), kmin/kmax (nj, ni, "
      "3); returns (nk, nj, ni, 3).");

  m.def(
      "extrude",
      [](const Face& base, const std::string& mode, const std::array<double, 3>& vector,
         const std::array<double, 3>& axis_point, double angle, double distance,
         const std::vector<double>& fractions) {
        const og::Grid2D b = face_from(base);
        og::Extruded e;
        {
          const nb::gil_scoped_release release;
          if (mode == "translate") {
            e = og::extrude_translate(b, vec(vector), fractions);
          } else if (mode == "rotate") {
            e = og::extrude_rotate(b, vec(axis_point), vec(vector), angle, fractions);
          } else if (mode == "normal") {
            e = og::extrude_normal(b, distance, fractions);
          } else {
            throw std::invalid_argument("mode must be 'translate', 'rotate' or 'normal'");
          }
        }
        return std::make_tuple(block_to(e.grid), e.mirrored);
      },
      nb::arg("base"), nb::arg("mode"), nb::arg("vector") = std::array<double, 3>{0, 0, 1},
      nb::arg("axis_point") = std::array<double, 3>{0, 0, 0}, nb::arg("angle") = 0.0,
      nb::arg("distance") = 0.0, nb::arg("fractions") = std::vector<double>{0.0, 1.0},
      "Extrude a surface grid (nj, ni, 3): mode 'translate' (vector), 'rotate' (axis_point, "
      "vector = axis, angle in radians) or 'normal' (distance); returns (block, mirrored).");

  m.def(
      "smooth_elliptic_3d",
      [](const Block& grid, const std::string& control,
         const std::vector<std::tuple<int, std::optional<double>, double, double, double>>& walls,
         double omega, int max_iterations, double tolerance, og::SmoothingProgress* progress) {
        og::Elliptic3DOptions options;
        if (control == "laplace") {
          options.control = og::ControlType::Laplace;
        } else if (control == "thomas_middlecoff") {
          options.control = og::ControlType::ThomasMiddlecoff;
        } else {
          throw std::invalid_argument("control must be 'laplace' or 'thomas_middlecoff'");
        }
        for (const auto& [face, spacing, decay, relaxation, max_source] : walls) {
          og::WallControl3D w;
          w.face = face;
          w.control.enabled = true;
          w.control.spacing = spacing;
          w.control.decay = decay;
          w.control.relaxation = relaxation;
          w.control.max_source = max_source;
          options.walls.push_back(w);
        }
        options.omega = omega;
        options.max_iterations = max_iterations;
        options.tolerance = tolerance;
        const og::Grid3D g = block_from(grid);
        og::Elliptic3DResult r;
        {
          const nb::gil_scoped_release release;
          r = og::smooth_elliptic_3d(g, options, progress);
        }
        return std::make_tuple(block_to(r.grid), r.residuals, r.iterations, r.converged,
                               r.cancelled, r.warnings);
      },
      nb::arg("grid"), nb::arg("control") = "laplace",
      nb::arg("walls") =
          std::vector<std::tuple<int, std::optional<double>, double, double, double>>{},
      nb::arg("omega") = 1.0, nb::arg("max_iterations") = 500, nb::arg("tolerance") = 1e-6,
      nb::arg("progress") = nb::none(),
      "3D elliptic smoothing; walls: (face 0..5 = imin..kmax, spacing or None, decay, "
      "relaxation, max_source). Returns (grid, residuals, iterations, converged, cancelled, "
      "warnings).");

  m.def(
      "smooth_on_surface",
      [](const og::Surface& surface, const Face& grid, const std::string& control, double omega,
         int max_iterations, double tolerance, og::SmoothingProgress* progress) {
        og::EllipticOptions options;
        if (control == "laplace") {
          options.control = og::ControlType::Laplace;
        } else if (control == "thomas_middlecoff") {
          options.control = og::ControlType::ThomasMiddlecoff;
        } else {
          throw std::invalid_argument("control must be 'laplace' or 'thomas_middlecoff'");
        }
        options.omega = omega;
        options.max_iterations = max_iterations;
        options.tolerance = tolerance;
        const og::Grid2D g = face_from(grid);
        og::EllipticResult r;
        {
          const nb::gil_scoped_release release;
          r = og::smooth_on_surface(surface, g, options, progress);
        }
        return std::make_tuple(face_to(r.grid), r.residuals, r.iterations, r.converged, r.cancelled,
                               r.warnings);
      },
      nb::arg("surface"), nb::arg("grid"), nb::arg("control") = "laplace", nb::arg("omega") = 1.0,
      nb::arg("max_iterations") = 1000, nb::arg("tolerance") = 1e-6,
      nb::arg("progress") = nb::none(),
      "Elliptic smoothing of a grid (nj, ni, 3) on a surface (nodes projected back each update).");

  m.def(
      "extrude_hyperbolic",
      [](const Face& base, double distance, const std::vector<double>& fractions,
         const std::vector<std::string>& edges, double dissipation, int normal_smoothing) {
        og::HyperbolicOptions options;
        if (edges.size() != 4) {
          throw std::invalid_argument("edges: four conditions (imin, imax, jmin, jmax)");
        }
        for (std::size_t k = 0; k < 4; ++k) {
          const std::string& e = edges[k];
          if (e == "free") {
            options.edges[k] = og::HyperbolicEdge::Free;
          } else if (e == "constant_x") {
            options.edges[k] = og::HyperbolicEdge::ConstantX;
          } else if (e == "constant_y") {
            options.edges[k] = og::HyperbolicEdge::ConstantY;
          } else if (e == "constant_z") {
            options.edges[k] = og::HyperbolicEdge::ConstantZ;
          } else if (e == "periodic") {
            options.edges[k] = og::HyperbolicEdge::Periodic;
          } else {
            throw std::invalid_argument(
                "edge conditions: free, constant_x, constant_y, constant_z or periodic");
          }
        }
        options.dissipation = dissipation;
        options.normal_smoothing = normal_smoothing;
        const og::Grid2D b = face_from(base);
        og::HyperbolicResult r;
        {
          const nb::gil_scoped_release release;
          r = og::extrude_hyperbolic(b, distance, fractions, options);
        }
        return std::make_tuple(block_to(r.grid), r.mirrored, r.inverted_layers, r.warnings);
      },
      nb::arg("base"), nb::arg("distance"), nb::arg("fractions"),
      nb::arg("edges") = std::vector<std::string>{"free", "free", "free", "free"},
      nb::arg("dissipation") = 0.2, nb::arg("normal_smoothing") = 4,
      "Hyperbolic extrusion of a surface grid (nj, ni, 3); returns (block, mirrored, "
      "inverted_layers, warnings).");

  m.def(
      "quality_3d",
      [](const Block& grid) {
        const og::Grid3D g = block_from(grid);
        og::Quality3D q;
        {
          const nb::gil_scoped_release release;
          q = og::quality_3d(g);
        }
        const int ci = g.ni - 1;
        const int cj = g.nj - 1;
        const int ck = g.nk - 1;
        std::map<std::string, CellsOut> out;
        out.emplace("scaled_jacobian", cells3_to(q.scaled_jacobian, ci, cj, ck));
        out.emplace("skewness", cells3_to(q.skewness, ci, cj, ck));
        out.emplace("aspect_ratio", cells3_to(q.aspect_ratio, ci, cj, ck));
        out.emplace("min_angle", cells3_to(q.min_angle, ci, cj, ck));
        out.emplace("max_angle", cells3_to(q.max_angle, ci, cj, ck));
        out.emplace("non_orthogonality", cells3_to(q.non_orthogonality, ci, cj, ck));
        out.emplace("volume_ratio", cells3_to(q.volume_ratio, ci, cj, ck));
        out.emplace("volume", cells3_to(q.volume, ci, cj, ck));
        return out;
      },
      nb::arg("grid"), "Per-cell quality of a block (nk, nj, ni, 3); arrays (nk-1, nj-1, ni-1).");

  m.def(
      "write_cgns_3d",
      [](const std::string& path, const std::string& base_name,
         const std::vector<std::tuple<std::string, Block, std::string>>& blocks,
         const std::vector<
             std::tuple<std::string, std::vector<double>, std::vector<double>, std::vector<double>,
                        std::vector<std::tuple<std::string, std::string, std::vector<std::int64_t>,
                                               std::string>>,
                        std::string>>& zones,
         const std::vector<std::tuple<std::string, std::string>>& families,
         const std::vector<std::tuple<int, std::string, std::string, std::array<int, 6>>>& bcs,
         const std::vector<std::tuple<int, std::string, int, std::array<int, 6>, std::array<int, 6>,
                                      std::array<int, 3>>>& connections,
         const std::vector<std::tuple<int, std::string, int, std::vector<std::int64_t>,
                                      std::vector<std::int64_t>>>& vertex_connections) {
        std::vector<og::CgnsBlock> bl;
        bl.reserve(blocks.size());
        for (const auto& [name, grid, family] : blocks) {
          const og::Grid3D g = block_from(grid);
          og::CgnsBlock b{name, g.ni, g.nj, g.nk, {}, {}, {}, family};
          b.x.reserve(g.x.size());
          b.y.reserve(g.x.size());
          b.z.reserve(g.x.size());
          for (const og::Vec3& p : g.x) {
            b.x.push_back(p.x);
            b.y.push_back(p.y);
            b.z.push_back(p.z);
          }
          bl.push_back(std::move(b));
        }
        std::vector<og::CgnsUnstructuredZone3D> z;
        z.reserve(zones.size());
        for (const auto& [name, x, y, zz, sections, family] : zones) {
          og::CgnsUnstructuredZone3D zone{name, x, y, zz, {}, family};
          for (const auto& [sname, type, nodes, bc] : sections) {
            zone.sections.push_back({sname, type, nodes, bc});
          }
          z.push_back(std::move(zone));
        }
        std::vector<og::CgnsFamily> f;
        f.reserve(families.size());
        for (const auto& [name, bc_type] : families) {
          f.push_back({name, bc_type});
        }
        std::vector<og::CgnsBoundary3D> b;
        b.reserve(bcs.size());
        for (const auto& [zone, name, family, range] : bcs) {
          b.push_back({zone, name, family, range});
        }
        std::vector<og::CgnsConnection3D> c;
        c.reserve(connections.size());
        for (const auto& [zone, name, donor, range, donor_range, transform] : connections) {
          c.push_back({zone, name, donor, range, donor_range, transform});
        }
        std::vector<og::CgnsVertexConnection> vc;
        vc.reserve(vertex_connections.size());
        for (const auto& [zone, name, donor, points, donor_points] : vertex_connections) {
          vc.push_back({zone, name, donor, points, donor_points});
        }
        const nb::gil_scoped_release release;
        og::write_cgns_3d(path, base_name, bl, z, f, b, c, vc);
      },
      nb::arg("path"), nb::arg("base_name"), nb::arg("blocks"), nb::arg("zones"),
      nb::arg("families"), nb::arg("boundaries"), nb::arg("connections"),
      nb::arg("vertex_connections"),
      "Write a 3D CGNS file. Zones are numbered blocks first, then unstructured zones. blocks: "
      "(name, grid (nk, nj, ni, 3), family); zones: (name, x, y, z, sections [(name, type, "
      "0-based nodes, bc family or '')], family); families: (name, bc_type or ''); boundaries "
      "(blocks only): (zone, name, family, (imin, jmin, kmin, imax, jmax, kmax)); connections "
      "(blocks only): (zone, name, donor zone, range, donor range, transform); "
      "vertex_connections: (zone, name, donor zone, points, donor points), points as (i, j, k) "
      "triples on blocks and vertex numbers on unstructured zones. All indices 1-based except "
      "section nodes.");
  m.def(
      "inflate_3d",
      [](const nb::ndarray<const double, nb::shape<-1, 3>, nb::c_contig, nb::device::cpu>& points,
         const std::vector<std::vector<std::int64_t>>& faces, const std::vector<bool>& wall,
         const std::vector<std::int64_t>& face_group, const std::vector<std::int64_t>& node_group,
         double first_height, double growth_rate, int layers, std::optional<double> total_thickness,
         bool transition, std::optional<double> max_size, double min_jacobian, double max_scale,
         double max_slope, double proximity, double fan_threshold, double fan_step) {
        og::Inflation3DInput input;
        input.points.resize(points.shape(0));
        const double* d = points.data();
        for (std::size_t i = 0; i < input.points.size(); ++i) {
          input.points[i] = {d[3 * i], d[3 * i + 1], d[3 * i + 2]};
        }
        input.faces = faces;
        input.wall = wall;
        input.face_group = face_group;
        input.node_group = node_group;
        og::Inflation3DOptions options;
        options.first_height = first_height;
        options.growth_rate = growth_rate;
        options.max_layers = layers;
        options.total_thickness = total_thickness;
        options.transition = transition;
        options.max_size = max_size;
        options.min_jacobian = min_jacobian;
        options.max_scale = max_scale;
        options.max_slope = max_slope;
        options.proximity = proximity;
        options.fan_threshold_deg = fan_threshold;
        options.fan_step_deg = fan_step;
        og::Inflation3DResult r;
        {
          const nb::gil_scoped_release release;
          r = og::inflate_3d(input, options);
        }
        auto* xyz = new double[r.points.size() * 3];
        for (std::size_t i = 0; i < r.points.size(); ++i) {
          xyz[3 * i] = r.points[i][0];
          xyz[3 * i + 1] = r.points[i][1];
          xyz[3 * i + 2] = r.points[i][2];
        }
        const nb::capsule own(xyz, [](void* ptr) noexcept { delete[] static_cast<double*>(ptr); });
        nb::dict out;
        out["points"] = nb::ndarray<nb::numpy, double, nb::ndim<2>>(xyz, {r.points.size(), 3}, own);
        out["tets"] = index_array(r.tets);
        out["pyramids"] = index_array(r.pyramids);
        out["prisms"] = index_array(r.prisms);
        out["hexes"] = index_array(r.hexes);
        out["front"] = r.front;
        out["front_face"] = r.front_face;
        out["sides"] = r.sides;
        out["fan_columns"] = r.fan_columns;
        out["side_face"] = r.side_face;
        out["top"] = r.top;
        out["height"] = r.height;
        out["layers"] = r.layers;
        out["reason"] = r.reason;
        out["max_layers"] = r.max_layers;
        out["stops"] = r.stops;
        out["warnings"] = r.warnings;
        return out;
      },
      nb::arg("points"), nb::arg("faces"), nb::arg("wall"), nb::arg("face_group"),
      nb::arg("node_group"), nb::arg("first_height"), nb::arg("growth_rate"), nb::arg("layers"),
      nb::arg("total_thickness").none(), nb::arg("transition"), nb::arg("max_size").none(),
      nb::arg("min_jacobian") = 0.2, nb::arg("max_scale") = 2.0, nb::arg("max_slope") = 0.5,
      nb::arg("proximity") = 1.0, nb::arg("fan_threshold") = 60.0, nb::arg("fan_step") = 30.0,
      "3D boundary-layer inflation (see inflation3d.hpp). faces: closed surface oriented out of "
      "the region; wall / face_group per face; node_group per point (plane group a wall node "
      "may slide in, or -1). Returns a dict: points (input first), tets, pyramids, prisms, "
      "hexes, front / front_face, sides / side_face, top, height, layers, max_layers, stops, "
      "warnings.");

  m.def(
      "optimize_3d",
      [](const Points& points, const Cells<4>& tets, const Cells<5>& pyramids,
         const Cells<6>& prisms, const Cells<8>& hexes, const std::vector<bool>& movable,
         std::int64_t n_fixed_tets, double threshold, int passes, bool move_nodes, bool swap) {
        og::Optimize3DInput input;
        input.points.resize(points.shape(0));
        const double* d = points.data();
        for (std::size_t i = 0; i < input.points.size(); ++i) {
          input.points[i] = {d[3 * i], d[3 * i + 1], d[3 * i + 2]};
        }
        input.tets = index_input<4>(tets);
        input.pyramids = index_input<5>(pyramids);
        input.prisms = index_input<6>(prisms);
        input.hexes = index_input<8>(hexes);
        input.movable = movable;
        input.n_fixed_tets = n_fixed_tets;
        og::Optimize3DOptions options;
        options.threshold = threshold;
        options.passes = passes;
        options.move_nodes = move_nodes;
        options.swap = swap;
        og::Optimize3DResult r;
        {
          const nb::gil_scoped_release release;
          r = og::optimize_3d(input, options);
        }
        auto* xyz = new double[r.points.size() * 3];
        for (std::size_t i = 0; i < r.points.size(); ++i) {
          xyz[3 * i] = r.points[i][0];
          xyz[3 * i + 1] = r.points[i][1];
          xyz[3 * i + 2] = r.points[i][2];
        }
        const nb::capsule own(xyz, [](void* ptr) noexcept { delete[] static_cast<double*>(ptr); });
        nb::dict out;
        out["points"] = nb::ndarray<nb::numpy, double, nb::ndim<2>>(xyz, {r.points.size(), 3}, own);
        out["tets"] = index_array(r.tets);
        out["moved"] = r.moved;
        out["swaps_23"] = r.swaps_23;
        out["swaps_32"] = r.swaps_32;
        out["worst_before"] = r.worst_before;
        out["worst_after"] = r.worst_after;
        out["below_before"] = r.below_before;
        out["below_after"] = r.below_after;
        out["inverted_after"] = r.inverted_after;
        return out;
      },
      nb::arg("points"), nb::arg("tets"), nb::arg("pyramids"), nb::arg("prisms"), nb::arg("hexes"),
      nb::arg("movable"), nb::arg("n_fixed_tets") = 0, nb::arg("threshold") = 0.3,
      nb::arg("passes") = 4, nb::arg("move_nodes") = true, nb::arg("swap") = true,
      "Quality optimisation of a mixed 3D mesh (see optimize3d.hpp): node moves of movable points "
      "and 2-3 / 3-2 swaps of the tetrahedra after the first n_fixed_tets. Cell arrays are "
      "int64 (n, 4 / 5 / 6 / 8). Returns a dict: points, tets (fixed ones first), moved, "
      "swaps_23, swaps_32, worst_before, worst_after, below_before, below_after, "
      "inverted_after.");

  m.def(
      "cell_quality_3d",
      [](const Points& points, const std::vector<std::int64_t>& nodes) {
        std::vector<std::array<double, 3>> p(points.shape(0));
        const double* d = points.data();
        for (std::size_t i = 0; i < p.size(); ++i) {
          p[i] = {d[3 * i], d[3 * i + 1], d[3 * i + 2]};
        }
        for (const std::int64_t v : nodes) {
          if (v < 0 || std::cmp_greater_equal(v, p.size())) {
            throw std::invalid_argument("cell_quality_3d: node out of range");
          }
        }
        return og::cell_quality_3d(p, nodes.data(), static_cast<int>(nodes.size()));
      },
      nb::arg("points"), nb::arg("nodes"),
      "Signed scaled Jacobian of one cell (4, 5, 6 or 8 nodes; see optimize3d.hpp).");
}
