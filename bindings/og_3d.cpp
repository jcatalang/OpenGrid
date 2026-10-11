// Bindings of the 3D structured core: surfaces, surface TFI, block TFI, extrusion, 3D quality
// and the 3D CGNS writer. Grids cross the boundary as NumPy arrays of shape (nk, nj, ni, 3)
// (blocks) and (nj, ni, 3) (faces / surface grids).
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <array>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "opengrid/core/geom/curve.hpp"
#include "opengrid/core/geom/surface.hpp"
#include "opengrid/core/io/cgns_structured.hpp"
#include "opengrid/core/quality/quality3d.hpp"
#include "opengrid/core/structured/extrude.hpp"
#include "opengrid/core/structured/surface_tfi.hpp"
#include "opengrid/core/structured/tfi3d.hpp"

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
      "write_cgns_structured_3d",
      [](const std::string& path, const std::string& base_name,
         const std::vector<std::tuple<std::string, Block, std::string>>& blocks,
         const std::vector<std::tuple<std::string, std::string>>& families,
         const std::vector<std::tuple<int, std::string, std::string, std::array<int, 6>>>& bcs,
         const std::vector<std::tuple<int, std::string, int, std::array<int, 6>, std::array<int, 6>,
                                      std::array<int, 3>>>& connections) {
        std::vector<og::CgnsBlock> z;
        z.reserve(blocks.size());
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
          z.push_back(std::move(b));
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
        const nb::gil_scoped_release release;
        og::write_cgns_structured_3d(path, base_name, z, f, b, c);
      },
      nb::arg("path"), nb::arg("base_name"), nb::arg("blocks"), nb::arg("families"),
      nb::arg("boundaries"), nb::arg("connections"),
      "Write a 3D structured CGNS file. blocks: (name, grid (nk, nj, ni, 3), family); families: "
      "(name, bc_type or ''); boundaries: (zone, name, family, (imin, jmin, kmin, imax, jmax, "
      "kmax)); connections: (zone, name, donor zone, range, donor range, transform). 1-based.");
}
