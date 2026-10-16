#include "opengrid/core/io/cgns_structured.hpp"

#include <array>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>

#ifdef OG_WITH_CGNS
#include <cgnslib.h>
#endif

namespace og {

#ifndef OG_WITH_CGNS

bool cgns_available() {
  return false;
}

void write_cgns_2d(const std::string&, const std::string&, const std::vector<CgnsZone>&,
                   const std::vector<CgnsUnstructuredZone>&, const std::vector<CgnsFamily>&,
                   const std::vector<CgnsBoundary>&, const std::vector<CgnsEdgeBoundary>&,
                   const std::vector<CgnsConnection>&, const std::vector<CgnsVertexConnection>&) {
  throw std::runtime_error("cgns: OpenGrid was built without the CGNS library");
}

void write_cgns_3d(const std::string&, const std::string&, const std::vector<CgnsBlock>&,
                   const std::vector<CgnsUnstructuredZone3D>&, const std::vector<CgnsFamily>&,
                   const std::vector<CgnsBoundary3D>&, const std::vector<CgnsConnection3D>&,
                   const std::vector<CgnsVertexConnection>&) {
  throw std::runtime_error("cgns: OpenGrid was built without the CGNS library");
}

void write_cgns_elements(const std::string&, const std::string&, int,
                         const std::vector<CgnsUnstructuredZone3D>&, const std::vector<CgnsFamily>&,
                         const std::vector<CgnsVertexConnection>&) {
  throw std::runtime_error("cgns: OpenGrid was built without the CGNS library");
}

#else

namespace {

void check(int status, const char* what) {
  if (status != CG_OK) {
    throw std::runtime_error(std::string("cgns: ") + what + ": " + cg_get_error());
  }
}

void check_name(const std::string& name, const char* what) {
  if (name.empty() || name.size() > 32) {
    throw std::invalid_argument(std::string("cgns: ") + what + " name '" + name +
                                "' must have 1 to 32 characters");
  }
}

CGNS_ENUMT(BCType_t) bc_type(const std::string& name) {
  if (name == "BCWall") {
    return CGNS_ENUMV(BCWall);
  }
  if (name == "BCInflow") {
    return CGNS_ENUMV(BCInflow);
  }
  if (name == "BCOutflow") {
    return CGNS_ENUMV(BCOutflow);
  }
  if (name == "BCSymmetryPlane") {
    return CGNS_ENUMV(BCSymmetryPlane);
  }
  if (name == "BCFarfield") {
    return CGNS_ENUMV(BCFarfield);
  }
  if (name == "BCTypeUserDefined") {
    return CGNS_ENUMV(BCTypeUserDefined);
  }
  throw std::invalid_argument("cgns: unknown BC type '" + name + "'");
}

// Element types by SIDS name: linear and the complete Lagrange elements up to order 4.
CGNS_ENUMT(ElementType_t) element_type(const std::string& name, int& nodes) {
  static const std::array<std::pair<const char*, CGNS_ENUMT(ElementType_t)>, 28> types = {{
      {"BAR_2", CGNS_ENUMV(BAR_2)},       {"BAR_3", CGNS_ENUMV(BAR_3)},
      {"BAR_4", CGNS_ENUMV(BAR_4)},       {"BAR_5", CGNS_ENUMV(BAR_5)},
      {"TRI_3", CGNS_ENUMV(TRI_3)},       {"TRI_6", CGNS_ENUMV(TRI_6)},
      {"TRI_10", CGNS_ENUMV(TRI_10)},     {"TRI_15", CGNS_ENUMV(TRI_15)},
      {"QUAD_4", CGNS_ENUMV(QUAD_4)},     {"QUAD_9", CGNS_ENUMV(QUAD_9)},
      {"QUAD_16", CGNS_ENUMV(QUAD_16)},   {"QUAD_25", CGNS_ENUMV(QUAD_25)},
      {"TETRA_4", CGNS_ENUMV(TETRA_4)},   {"TETRA_10", CGNS_ENUMV(TETRA_10)},
      {"TETRA_20", CGNS_ENUMV(TETRA_20)}, {"TETRA_35", CGNS_ENUMV(TETRA_35)},
      {"PYRA_5", CGNS_ENUMV(PYRA_5)},     {"PYRA_14", CGNS_ENUMV(PYRA_14)},
      {"PYRA_30", CGNS_ENUMV(PYRA_30)},   {"PYRA_55", CGNS_ENUMV(PYRA_55)},
      {"PENTA_6", CGNS_ENUMV(PENTA_6)},   {"PENTA_18", CGNS_ENUMV(PENTA_18)},
      {"PENTA_40", CGNS_ENUMV(PENTA_40)}, {"PENTA_75", CGNS_ENUMV(PENTA_75)},
      {"HEXA_8", CGNS_ENUMV(HEXA_8)},     {"HEXA_27", CGNS_ENUMV(HEXA_27)},
      {"HEXA_64", CGNS_ENUMV(HEXA_64)},   {"HEXA_125", CGNS_ENUMV(HEXA_125)},
  }};
  for (const auto& [label, type] : types) {
    if (name == label) {
      check(cg_npe(type, &nodes), "nodes per element");
      return type;
    }
  }
  throw std::invalid_argument("cgns: unknown element type '" + name + "'");
}

// RAII guard closing the file on every path.
struct File {
  int handle = -1;
  ~File() {
    if (handle >= 0) {
      cg_close(handle);
    }
  }
};

}  // namespace

bool cgns_available() {
  return true;
}

void write_cgns_2d(const std::string& path, const std::string& base_name,
                   const std::vector<CgnsZone>& zones,
                   const std::vector<CgnsUnstructuredZone>& unstructured_zones,
                   const std::vector<CgnsFamily>& families,
                   const std::vector<CgnsBoundary>& boundaries,
                   const std::vector<CgnsEdgeBoundary>& edge_boundaries,
                   const std::vector<CgnsConnection>& connections,
                   const std::vector<CgnsVertexConnection>& vertex_connections) {
  check_name(base_name, "base");
  for (const CgnsUnstructuredZone& zone : unstructured_zones) {
    check_name(zone.name, "zone");
    const std::size_t n = zone.x.size();
    bool ok = zone.y.size() == n && n >= 3 && !(zone.triangles.empty() && zone.quads.empty());
    for (const auto& t : zone.triangles) {
      for (const std::int64_t v : t) {
        ok = ok && v >= 0 && std::cmp_less(v, n);
      }
    }
    for (const auto& q : zone.quads) {
      for (const std::int64_t v : q) {
        ok = ok && v >= 0 && std::cmp_less(v, n);
      }
    }
    if (!ok) {
      throw std::invalid_argument("cgns: unstructured zone '" + zone.name +
                                  "' has inconsistent sizes or vertex ids");
    }
  }
  std::vector<std::string> zone_names;
  zone_names.reserve(zones.size() + unstructured_zones.size());
  for (const CgnsZone& zone : zones) {
    zone_names.push_back(zone.name);
  }
  for (const CgnsUnstructuredZone& zone : unstructured_zones) {
    zone_names.push_back(zone.name);
  }
  for (const CgnsZone& zone : zones) {
    check_name(zone.name, "zone");
    const auto n = static_cast<std::size_t>(zone.ni) * static_cast<std::size_t>(zone.nj);
    if (zone.ni < 2 || zone.nj < 2 || zone.x.size() != n || zone.y.size() != n) {
      throw std::invalid_argument("cgns: zone '" + zone.name + "' has inconsistent sizes");
    }
  }
  for (const CgnsFamily& family : families) {
    check_name(family.name, "family");
    if (!family.bc_type.empty()) {
      bc_type(family.bc_type);
    }
  }

  File file;
  check(cg_open(path.c_str(), CG_MODE_WRITE, &file.handle), "open");
  int base = 0;
  check(cg_base_write(file.handle, base_name.c_str(), 2, 2, &base), "base");
  // Project coordinates are in metres (".ogp" units); declared once at base level.
  check(cg_goto(file.handle, base, "end"), "goto base");
  check(cg_dataclass_write(CGNS_ENUMV(Dimensional)), "data class");
  check(cg_units_write(CGNS_ENUMV(Kilogram), CGNS_ENUMV(Meter), CGNS_ENUMV(Second),
                       CGNS_ENUMV(Kelvin), CGNS_ENUMV(Radian)),
        "units");
  for (const CgnsFamily& family : families) {
    int fam = 0;
    check(cg_family_write(file.handle, base, family.name.c_str(), &fam), "family");
    if (!family.bc_type.empty()) {
      int fambc = 0;
      check(cg_fambc_write(file.handle, base, fam, "FamBC", bc_type(family.bc_type), &fambc),
            "family BC");
    }
  }
  std::vector<int> zone_ids;
  auto write_coordinates = [&](int z, const std::vector<double>& x, const std::vector<double>& y,
                               const std::string& family) {
    for (const auto& [label, values] :
         {std::pair{"CoordinateX", &x}, std::pair{"CoordinateY", &y}}) {
      int coord = 0;
      check(cg_coord_write(file.handle, base, z, CGNS_ENUMV(RealDouble), label, values->data(),
                           &coord),
            label);
      // Lengths: exponents (mass, length, time, temperature, angle) = (0, 1, 0, 0, 0).
      const std::array<double, 5> exponents = {0.0, 1.0, 0.0, 0.0, 0.0};
      check(cg_goto(file.handle, base, "Zone_t", z, "GridCoordinates_t", 1, "DataArray_t", coord,
                    "end"),
            "goto coordinate");
      check(cg_exponents_write(CGNS_ENUMV(RealDouble), exponents.data()), "exponents");
    }
    if (!family.empty()) {
      check(cg_goto(file.handle, base, "Zone_t", z, "end"), "goto zone");
      check(cg_famname_write(family.c_str()), "zone family");
    }
  };
  for (const CgnsZone& zone : zones) {
    std::array<cgsize_t, 6> size = {zone.ni, zone.nj, zone.ni - 1, zone.nj - 1, 0, 0};
    int z = 0;
    check(cg_zone_write(file.handle, base, zone.name.c_str(), size.data(), CGNS_ENUMV(Structured),
                        &z),
          "zone");
    write_coordinates(z, zone.x, zone.y, zone.family);
    zone_ids.push_back(z);
  }
  // Unstructured zones: cells, then the bar sections of their boundary conditions.
  for (std::size_t uz = 0; uz < unstructured_zones.size(); ++uz) {
    const CgnsUnstructuredZone& zone = unstructured_zones[uz];
    const auto n_cells =
        static_cast<cgsize_t>(zone.triangles.size()) + static_cast<cgsize_t>(zone.quads.size());
    std::array<cgsize_t, 3> size = {static_cast<cgsize_t>(zone.x.size()), n_cells, 0};
    int z = 0;
    check(cg_zone_write(file.handle, base, zone.name.c_str(), size.data(), CGNS_ENUMV(Unstructured),
                        &z),
          "zone");
    write_coordinates(z, zone.x, zone.y, zone.family);
    zone_ids.push_back(z);
    cgsize_t next = 1;
    int section = 0;
    if (!zone.triangles.empty()) {
      std::vector<cgsize_t> data;
      data.reserve(3 * zone.triangles.size());
      for (const auto& t : zone.triangles) {
        for (const std::int64_t v : t) {
          data.push_back(static_cast<cgsize_t>(v + 1));
        }
      }
      const auto last = next + static_cast<cgsize_t>(zone.triangles.size()) - 1;
      check(cg_section_write(file.handle, base, z, "Triangles", CGNS_ENUMV(TRI_3), next, last, 0,
                             data.data(), &section),
            "triangle section");
      next = last + 1;
    }
    if (!zone.quads.empty()) {
      std::vector<cgsize_t> data;
      data.reserve(4 * zone.quads.size());
      for (const auto& q : zone.quads) {
        for (const std::int64_t v : q) {
          data.push_back(static_cast<cgsize_t>(v + 1));
        }
      }
      const auto last = next + static_cast<cgsize_t>(zone.quads.size()) - 1;
      check(cg_section_write(file.handle, base, z, "Quadrilaterals", CGNS_ENUMV(QUAD_4), next, last,
                             0, data.data(), &section),
            "quad section");
      next = last + 1;
    }
    for (const CgnsEdgeBoundary& bc : edge_boundaries) {
      if (!std::cmp_equal(bc.zone, uz)) {
        continue;
      }
      check_name(bc.name, "BC");
      if (bc.edges.empty()) {
        throw std::invalid_argument("cgns: boundary '" + bc.name + "' has no edges");
      }
      std::vector<cgsize_t> data;
      data.reserve(2 * bc.edges.size());
      for (const auto& e : bc.edges) {
        for (const std::int64_t v : e) {
          if (v < 0 || std::cmp_greater_equal(v, zone.x.size())) {
            throw std::invalid_argument("cgns: boundary '" + bc.name +
                                        "' refers to a missing vertex");
          }
          data.push_back(static_cast<cgsize_t>(v + 1));
        }
      }
      const auto last = next + static_cast<cgsize_t>(bc.edges.size()) - 1;
      check(cg_section_write(file.handle, base, z, bc.name.c_str(), CGNS_ENUMV(BAR_2), next, last,
                             0, data.data(), &section),
            "boundary section");
      std::array<cgsize_t, 2> range = {next, last};
      int b = 0;
      check(cg_boco_write(file.handle, base, z, bc.name.c_str(), CGNS_ENUMV(FamilySpecified),
                          CGNS_ENUMV(PointRange), 2, range.data(), &b),
            "BC");
      check(cg_boco_gridlocation_write(file.handle, base, z, b, CGNS_ENUMV(EdgeCenter)),
            "BC location");
      check(cg_goto(file.handle, base, "Zone_t", z, "ZoneBC_t", 1, "BC_t", b, "end"), "goto BC");
      check(cg_famname_write(bc.family.c_str()), "BC family");
      next = last + 1;
    }
  }
  for (const CgnsBoundary& bc : boundaries) {
    check_name(bc.name, "BC");
    const int z = zone_ids.at(static_cast<std::size_t>(bc.zone));
    std::array<cgsize_t, 4> range = {bc.range[0], bc.range[1], bc.range[2], bc.range[3]};
    int b = 0;
    check(cg_boco_write(file.handle, base, z, bc.name.c_str(), CGNS_ENUMV(FamilySpecified),
                        CGNS_ENUMV(PointRange), 2, range.data(), &b),
          "BC");
    check(cg_goto(file.handle, base, "Zone_t", z, "ZoneBC_t", 1, "BC_t", b, "end"), "goto BC");
    check(cg_famname_write(bc.family.c_str()), "BC family");
  }
  for (const CgnsConnection& c : connections) {
    check_name(c.name, "connection");
    const int z = zone_ids.at(static_cast<std::size_t>(c.zone));
    const std::string& donor = zone_names.at(static_cast<std::size_t>(c.donor_zone));
    std::array<cgsize_t, 4> range = {c.range[0], c.range[1], c.range[2], c.range[3]};
    std::array<cgsize_t, 4> donor_range = {c.donor_range[0], c.donor_range[1], c.donor_range[2],
                                           c.donor_range[3]};
    std::array<int, 2> transform = {c.transform[0], c.transform[1]};
    int index = 0;
    check(cg_1to1_write(file.handle, base, z, c.name.c_str(), donor.c_str(), range.data(),
                        donor_range.data(), transform.data(), &index),
          "1to1");
    if (c.periodic) {
      // The library stores these as single precision (RealSingle) arrays.
      const std::array<float, 2> center = {0.0F, 0.0F};
      const std::array<float, 2> angle = {0.0F, 0.0F};
      const std::array<float, 2> translation = {static_cast<float>(c.translation[0]),
                                                static_cast<float>(c.translation[1])};
      check(cg_1to1_periodic_write(file.handle, base, z, index, center.data(), angle.data(),
                                   translation.data()),
            "periodic 1to1");
      // Dimensional exponents of the periodic data (mass, length, time, temperature, angle).
      const std::array<double, 5> length = {0.0, 1.0, 0.0, 0.0, 0.0};
      const std::array<double, 5> radian = {0.0, 0.0, 0.0, 0.0, 1.0};
      const std::string periodic_path =
          "/" + base_name + "/" + zones.at(static_cast<std::size_t>(c.zone)).name +
          "/ZoneGridConnectivity/" + c.name + "/GridConnectivityProperty/Periodic/";
      for (const auto& [array, exponents] :
           {std::pair{"RotationCenter", &length}, std::pair{"RotationAngle", &radian},
            std::pair{"Translation", &length}}) {
        check(cg_gopath(file.handle, (periodic_path + array).c_str()), "goto periodic data");
        check(cg_exponents_write(CGNS_ENUMV(RealDouble), exponents->data()), "periodic exponents");
      }
    }
  }
  for (const CgnsVertexConnection& c : vertex_connections) {
    check_name(c.name, "connection");
    const auto zi = static_cast<std::size_t>(c.zone);
    const auto di = static_cast<std::size_t>(c.donor_zone);
    const int z = zone_ids.at(zi);
    const bool donor_structured = di < zones.size();
    std::size_t dim = 1;
    if (zi < zones.size()) {
      dim = 2;
    }
    std::size_t donor_dim = 1;
    if (donor_structured) {
      donor_dim = 2;
    }
    if (c.points.empty() || c.points.size() % dim != 0 ||
        c.points.size() / dim != c.donor_points.size() / donor_dim ||
        c.donor_points.size() % donor_dim != 0) {
      throw std::invalid_argument("cgns: connection '" + c.name + "' has inconsistent point lists");
    }
    std::vector<cgsize_t> points(c.points.begin(), c.points.end());
    std::vector<cgsize_t> donor(c.donor_points.begin(), c.donor_points.end());
    CGNS_ENUMT(ZoneType_t) donor_type = CGNS_ENUMV(Unstructured);
    if (donor_structured) {
      donor_type = CGNS_ENUMV(Structured);
    }
    // The donor data are cgsize_t (32 or 64 bit depending on the library build).
    CGNS_ENUMT(DataType_t) donor_data_type = CGNS_ENUMV(Integer);
    if (sizeof(cgsize_t) == 8) {
      donor_data_type = CGNS_ENUMV(LongInteger);
    }
    int index = 0;
    check(cg_conn_write(
              file.handle, base, z, c.name.c_str(), CGNS_ENUMV(Vertex), CGNS_ENUMV(Abutting1to1),
              CGNS_ENUMV(PointList), static_cast<cgsize_t>(points.size() / dim), points.data(),
              zone_names.at(di).c_str(), donor_type, CGNS_ENUMV(PointListDonor), donor_data_type,
              static_cast<cgsize_t>(donor.size() / donor_dim), donor.data(), &index),
          "connectivity");
    if (c.periodic) {
      const std::array<float, 2> center = {0.0F, 0.0F};
      const std::array<float, 2> angle = {0.0F, 0.0F};
      const std::array<float, 2> translation = {static_cast<float>(c.translation[0]),
                                                static_cast<float>(c.translation[1])};
      check(cg_conn_periodic_write(file.handle, base, z, index, center.data(), angle.data(),
                                   translation.data()),
            "periodic connectivity");
      const std::array<double, 5> length = {0.0, 1.0, 0.0, 0.0, 0.0};
      const std::array<double, 5> radian = {0.0, 0.0, 0.0, 0.0, 1.0};
      const std::string periodic_path = "/" + base_name + "/" + zone_names.at(zi) +
                                        "/ZoneGridConnectivity/" + c.name +
                                        "/GridConnectivityProperty/Periodic/";
      for (const auto& [array, exponents] :
           {std::pair{"RotationCenter", &length}, std::pair{"RotationAngle", &radian},
            std::pair{"Translation", &length}}) {
        check(cg_gopath(file.handle, (periodic_path + array).c_str()), "goto periodic data");
        check(cg_exponents_write(CGNS_ENUMV(RealDouble), exponents->data()), "periodic exponents");
      }
    }
  }
  const int handle = file.handle;
  file.handle = -1;
  check(cg_close(handle), "close");
}

namespace {

void write_cgns_impl(const std::string& path, const std::string& base_name, int cell_dimension,
                     const std::vector<CgnsBlock>& blocks,
                     const std::vector<CgnsUnstructuredZone3D>& unstructured_zones,
                     const std::vector<CgnsFamily>& families,
                     const std::vector<CgnsBoundary3D>& boundaries,
                     const std::vector<CgnsConnection3D>& connections,
                     const std::vector<CgnsVertexConnection>& vertex_connections) {
  check_name(base_name, "base");
  if (cell_dimension != 2 && cell_dimension != 3) {
    throw std::invalid_argument("cgns: the cell dimension must be 2 or 3");
  }
  if (cell_dimension == 2 && !blocks.empty()) {
    throw std::invalid_argument("cgns: structured blocks need a 3D base");
  }
  std::vector<std::string> zone_names;
  std::vector<std::size_t> zone_vertices;
  for (const CgnsBlock& b : blocks) {
    check_name(b.name, "zone");
    const auto n = static_cast<std::size_t>(b.ni) * static_cast<std::size_t>(b.nj) *
                   static_cast<std::size_t>(b.nk);
    if (b.ni < 2 || b.nj < 2 || b.nk < 2 || b.x.size() != n || b.y.size() != n || b.z.size() != n) {
      throw std::invalid_argument("cgns: block '" + b.name + "' has inconsistent sizes");
    }
    zone_names.push_back(b.name);
    zone_vertices.push_back(n);
  }
  for (const CgnsUnstructuredZone3D& zone : unstructured_zones) {
    check_name(zone.name, "zone");
    const std::size_t nv = zone.x.size();
    if (zone.y.size() != nv || zone.z.size() != nv) {
      throw std::invalid_argument("cgns: zone '" + zone.name + "' has inconsistent sizes");
    }
    for (const CgnsSection& s : zone.sections) {
      int nodes = 0;
      (void)element_type(s.type, nodes);
      if (s.nodes.size() % static_cast<std::size_t>(nodes) != 0) {
        throw std::invalid_argument("cgns: section '" + s.name + "' has a partial element");
      }
      for (const std::int64_t v : s.nodes) {
        if (v < 0 || std::cmp_greater_equal(v, nv)) {
          throw std::invalid_argument("cgns: section '" + s.name + "' refers to a missing vertex");
        }
      }
    }
    zone_names.push_back(zone.name);
    zone_vertices.push_back(nv);
  }
  for (const CgnsFamily& family : families) {
    check_name(family.name, "family");
    if (!family.bc_type.empty()) {
      bc_type(family.bc_type);
    }
  }
  File file;
  check(cg_open(path.c_str(), CG_MODE_WRITE, &file.handle), "open");
  int base = 0;
  check(cg_base_write(file.handle, base_name.c_str(), cell_dimension, cell_dimension, &base),
        "base");
  check(cg_goto(file.handle, base, "end"), "goto base");
  check(cg_dataclass_write(CGNS_ENUMV(Dimensional)), "data class");
  check(cg_units_write(CGNS_ENUMV(Kilogram), CGNS_ENUMV(Meter), CGNS_ENUMV(Second),
                       CGNS_ENUMV(Kelvin), CGNS_ENUMV(Radian)),
        "units");
  for (const CgnsFamily& family : families) {
    int fam = 0;
    check(cg_family_write(file.handle, base, family.name.c_str(), &fam), "family");
    if (!family.bc_type.empty()) {
      int fambc = 0;
      check(cg_fambc_write(file.handle, base, fam, "FamBC", bc_type(family.bc_type), &fambc),
            "family BC");
    }
  }
  std::vector<int> zone_ids;
  auto write_coordinates = [&](int z, const std::vector<double>& x, const std::vector<double>& y,
                               const std::vector<double>& zc, const std::string& family) {
    const std::array<std::pair<const char*, const std::vector<double>*>, 3> axes = {
        {{"CoordinateX", &x}, {"CoordinateY", &y}, {"CoordinateZ", &zc}}};
    for (int axis = 0; axis < cell_dimension; ++axis) {
      const auto& [label, values] = axes.at(static_cast<std::size_t>(axis));
      int coord = 0;
      check(cg_coord_write(file.handle, base, z, CGNS_ENUMV(RealDouble), label, values->data(),
                           &coord),
            label);
      const std::array<double, 5> exponents = {0.0, 1.0, 0.0, 0.0, 0.0};
      check(cg_goto(file.handle, base, "Zone_t", z, "GridCoordinates_t", 1, "DataArray_t", coord,
                    "end"),
            "goto coordinate");
      check(cg_exponents_write(CGNS_ENUMV(RealDouble), exponents.data()), "exponents");
    }
    if (!family.empty()) {
      check(cg_goto(file.handle, base, "Zone_t", z, "end"), "goto zone");
      check(cg_famname_write(family.c_str()), "zone family");
    }
  };
  for (const CgnsBlock& b : blocks) {
    std::array<cgsize_t, 9> size = {b.ni, b.nj, b.nk, b.ni - 1, b.nj - 1, b.nk - 1, 0, 0, 0};
    int z = 0;
    check(cg_zone_write(file.handle, base, b.name.c_str(), size.data(), CGNS_ENUMV(Structured), &z),
          "zone");
    write_coordinates(z, b.x, b.y, b.z, b.family);
    zone_ids.push_back(z);
  }
  for (const CgnsUnstructuredZone3D& zone : unstructured_zones) {
    cgsize_t cells = 0;
    for (const CgnsSection& s : zone.sections) {
      int nodes = 0;
      (void)element_type(s.type, nodes);
      if (s.bc_family.empty()) {
        cells += static_cast<cgsize_t>(s.nodes.size()) / nodes;
      }
    }
    std::array<cgsize_t, 3> size = {static_cast<cgsize_t>(zone.x.size()), cells, 0};
    int z = 0;
    check(cg_zone_write(file.handle, base, zone.name.c_str(), size.data(), CGNS_ENUMV(Unstructured),
                        &z),
          "zone");
    write_coordinates(z, zone.x, zone.y, zone.z, zone.family);
    // Volume sections first, then boundary sections (CGNS numbers volume cells first).
    cgsize_t next = 1;
    for (const bool boundary : {false, true}) {
      for (const CgnsSection& s : zone.sections) {
        if (s.bc_family.empty() == boundary) {
          continue;
        }
        check_name(s.name, "section");
        int nodes = 0;
        const auto type = element_type(s.type, nodes);
        const auto count = static_cast<cgsize_t>(s.nodes.size()) / nodes;
        if (count == 0) {
          continue;
        }
        std::vector<cgsize_t> data;
        data.reserve(s.nodes.size());
        for (const std::int64_t v : s.nodes) {
          data.push_back(static_cast<cgsize_t>(v + 1));
        }
        const cgsize_t last = next + count - 1;
        int section = 0;
        check(cg_section_write(file.handle, base, z, s.name.c_str(), type, next, last, 0,
                               data.data(), &section),
              "section");
        if (boundary) {
          std::array<cgsize_t, 2> range = {next, last};
          int b = 0;
          check(cg_boco_write(file.handle, base, z, s.name.c_str(), CGNS_ENUMV(FamilySpecified),
                              CGNS_ENUMV(PointRange), 2, range.data(), &b),
                "BC");
          CGNS_ENUMT(GridLocation_t) location = CGNS_ENUMV(FaceCenter);
          if (cell_dimension == 2) {
            location = CGNS_ENUMV(EdgeCenter);
          }
          check(cg_boco_gridlocation_write(file.handle, base, z, b, location), "BC location");
          check(cg_goto(file.handle, base, "Zone_t", z, "ZoneBC_t", 1, "BC_t", b, "end"),
                "goto BC");
          check(cg_famname_write(s.bc_family.c_str()), "BC family");
        }
        next = last + 1;
      }
    }
    zone_ids.push_back(z);
  }
  for (const CgnsBoundary3D& bc : boundaries) {
    check_name(bc.name, "BC");
    if (bc.zone < 0 || std::cmp_greater_equal(bc.zone, blocks.size())) {
      throw std::invalid_argument("cgns: BC '" + bc.name + "' is not on a structured zone");
    }
    const int z = zone_ids.at(static_cast<std::size_t>(bc.zone));
    std::array<cgsize_t, 6> range{};
    for (std::size_t k = 0; k < 6; ++k) {
      range[k] = bc.range[k];
    }
    int b = 0;
    check(cg_boco_write(file.handle, base, z, bc.name.c_str(), CGNS_ENUMV(FamilySpecified),
                        CGNS_ENUMV(PointRange), 2, range.data(), &b),
          "BC");
    check(cg_goto(file.handle, base, "Zone_t", z, "ZoneBC_t", 1, "BC_t", b, "end"), "goto BC");
    check(cg_famname_write(bc.family.c_str()), "BC family");
  }
  for (const CgnsConnection3D& c : connections) {
    check_name(c.name, "connection");
    if (c.zone < 0 || c.donor_zone < 0 || std::cmp_greater_equal(c.zone, blocks.size()) ||
        std::cmp_greater_equal(c.donor_zone, blocks.size())) {
      throw std::invalid_argument("cgns: 1-to-1 connection '" + c.name +
                                  "' is not between structured zones");
    }
    const int z = zone_ids.at(static_cast<std::size_t>(c.zone));
    const std::string& donor = blocks.at(static_cast<std::size_t>(c.donor_zone)).name;
    std::array<cgsize_t, 6> range{};
    std::array<cgsize_t, 6> donor_range{};
    for (std::size_t k = 0; k < 6; ++k) {
      range[k] = c.range[k];
      donor_range[k] = c.donor_range[k];
    }
    std::array<int, 3> transform = c.transform;
    int index = 0;
    check(cg_1to1_write(file.handle, base, z, c.name.c_str(), donor.c_str(), range.data(),
                        donor_range.data(), transform.data(), &index),
          "1to1");
  }
  // The donor data are cgsize_t (32 or 64 bit depending on the library build).
  CGNS_ENUMT(DataType_t) donor_data_type = CGNS_ENUMV(Integer);
  if (sizeof(cgsize_t) == 8) {
    donor_data_type = CGNS_ENUMV(LongInteger);
  }
  for (const CgnsVertexConnection& c : vertex_connections) {
    check_name(c.name, "connection");
    if (c.zone < 0 || c.donor_zone < 0 || std::cmp_greater_equal(c.zone, zone_ids.size()) ||
        std::cmp_greater_equal(c.donor_zone, zone_ids.size())) {
      throw std::invalid_argument("cgns: connection '" + c.name + "' refers to a missing zone");
    }
    const auto zi = static_cast<std::size_t>(c.zone);
    const auto di = static_cast<std::size_t>(c.donor_zone);
    // Structured zones: (i, j, k) triples; unstructured zones: vertex numbers.
    std::size_t dim = 1;
    if (zi < blocks.size()) {
      dim = 3;
    }
    std::size_t donor_dim = 1;
    CGNS_ENUMT(ZoneType_t) donor_type = CGNS_ENUMV(Unstructured);
    if (di < blocks.size()) {
      donor_dim = 3;
      donor_type = CGNS_ENUMV(Structured);
    }
    if (c.points.empty() || c.points.size() % dim != 0 || c.donor_points.size() % donor_dim != 0 ||
        c.points.size() / dim != c.donor_points.size() / donor_dim) {
      throw std::invalid_argument("cgns: connection '" + c.name + "' has inconsistent points");
    }
    if (dim == 1) {
      for (const std::int64_t v : c.points) {
        if (v < 1 || std::cmp_greater(v, zone_vertices[zi])) {
          throw std::invalid_argument("cgns: connection '" + c.name +
                                      "' refers to a missing vertex");
        }
      }
    }
    std::vector<cgsize_t> points(c.points.begin(), c.points.end());
    std::vector<cgsize_t> donor(c.donor_points.begin(), c.donor_points.end());
    int index = 0;
    check(cg_conn_write(file.handle, base, zone_ids[zi], c.name.c_str(), CGNS_ENUMV(Vertex),
                        CGNS_ENUMV(Abutting1to1), CGNS_ENUMV(PointList),
                        static_cast<cgsize_t>(points.size() / dim), points.data(),
                        zone_names[di].c_str(), donor_type, CGNS_ENUMV(PointListDonor),
                        donor_data_type, static_cast<cgsize_t>(donor.size() / donor_dim),
                        donor.data(), &index),
          "connectivity");
    if (c.periodic) {
      // Translational periodicity (in the base's physical dimension).
      const auto n = static_cast<std::size_t>(cell_dimension);
      const std::vector<float> center(n, 0.0F);
      const std::vector<float> angle(n, 0.0F);
      std::vector<float> translation(n, 0.0F);
      translation[0] = static_cast<float>(c.translation[0]);
      translation[1] = static_cast<float>(c.translation[1]);
      check(cg_conn_periodic_write(file.handle, base, zone_ids[zi], index, center.data(),
                                   angle.data(), translation.data()),
            "periodic connectivity");
      const std::array<double, 5> length = {0.0, 1.0, 0.0, 0.0, 0.0};
      const std::array<double, 5> radian = {0.0, 0.0, 0.0, 0.0, 1.0};
      const std::string periodic_path = "/" + base_name + "/" + zone_names.at(zi) +
                                        "/ZoneGridConnectivity/" + c.name +
                                        "/GridConnectivityProperty/Periodic/";
      for (const auto& [array, exponents] :
           {std::pair{"RotationCenter", &length}, std::pair{"RotationAngle", &radian},
            std::pair{"Translation", &length}}) {
        check(cg_gopath(file.handle, (periodic_path + array).c_str()), "goto periodic data");
        check(cg_exponents_write(CGNS_ENUMV(RealDouble), exponents->data()), "periodic exponents");
      }
    }
  }
  const int handle = file.handle;
  file.handle = -1;
  check(cg_close(handle), "close");
}

}  // namespace

void write_cgns_3d(const std::string& path, const std::string& base_name,
                   const std::vector<CgnsBlock>& blocks,
                   const std::vector<CgnsUnstructuredZone3D>& unstructured_zones,
                   const std::vector<CgnsFamily>& families,
                   const std::vector<CgnsBoundary3D>& boundaries,
                   const std::vector<CgnsConnection3D>& connections,
                   const std::vector<CgnsVertexConnection>& vertex_connections) {
  write_cgns_impl(path, base_name, 3, blocks, unstructured_zones, families, boundaries, connections,
                  vertex_connections);
}

void write_cgns_elements(const std::string& path, const std::string& base_name, int cell_dimension,
                         const std::vector<CgnsUnstructuredZone3D>& zones,
                         const std::vector<CgnsFamily>& families,
                         const std::vector<CgnsVertexConnection>& vertex_connections) {
  write_cgns_impl(path, base_name, cell_dimension, {}, zones, families, {}, {}, vertex_connections);
}

#endif

void write_cgns_structured_2d(const std::string& path, const std::string& base_name,
                              const std::vector<CgnsZone>& zones,
                              const std::vector<CgnsFamily>& families,
                              const std::vector<CgnsBoundary>& boundaries,
                              const std::vector<CgnsConnection>& connections) {
  write_cgns_2d(path, base_name, zones, {}, families, boundaries, {}, connections, {});
}

void write_cgns_structured_3d(const std::string& path, const std::string& base_name,
                              const std::vector<CgnsBlock>& blocks,
                              const std::vector<CgnsFamily>& families,
                              const std::vector<CgnsBoundary3D>& boundaries,
                              const std::vector<CgnsConnection3D>& connections) {
  write_cgns_3d(path, base_name, blocks, {}, families, boundaries, connections, {});
}

void write_cgns_unstructured_3d(const std::string& path, const std::string& base_name,
                                const std::vector<CgnsUnstructuredZone3D>& zones,
                                const std::vector<CgnsFamily>& families,
                                const std::vector<CgnsVertexConnection>& vertex_connections) {
  write_cgns_3d(path, base_name, {}, zones, families, {}, {}, vertex_connections);
}

}  // namespace og
