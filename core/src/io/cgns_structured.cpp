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

void write_cgns_structured_3d(const std::string&, const std::string&, const std::vector<CgnsBlock>&,
                              const std::vector<CgnsFamily>&, const std::vector<CgnsBoundary3D>&,
                              const std::vector<CgnsConnection3D>&) {
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

void write_cgns_structured_3d(const std::string& path, const std::string& base_name,
                              const std::vector<CgnsBlock>& blocks,
                              const std::vector<CgnsFamily>& families,
                              const std::vector<CgnsBoundary3D>& boundaries,
                              const std::vector<CgnsConnection3D>& connections) {
  check_name(base_name, "base");
  for (const CgnsBlock& b : blocks) {
    check_name(b.name, "zone");
    const auto n = static_cast<std::size_t>(b.ni) * static_cast<std::size_t>(b.nj) *
                   static_cast<std::size_t>(b.nk);
    if (b.ni < 2 || b.nj < 2 || b.nk < 2 || b.x.size() != n || b.y.size() != n || b.z.size() != n) {
      throw std::invalid_argument("cgns: block '" + b.name + "' has inconsistent sizes");
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
  check(cg_base_write(file.handle, base_name.c_str(), 3, 3, &base), "base");
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
  for (const CgnsBlock& b : blocks) {
    std::array<cgsize_t, 9> size = {b.ni, b.nj, b.nk, b.ni - 1, b.nj - 1, b.nk - 1, 0, 0, 0};
    int z = 0;
    check(cg_zone_write(file.handle, base, b.name.c_str(), size.data(), CGNS_ENUMV(Structured), &z),
          "zone");
    for (const auto& [label, values] :
         {std::pair{"CoordinateX", &b.x}, std::pair{"CoordinateY", &b.y},
          std::pair{"CoordinateZ", &b.z}}) {
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
    if (!b.family.empty()) {
      check(cg_goto(file.handle, base, "Zone_t", z, "end"), "goto zone");
      check(cg_famname_write(b.family.c_str()), "zone family");
    }
    zone_ids.push_back(z);
  }
  for (const CgnsBoundary3D& bc : boundaries) {
    check_name(bc.name, "BC");
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
  const int handle = file.handle;
  file.handle = -1;
  check(cg_close(handle), "close");
}

#endif

void write_cgns_structured_2d(const std::string& path, const std::string& base_name,
                              const std::vector<CgnsZone>& zones,
                              const std::vector<CgnsFamily>& families,
                              const std::vector<CgnsBoundary>& boundaries,
                              const std::vector<CgnsConnection>& connections) {
  write_cgns_2d(path, base_name, zones, {}, families, boundaries, {}, connections, {});
}

}  // namespace og
