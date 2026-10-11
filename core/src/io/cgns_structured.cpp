#include "opengrid/core/io/cgns_structured.hpp"

#include <array>
#include <stdexcept>
#include <utility>

#ifdef OG_WITH_CGNS
#include <cgnslib.h>
#endif

namespace og {

#ifndef OG_WITH_CGNS

bool cgns_available() {
  return false;
}

void write_cgns_structured_2d(const std::string&, const std::string&, const std::vector<CgnsZone>&,
                              const std::vector<CgnsFamily>&, const std::vector<CgnsBoundary>&,
                              const std::vector<CgnsConnection>&) {
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

void write_cgns_structured_2d(const std::string& path, const std::string& base_name,
                              const std::vector<CgnsZone>& zones,
                              const std::vector<CgnsFamily>& families,
                              const std::vector<CgnsBoundary>& boundaries,
                              const std::vector<CgnsConnection>& connections) {
  check_name(base_name, "base");
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
  for (const CgnsZone& zone : zones) {
    std::array<cgsize_t, 6> size = {zone.ni, zone.nj, zone.ni - 1, zone.nj - 1, 0, 0};
    int z = 0;
    check(cg_zone_write(file.handle, base, zone.name.c_str(), size.data(), CGNS_ENUMV(Structured),
                        &z),
          "zone");
    for (const auto& [label, values] :
         {std::pair{"CoordinateX", &zone.x}, std::pair{"CoordinateY", &zone.y}}) {
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
    if (!zone.family.empty()) {
      check(cg_goto(file.handle, base, "Zone_t", z, "end"), "goto zone");
      check(cg_famname_write(zone.family.c_str()), "zone family");
    }
    zone_ids.push_back(z);
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
    const std::string& donor = zones.at(static_cast<std::size_t>(c.donor_zone)).name;
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
  const int handle = file.handle;
  file.handle = -1;
  check(cg_close(handle), "close");
}

#endif

}  // namespace og
