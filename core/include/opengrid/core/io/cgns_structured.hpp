// CGNS export of structured 2D multi-block grids (one base, one structured zone per block)
// with Family_t nodes, boundary conditions (BC_t with FamilyName) and 1-to-1 connectivity
// (GridConnectivity1to1_t). Index ranges are 1-based (CGNS convention): {imin, jmin, imax, jmax}.
// Written with the CGNS mid-level library. Reference: CFD General Notation System, Standard
// Interface Data Structures (SIDS), and the CGNS Mid-Level Library documentation.
#pragma once

#include <array>
#include <string>
#include <vector>

namespace og {

struct CgnsZone {
  std::string name;
  int ni = 0;
  int nj = 0;
  std::vector<double> x;  // ni * nj, i fastest
  std::vector<double> y;
  std::string family;  // optional zone family (volume tag)
};

struct CgnsFamily {
  std::string name;
  std::string bc_type;  // e.g. "BCWall"; empty for a zone (volume) family
};

struct CgnsBoundary {
  int zone = 0;  // index into the zone list
  std::string name;
  std::string family;  // BCType FamilySpecified + FamilyName
  std::array<int, 4> range{};
};

struct CgnsConnection {
  int zone = 0;
  std::string name;
  int donor_zone = 0;
  std::array<int, 4> range{};
  std::array<int, 4> donor_range{};
  std::array<int, 2> transform{};
  // Translational periodicity (GridConnectivityProperty/Periodic): the donor is this face
  // shifted by `translation`. A 2D base (PhysicalDimension 2) cannot describe a rotation about
  // z, so rotational periodicity is rejected by the caller.
  bool periodic = false;
  std::array<double, 2> translation{};
};

// True when OpenGrid was built with the CGNS library.
bool cgns_available();

// Throws std::runtime_error (with the CGNS library message) on failure, std::invalid_argument
// for unknown BC type names, names longer than 32 characters or inconsistent sizes.
void write_cgns_structured_2d(const std::string& path, const std::string& base_name,
                              const std::vector<CgnsZone>& zones,
                              const std::vector<CgnsFamily>& families,
                              const std::vector<CgnsBoundary>& boundaries,
                              const std::vector<CgnsConnection>& connections);

}  // namespace og
