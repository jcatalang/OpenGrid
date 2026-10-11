// CGNS export of 2D grids: one base (cell and physical dimension 2) with structured zones (one
// per block) and unstructured zones (one per unstructured domain), Family_t nodes, boundary
// conditions (BC_t, FamilySpecified + FamilyName) and connectivity.
//  * Structured zones: BCs as PointRange; 1-to-1 connectivity (GridConnectivity1to1_t). Index
//    ranges are 1-based (CGNS convention): {imin, jmin, imax, jmax}.
//  * Unstructured zones: element sections "Triangles" (TRI_3) and "Quadrilaterals" (QUAD_4),
//    numbered 1.. in that order, then one BAR_2 section per boundary condition; the BC_t refers
//    to its bar elements by PointRange with GridLocation EdgeCenter (faces of a 2D zone are
//    edges; SIDS "unstructured BCs as element ranges", CGNS >= 3.1).
//  * Interfaces involving an unstructured zone: general connectivity (GridConnectivity_t,
//    Abutting1to1, GridLocation Vertex) with PointList / PointListDonor of matching vertices
//    (unstructured zones: vertex numbers; structured zones: (i, j) pairs), optionally periodic.
// Written with the CGNS mid-level library. Reference: CFD General Notation System, Standard
// Interface Data Structures (SIDS), and the CGNS Mid-Level Library documentation.
#pragma once

#include <array>
#include <cstdint>
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

struct CgnsUnstructuredZone {
  std::string name;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<std::array<std::int64_t, 3>> triangles;  // 0-based vertex ids, counter-clockwise
  std::vector<std::array<std::int64_t, 4>> quads;
  std::string family;
};

struct CgnsEdgeBoundary {
  int zone = 0;  // index into the unstructured zones
  std::string name;
  std::string family;
  std::vector<std::array<std::int64_t, 2>> edges;  // 0-based vertex ids
};

// Zones are numbered globally: structured zones first, then unstructured ones.
struct CgnsVertexConnection {
  int zone = 0;
  std::string name;
  int donor_zone = 0;
  // 1-based point indices: one per vertex for unstructured zones, (i, j) pairs for structured
  // zones; points[k] matches donor_points[k].
  std::vector<std::int64_t> points;
  std::vector<std::int64_t> donor_points;
  bool periodic = false;
  std::array<double, 2> translation{};
};

// 3D structured blocks (base with cell and physical dimension 3): ranges are
// {imin, jmin, kmin, imax, jmax, kmax}, 1-based; the 1-to-1 transform has three signed entries.
struct CgnsBlock {
  std::string name;
  int ni = 0;
  int nj = 0;
  int nk = 0;
  std::vector<double> x;  // ni * nj * nk, i fastest, then j
  std::vector<double> y;
  std::vector<double> z;
  std::string family;
};

struct CgnsBoundary3D {
  int zone = 0;
  std::string name;
  std::string family;
  std::array<int, 6> range{};
};

struct CgnsConnection3D {
  int zone = 0;
  std::string name;
  int donor_zone = 0;
  std::array<int, 6> range{};
  std::array<int, 6> donor_range{};
  std::array<int, 3> transform{};
};

// 3D unstructured zones (base with cell and physical dimension 3): element sections numbered 1..
// in the given order; a section with a boundary family is a boundary section (TRI_3 / QUAD_4)
// referenced by a BC_t (PointRange of its elements, GridLocation FaceCenter).
struct CgnsSection {
  std::string name;
  std::string type;                 // TETRA_4, PYRA_5, PENTA_6, HEXA_8, TRI_3, QUAD_4
  std::vector<std::int64_t> nodes;  // 0-based vertex ids, element after element
  std::string bc_family;            // non-empty: boundary condition of this section
};

struct CgnsUnstructuredZone3D {
  std::string name;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<CgnsSection> sections;
  std::string family;
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

void write_cgns_2d(const std::string& path, const std::string& base_name,
                   const std::vector<CgnsZone>& zones,
                   const std::vector<CgnsUnstructuredZone>& unstructured_zones,
                   const std::vector<CgnsFamily>& families,
                   const std::vector<CgnsBoundary>& boundaries,
                   const std::vector<CgnsEdgeBoundary>& edge_boundaries,
                   const std::vector<CgnsConnection>& connections,
                   const std::vector<CgnsVertexConnection>& vertex_connections);

void write_cgns_structured_3d(const std::string& path, const std::string& base_name,
                              const std::vector<CgnsBlock>& blocks,
                              const std::vector<CgnsFamily>& families,
                              const std::vector<CgnsBoundary3D>& boundaries,
                              const std::vector<CgnsConnection3D>& connections);

// vertex_connections as in write_cgns_2d (1-based vertex numbers on both sides).
void write_cgns_unstructured_3d(const std::string& path, const std::string& base_name,
                                const std::vector<CgnsUnstructuredZone3D>& zones,
                                const std::vector<CgnsFamily>& families,
                                const std::vector<CgnsVertexConnection>& vertex_connections);

}  // namespace og
