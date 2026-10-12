// STEP / IGES import through OpenCASCADE (OCCT, LGPL-2.1 with exception, linked dynamically).
// This optional backend is built only with OG_WITH_OCCT; og_core never depends on it.
//
// Method:
//  1. Read: STEPControl_Reader / IGESControl_Reader (TransferRoots, OneShape). Coordinates are
//     OCCT's session unit (millimetres unless the reader is told otherwise), then scaled by
//     `scale`.
//  2. Heal (optional, each step logged):
//     * ShapeFix_Shape at precision `tolerance` (wires, edges, gaps, orientation);
//     * BRepBuilderAPI_Sewing at `tolerance`, which merges coincident edges of separate faces
//       (IGES files usually arrive as unsewn faces);
//     * ShapeFix_Wireframe: drops edges shorter than `tolerance` and closes wire gaps.
//  3. Split (optional): ShapeUpgrade_ShapeDivideClosed splits periodic faces (a full cylinder
//     becomes two halves, so the seam turns into an ordinary edge), and
//     ShapeUpgrade_ShapeDivideClosedEdges splits closed edges (a full circle becomes two arcs).
//     The OpenGrid meshers need simply parametrised faces and open connectors.
//  4. Topology: unique vertices, edges and faces (TopExp::MapShapes; edges shared by sewn faces
//     are one entry). Each face lists its wires as ordered edges (BRepTools_WireExplorer) with
//     their orientation; degenerate edges (sphere poles, cone apexes) are left out and flagged.
//  5. Geometry: every edge's 3D curve, trimmed to its range, becomes a clamped B-spline
//     (GeomConvert::CurveToBSplineCurve, falling back to GeomConvert_ApproxCurve); its end poles
//     are set to the vertex points so that adjacent edges meet exactly. Every face's surface,
//     trimmed to the face's (u, v) bounds, becomes a B-spline surface
//     (GeomConvert::SurfaceToBSplineSurface, falling back to GeomConvert_ApproxSurface).
// Reference: OCCT 8.0 user guides "Shape Healing" and "Modeling Data" (Open CASCADE SAS);
// L. Piegl & W. Tiller, The NURBS Book, 2nd ed. (1997) for the B-spline forms.
#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace og::cad {

// Clamped B-spline curve: flat knot vector (size poles + degree + 1), weights all 1 when the
// curve is polynomial.
struct BSplineCurveData {
  int degree = 1;
  std::vector<double> knots;
  std::vector<std::array<double, 3>> poles;
  std::vector<double> weights;
};

// Tensor-product B-spline surface: poles[iu * nv + iv].
struct BSplineSurfaceData {
  int degree_u = 1;
  int degree_v = 1;
  std::vector<double> knots_u;
  std::vector<double> knots_v;
  int nu = 0;
  int nv = 0;
  std::vector<std::array<double, 3>> poles;
  std::vector<double> weights;
};

struct CadEdge {
  BSplineCurveData curve;
  std::int64_t v0 = -1;  // vertex at the curve's start
  std::int64_t v1 = -1;  // vertex at the curve's end
  double length = 0.0;
  double tolerance = 0.0;  // OCCT edge tolerance
  int faces = 0;           // faces using the edge (1: free edge of an open shell)
  std::string type;        // OCCT curve type: line, circle, ellipse, bspline, ...
};

struct CadWire {
  std::vector<std::int64_t> edges;  // indices into CadModel::edges, in wire order
  std::vector<bool> reversed;       // edge traversed against its curve
  bool outer = false;
};

struct CadFace {
  BSplineSurfaceData surface;
  bool has_surface = false;  // false when the surface could not be converted
  std::string type;          // OCCT surface type: plane, cylinder, cone, sphere, torus, ...
  std::vector<CadWire> wires;
  bool degenerate = false;  // has a degenerate edge (pole / apex)
};

struct CadOptions {
  double tolerance = 1e-6;   // healing / sewing precision, in the file's (scaled) units
  double scale = 1.0;        // applied to the shape after reading
  bool heal = true;          // ShapeFix_Shape + ShapeFix_Wireframe
  bool sew = true;           // BRepBuilderAPI_Sewing
  bool split_closed = true;  // split periodic faces and closed edges
};

struct CadModel {
  std::vector<std::array<double, 3>> vertices;
  std::vector<CadEdge> edges;
  std::vector<CadFace> faces;
  int solids = 0;
  int shells = 0;
  double max_tolerance = 0.0;    // largest vertex / edge tolerance after healing
  std::vector<std::string> log;  // what was read, healed and split, and conversion warnings
};

// Throws std::runtime_error when the file cannot be read or holds no shape, and
// std::invalid_argument for an unknown extension (.step / .stp / .iges / .igs) or bad options.
CadModel read_cad(const std::string& path, const CadOptions& options = {});

}  // namespace og::cad
