#include "opengrid/backends/occt/cad_import.hpp"

#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepLib.hxx>
#include <BRepTools.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Tool.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <GeomConvert.hxx>
#include <GeomConvert_ApproxCurve.hxx>
#include <GeomConvert_ApproxSurface.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <IGESControl_Reader.hxx>
#include <Message.hxx>
#include <Message_Gravity.hxx>
#include <Message_Messenger.hxx>
#include <Message_Printer.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <NCollection_IndexedMap.hxx>
#include <NCollection_List.hxx>
#include <Precision.hxx>
#include <STEPControl_Reader.hxx>
#include <ShapeAnalysis.hxx>
#include <ShapeExtend_Status.hxx>
#include <ShapeFix_Shape.hxx>
#include <ShapeFix_Wireframe.hxx>
#include <ShapeUpgrade_ShapeDivideClosed.hxx>
#include <ShapeUpgrade_ShapeDivideClosedEdges.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace og::cad {

namespace {

using ShapeMap = NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>;
using AncestorMap = NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
                                               TopTools_ShapeMapHasher>;

// Mutes OCCT's console output (reader statistics) while it lives; restores the trace levels.
class QuietMessenger {
 public:
  QuietMessenger() {
    for (const auto& printer : Message::DefaultMessenger()->Printers()) {
      levels_.push_back(printer->GetTraceLevel());
      printer->SetTraceLevel(Message_Fail);
    }
  }
  QuietMessenger(const QuietMessenger&) = delete;
  QuietMessenger& operator=(const QuietMessenger&) = delete;
  ~QuietMessenger() {
    std::size_t k = 0;
    for (const auto& printer : Message::DefaultMessenger()->Printers()) {
      if (k < levels_.size()) {
        printer->SetTraceLevel(levels_[k]);
      }
      ++k;
    }
  }

 private:
  std::vector<Message_Gravity> levels_;
};

std::string lower_extension(const std::string& path) {
  const auto dot = path.find_last_of('.');
  if (dot == std::string::npos) {
    return "";
  }
  std::string ext = path.substr(dot + 1);
  std::transform(ext.begin(), ext.end(), ext.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return ext;
}

TopoDS_Shape read_shape(const std::string& path, std::vector<std::string>& log) {
  const std::string ext = lower_extension(path);
  if (ext == "step" || ext == "stp") {
    STEPControl_Reader reader;
    if (reader.ReadFile(path.c_str()) != IFSelect_RetDone) {
      throw std::runtime_error("read_cad: cannot read STEP file " + path);
    }
    const int roots = reader.TransferRoots();
    log.emplace_back("read STEP: " + std::to_string(roots) + " root(s) transferred");
    return reader.OneShape();
  }
  if (ext == "iges" || ext == "igs") {
    IGESControl_Reader reader;
    if (reader.ReadFile(path.c_str()) != IFSelect_RetDone) {
      throw std::runtime_error("read_cad: cannot read IGES file " + path);
    }
    const int roots = reader.TransferRoots();
    log.emplace_back("read IGES: " + std::to_string(roots) + " root(s) transferred");
    return reader.OneShape();
  }
  throw std::invalid_argument("read_cad: unknown extension of " + path +
                              " (expected .step, .stp, .iges or .igs)");
}

int count(const TopoDS_Shape& shape, TopAbs_ShapeEnum type) {
  ShapeMap map;
  TopExp::MapShapes(shape, type, map);
  return map.Extent();
}

std::string counts(const TopoDS_Shape& shape) {
  std::ostringstream text;
  text << count(shape, TopAbs_SOLID) << " solid(s), " << count(shape, TopAbs_SHELL) << " shell(s), "
       << count(shape, TopAbs_FACE) << " face(s), " << count(shape, TopAbs_EDGE) << " edge(s), "
       << count(shape, TopAbs_VERTEX) << " vertex(es)";
  return text.str();
}

TopoDS_Shape heal(const TopoDS_Shape& input, const CadOptions& options,
                  std::vector<std::string>& log) {
  TopoDS_Shape shape = input;
  if (options.heal) {
    ShapeFix_Shape fix(shape);
    fix.SetPrecision(options.tolerance);
    fix.SetMaxTolerance(std::max(options.tolerance * 1000.0, options.tolerance));
    fix.Perform();
    shape = fix.Shape();
    if (fix.Status(ShapeExtend_DONE)) {
      log.emplace_back("ShapeFix_Shape: fixes applied");
    } else {
      log.emplace_back("ShapeFix_Shape: nothing to fix");
    }
  }
  int free_edges = 0;
  AncestorMap edge_faces;
  TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, edge_faces);
  for (int i = 1; i <= edge_faces.Extent(); ++i) {
    if (edge_faces(i).Extent() == 1 &&
        !BRep_Tool::Degenerated(TopoDS::Edge(edge_faces.FindKey(i)))) {
      ++free_edges;
    }
  }
  if (options.sew && free_edges == 0) {
    log.emplace_back("sewing skipped: no free edges (faces already share their edges)");
  }
  if (options.sew && free_edges > 0) {
    // Sew at the larger of the requested tolerance and the file's own tolerances (IGES
    // files commonly declare 1e-5).
    const double tolerance =
        std::max({options.tolerance, BRep_Tool::MaxTolerance(shape, TopAbs_EDGE),
                  BRep_Tool::MaxTolerance(shape, TopAbs_VERTEX)});
    BRepBuilderAPI_Sewing sewing(tolerance);
    sewing.Add(shape);
    sewing.Perform();
    const TopoDS_Shape& sewed = sewing.SewedShape();
    if (!sewed.IsNull()) {
      shape = sewed;
    }
    std::ostringstream text;
    text << "sewing " << free_edges << " free edge(s) at " << tolerance << ": "
         << sewing.NbFreeEdges() << " free edge(s), " << sewing.NbMultipleEdges()
         << " multiple edge(s), " << sewing.NbDegeneratedShapes() << " degenerated shape(s)";
    log.emplace_back(text.str());
  }
  if (options.heal) {
    ShapeFix_Wireframe wireframe(shape);
    wireframe.SetPrecision(options.tolerance);
    wireframe.ModeDropSmallEdges() = true;
    const bool gaps = wireframe.FixWireGaps();
    const bool small = wireframe.FixSmallEdges();
    shape = wireframe.Shape();
    if (gaps) {
      log.emplace_back("ShapeFix_Wireframe: wire gaps closed");
    }
    if (small) {
      log.emplace_back("ShapeFix_Wireframe: edges shorter than the tolerance removed");
    }
  }
  return shape;
}

TopoDS_Shape split_closed(const TopoDS_Shape& input, double tolerance,
                          std::vector<std::string>& log) {
  TopoDS_Shape shape = input;
  ShapeUpgrade_ShapeDivideClosed faces(shape);
  faces.SetNbSplitPoints(1);
  faces.SetPrecision(tolerance);
  if (faces.Perform() && !faces.Result().IsNull()) {
    shape = faces.Result();
    log.emplace_back("closed (periodic) faces split in two");
  }
  ShapeUpgrade_ShapeDivideClosedEdges edges(shape);
  edges.SetNbSplitPoints(1);
  edges.SetPrecision(tolerance);
  if (edges.Perform() && !edges.Result().IsNull()) {
    shape = edges.Result();
    log.emplace_back("closed edges split in two");
  }
  return shape;
}

std::string_view curve_type(GeomAbs_CurveType type) {
  switch (type) {
    case GeomAbs_Line:
      return "line";
    case GeomAbs_Circle:
      return "circle";
    case GeomAbs_Ellipse:
      return "ellipse";
    case GeomAbs_Hyperbola:
      return "hyperbola";
    case GeomAbs_Parabola:
      return "parabola";
    case GeomAbs_BezierCurve:
      return "bezier";
    case GeomAbs_BSplineCurve:
      return "bspline";
    case GeomAbs_OffsetCurve:
      return "offset";
    default:
      return "other";
  }
}

std::string_view surface_type(GeomAbs_SurfaceType type) {
  switch (type) {
    case GeomAbs_Plane:
      return "plane";
    case GeomAbs_Cylinder:
      return "cylinder";
    case GeomAbs_Cone:
      return "cone";
    case GeomAbs_Sphere:
      return "sphere";
    case GeomAbs_Torus:
      return "torus";
    case GeomAbs_BezierSurface:
      return "bezier";
    case GeomAbs_BSplineSurface:
      return "bspline";
    case GeomAbs_SurfaceOfRevolution:
      return "revolution";
    case GeomAbs_SurfaceOfExtrusion:
      return "extrusion";
    case GeomAbs_OffsetSurface:
      return "offset";
    default:
      return "other";
  }
}

std::array<double, 3> xyz(const gp_Pnt& p) {
  return {p.X(), p.Y(), p.Z()};
}

occ::handle<Geom_BSplineCurve> to_bspline(const occ::handle<Geom_Curve>& curve, double first,
                                          double last, double tolerance) {
  occ::handle<Geom_TrimmedCurve> trimmed = new Geom_TrimmedCurve(curve, first, last);
  occ::handle<Geom_BSplineCurve> result;
  try {
    result = GeomConvert::CurveToBSplineCurve(trimmed);
  } catch (const Standard_Failure&) {
    result.Nullify();
  }
  if (result.IsNull()) {
    GeomConvert_ApproxCurve approx(trimmed, tolerance, GeomAbs_C2, 200, 9);
    if (approx.HasResult()) {
      result = approx.Curve();
    }
  }
  if (!result.IsNull() && result->IsPeriodic()) {
    result->SetNotPeriodic();
  }
  return result;
}

BSplineCurveData curve_data(const occ::handle<Geom_BSplineCurve>& b) {
  BSplineCurveData data;
  data.degree = b->Degree();
  const auto& knots = b->KnotSequence();
  for (int i = knots.Lower(); i <= knots.Upper(); ++i) {
    data.knots.push_back(knots(i));
  }
  for (int i = 1; i <= b->NbPoles(); ++i) {
    data.poles.push_back(xyz(b->Pole(i)));
    data.weights.push_back(b->Weight(i));
  }
  return data;
}

occ::handle<Geom_BSplineSurface> to_bspline(const occ::handle<Geom_Surface>& surface,
                                            const TopoDS_Face& face, double tolerance) {
  double u0 = 0.0;
  double u1 = 0.0;
  double v0 = 0.0;
  double v1 = 0.0;
  BRepTools::UVBounds(face, u0, u1, v0, v1);
  double su0 = 0.0;
  double su1 = 0.0;
  double sv0 = 0.0;
  double sv1 = 0.0;
  surface->Bounds(su0, su1, sv0, sv1);
  // Pcurves may overshoot the surface's own bounds by rounding.
  if (!surface->IsUPeriodic() && !Precision::IsInfinite(su0)) {
    u0 = std::max(u0, su0);
  }
  if (!surface->IsUPeriodic() && !Precision::IsInfinite(su1)) {
    u1 = std::min(u1, su1);
  }
  if (!surface->IsVPeriodic() && !Precision::IsInfinite(sv0)) {
    v0 = std::max(v0, sv0);
  }
  if (!surface->IsVPeriodic() && !Precision::IsInfinite(sv1)) {
    v1 = std::min(v1, sv1);
  }
  occ::handle<Geom_BSplineSurface> result;
  try {
    occ::handle<Geom_RectangularTrimmedSurface> trimmed =
        new Geom_RectangularTrimmedSurface(surface, u0, u1, v0, v1);
    try {
      result = GeomConvert::SurfaceToBSplineSurface(trimmed);
    } catch (const Standard_Failure&) {
      result.Nullify();
    }
    if (result.IsNull()) {
      GeomConvert_ApproxSurface approx(trimmed, tolerance, GeomAbs_C1, GeomAbs_C1, 9, 9, 100, 0);
      if (approx.HasResult()) {
        result = approx.Surface();
      }
    }
  } catch (const Standard_Failure&) {
    result.Nullify();
  }
  if (!result.IsNull()) {
    if (result->IsUPeriodic()) {
      result->SetUNotPeriodic();
    }
    if (result->IsVPeriodic()) {
      result->SetVNotPeriodic();
    }
  }
  return result;
}

BSplineSurfaceData surface_data(const occ::handle<Geom_BSplineSurface>& b) {
  BSplineSurfaceData data;
  data.degree_u = b->UDegree();
  data.degree_v = b->VDegree();
  const auto& ku = b->UKnotSequence();
  for (int i = ku.Lower(); i <= ku.Upper(); ++i) {
    data.knots_u.push_back(ku(i));
  }
  const auto& kv = b->VKnotSequence();
  for (int i = kv.Lower(); i <= kv.Upper(); ++i) {
    data.knots_v.push_back(kv(i));
  }
  data.nu = b->NbUPoles();
  data.nv = b->NbVPoles();
  for (int i = 1; i <= data.nu; ++i) {
    for (int j = 1; j <= data.nv; ++j) {
      data.poles.push_back(xyz(b->Pole(i, j)));
      data.weights.push_back(b->Weight(i, j));
    }
  }
  return data;
}

double distance(const std::array<double, 3>& a, const std::array<double, 3>& b) {
  return std::hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

}  // namespace

CadModel read_cad(const std::string& path, const CadOptions& options) {
  if (!(options.tolerance > 0.0) || !(options.scale > 0.0)) {
    throw std::invalid_argument("read_cad: tolerance and scale must be positive");
  }
  const QuietMessenger quiet;
  CadModel model;
  TopoDS_Shape shape;
  try {
    shape = read_shape(path, model.log);
    if (shape.IsNull()) {
      throw std::runtime_error("read_cad: " + path + " holds no shape");
    }
    model.log.emplace_back("as read: " + counts(shape));
    if (options.scale != 1.0) {
      gp_Trsf scaling;
      scaling.SetScale(gp_Pnt(0.0, 0.0, 0.0), options.scale);
      shape = BRepBuilderAPI_Transform(shape, scaling, true).Shape();
      std::ostringstream text;
      text << "scaled by " << options.scale;
      model.log.emplace_back(text.str());
    }
    shape = heal(shape, options, model.log);
    if (options.split_closed) {
      shape = split_closed(shape, options.tolerance, model.log);
    }
    BRepLib::BuildCurves3d(shape);
    model.log.emplace_back("for meshing: " + counts(shape));
  } catch (const Standard_Failure& failure) {
    throw std::runtime_error(std::string("read_cad: OpenCASCADE failed: ") + failure.what());
  }

  model.solids = count(shape, TopAbs_SOLID);
  model.shells = count(shape, TopAbs_SHELL);

  ShapeMap vertices;
  TopExp::MapShapes(shape, TopAbs_VERTEX, vertices);
  for (int i = 1; i <= vertices.Extent(); ++i) {
    const TopoDS_Vertex& vertex = TopoDS::Vertex(vertices(i));
    model.vertices.push_back(xyz(BRep_Tool::Pnt(vertex)));
    model.max_tolerance = std::max(model.max_tolerance, BRep_Tool::Tolerance(vertex));
  }

  ShapeMap edges;
  TopExp::MapShapes(shape, TopAbs_EDGE, edges);
  AncestorMap edge_faces;
  TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, edge_faces);
  std::vector<std::int64_t> compact(static_cast<std::size_t>(edges.Extent()) + 1, -1);
  for (int i = 1; i <= edges.Extent(); ++i) {
    TopoDS_Edge edge = TopoDS::Edge(edges(i));
    edge.Orientation(TopAbs_FORWARD);
    if (BRep_Tool::Degenerated(edge)) {
      continue;
    }
    double first = 0.0;
    double last = 0.0;
    occ::handle<Geom_Curve> curve = BRep_Tool::Curve(edge, first, last);
    if (curve.IsNull()) {
      model.log.emplace_back("edge " + std::to_string(i) + ": no 3D curve, skipped");
      continue;
    }
    occ::handle<Geom_BSplineCurve> bspline;
    try {
      bspline = to_bspline(curve, first, last, options.tolerance);
    } catch (const Standard_Failure&) {
      bspline.Nullify();
    }
    if (bspline.IsNull()) {
      model.log.emplace_back("edge " + std::to_string(i) + ": B-spline conversion failed, skipped");
      continue;
    }
    CadEdge out;
    TopoDS_Vertex start;
    TopoDS_Vertex end;
    TopExp::Vertices(edge, start, end);
    if (!start.IsNull()) {
      out.v0 = vertices.FindIndex(start) - 1;
    }
    if (!end.IsNull()) {
      out.v1 = vertices.FindIndex(end) - 1;
    }
    out.curve = curve_data(bspline);
    // Clamped B-splines interpolate their end poles: move them onto the vertices so that
    // adjacent edges meet exactly (the gap is within the CAD tolerances).
    if (out.v0 >= 0) {
      const auto& p = model.vertices[static_cast<std::size_t>(out.v0)];
      const double gap = distance(out.curve.poles.front(), p);
      if (gap > BRep_Tool::Tolerance(start) + options.tolerance) {
        std::ostringstream text;
        text << "edge " << i << ": start is " << gap << " from its vertex";
        model.log.emplace_back(text.str());
      }
      out.curve.poles.front() = p;
    }
    if (out.v1 >= 0) {
      const auto& p = model.vertices[static_cast<std::size_t>(out.v1)];
      const double gap = distance(out.curve.poles.back(), p);
      if (gap > BRep_Tool::Tolerance(end) + options.tolerance) {
        std::ostringstream text;
        text << "edge " << i << ": end is " << gap << " from its vertex";
        model.log.emplace_back(text.str());
      }
      out.curve.poles.back() = p;
    }
    BRepAdaptor_Curve adaptor(edge);
    out.type = std::string(curve_type(adaptor.GetType()));
    out.length = GCPnts_AbscissaPoint::Length(adaptor);
    out.tolerance = BRep_Tool::Tolerance(edge);
    model.max_tolerance = std::max(model.max_tolerance, out.tolerance);
    if (edge_faces.Contains(edge)) {
      out.faces = edge_faces.FindFromKey(edge).Extent();
    }
    compact[static_cast<std::size_t>(i)] = static_cast<std::int64_t>(model.edges.size());
    model.edges.push_back(std::move(out));
  }

  ShapeMap faces;
  TopExp::MapShapes(shape, TopAbs_FACE, faces);
  for (int i = 1; i <= faces.Extent(); ++i) {
    TopoDS_Face face = TopoDS::Face(faces(i));
    face.Orientation(TopAbs_FORWARD);
    CadFace out;
    BRepAdaptor_Surface adaptor(face, false);
    out.type = std::string(surface_type(adaptor.GetType()));
    occ::handle<Geom_Surface> surface = BRep_Tool::Surface(face);
    if (!surface.IsNull()) {
      occ::handle<Geom_BSplineSurface> bspline = to_bspline(surface, face, options.tolerance);
      if (!bspline.IsNull()) {
        out.surface = surface_data(bspline);
        out.has_surface = true;
      }
    }
    if (!out.has_surface) {
      model.log.emplace_back("face " + std::to_string(i) + ": B-spline conversion failed");
    }
    const TopoDS_Wire outer = ShapeAnalysis::OuterWire(face);
    for (TopExp_Explorer w(face, TopAbs_WIRE); w.More(); w.Next()) {
      const TopoDS_Wire& wire = TopoDS::Wire(w.Current());
      CadWire cw;
      cw.outer = wire.IsSame(outer);
      for (BRepTools_WireExplorer e(wire, face); e.More(); e.Next()) {
        const TopoDS_Edge& edge = e.Current();
        if (BRep_Tool::Degenerated(edge)) {
          out.degenerate = true;
          continue;
        }
        const int index = edges.FindIndex(edge);
        std::int64_t mapped = -1;
        if (index > 0) {
          mapped = compact[static_cast<std::size_t>(index)];
        }
        if (mapped < 0) {
          model.log.emplace_back("face " + std::to_string(i) + ": a wire edge has no curve");
          continue;
        }
        cw.edges.push_back(mapped);
        cw.reversed.push_back(e.Orientation() == TopAbs_REVERSED);
      }
      if (!cw.edges.empty()) {
        out.wires.push_back(std::move(cw));
      }
    }
    // Outer wire first.
    std::stable_partition(out.wires.begin(), out.wires.end(),
                          [](const CadWire& wire) { return wire.outer; });
    model.faces.push_back(std::move(out));
  }
  return model;
}

}  // namespace og::cad
