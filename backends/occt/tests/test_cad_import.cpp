// STEP / IGES import: a box with a cylindrical through hole is built with OCCT, written to
// STEP and IGES, and read back through og::cad::read_cad.
#include <gtest/gtest.h>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <IGESControl_Writer.hxx>
#include <STEPControl_Writer.hxx>
#include <gp_Ax2.hxx>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

#include "opengrid/backends/occt/cad_import.hpp"
#include "opengrid/core/geom/nurbs.hpp"
#include "opengrid/core/geom/surface.hpp"

namespace {

constexpr double kRadius = 0.4;

// Box [0, 2] x [0, 2] x [0, 1] minus the cylinder of radius 0.4 about (1, 1, z).
TopoDS_Shape part() {
  const TopoDS_Shape box = BRepPrimAPI_MakeBox(2.0, 2.0, 1.0).Shape();
  const gp_Ax2 axis(gp_Pnt(1.0, 1.0, -0.5), gp_Dir(0.0, 0.0, 1.0));
  const TopoDS_Shape cylinder = BRepPrimAPI_MakeCylinder(axis, kRadius, 2.0).Shape();
  return BRepAlgoAPI_Cut(box, cylinder).Shape();
}

// A file of the running test (ctest runs the tests in parallel: one name per test, or one test
// reads a file another is writing).
std::string temp_file(const std::string& name) {
  const auto dir = std::filesystem::temp_directory_path() / "og_occt_tests";
  std::filesystem::create_directories(dir);
  const auto* test = ::testing::UnitTest::GetInstance()->current_test_info();
  std::string prefix;
  if (test != nullptr) {
    prefix = std::string(test->test_suite_name()) + "." + test->name() + ".";
  }
  return (dir / (prefix + name)).string();
}

std::string write_step() {
  const std::string path = temp_file("part.step");
  STEPControl_Writer writer;
  writer.Transfer(part(), STEPControl_AsIs);
  writer.Write(path.c_str());
  return path;
}

std::string write_iges() {
  const std::string path = temp_file("part.igs");
  IGESControl_Writer writer("MM", 0);  // faces mode: trimmed surfaces, not sewn
  writer.AddShape(part());
  writer.ComputeModel();
  writer.Write(path.c_str());
  return path;
}

og::Vec3 vec(const std::array<double, 3>& p) {
  return {p[0], p[1], p[2]};
}

// Distance from a point to the part's boundary surfaces (box planes and the hole's cylinder).
double off_part(const og::Vec3& p) {
  double best = std::abs(std::hypot(p.x - 1.0, p.y - 1.0) - kRadius);
  for (const double v : {p.x, 2.0 - p.x, p.y, 2.0 - p.y, p.z, 1.0 - p.z}) {
    best = std::min(best, std::abs(v));
  }
  return best;
}

og::Nurbs curve_of(const og::cad::CadEdge& edge) {
  std::vector<og::Vec3> poles;
  poles.reserve(edge.curve.poles.size());
  for (const auto& p : edge.curve.poles) {
    poles.push_back(vec(p));
  }
  return {edge.curve.degree, edge.curve.knots, poles, edge.curve.weights};
}

og::NurbsSurface surface_of(const og::cad::CadFace& face) {
  std::vector<og::Vec3> poles;
  poles.reserve(face.surface.poles.size());
  for (const auto& p : face.surface.poles) {
    poles.push_back(vec(p));
  }
  return {face.surface.degree_u,
          face.surface.degree_v,
          face.surface.knots_u,
          face.surface.knots_v,
          face.surface.nu,
          face.surface.nv,
          poles,
          face.surface.weights};
}

void expect_closed_part(const og::cad::CadModel& model) {
  ASSERT_FALSE(model.faces.empty());
  // 4 sides + top and bottom (each with a hole) + the hole's cylinder split in two halves.
  EXPECT_EQ(model.faces.size(), 8U);
  int cylinders = 0;
  int holed = 0;
  for (const auto& face : model.faces) {
    EXPECT_TRUE(face.has_surface);
    EXPECT_FALSE(face.degenerate);
    ASSERT_FALSE(face.wires.empty());
    EXPECT_TRUE(face.wires.front().outer);
    if (face.type != "plane") {  // a cylinder (STEP) or surface of revolution (IGES)
      ++cylinders;
    }
    if (face.wires.size() == 2) {
      ++holed;
    }
  }
  EXPECT_EQ(cylinders, 2);
  EXPECT_EQ(holed, 2);
  for (const auto& edge : model.edges) {
    // Closed solid: every edge bounds two faces; no closed edges remain after splitting.
    EXPECT_EQ(edge.faces, 2);
    EXPECT_NE(edge.v0, edge.v1);
    EXPECT_GT(edge.length, 0.0);
  }
}

TEST(CadImport, StepPartTopologyAndGeometry) {
  const og::cad::CadModel model = og::cad::read_cad(write_step());
  expect_closed_part(model);
  EXPECT_EQ(model.solids, 1);
  for (const auto& edge : model.edges) {
    const og::Nurbs curve = curve_of(edge);
    // End poles are the vertices, so the curve ends are exactly the vertex points.
    const og::Vec3 start = curve.eval(0.0);
    const og::Vec3 end = curve.eval(1.0);
    const og::Vec3 v0 = vec(model.vertices[static_cast<std::size_t>(edge.v0)]);
    const og::Vec3 v1 = vec(model.vertices[static_cast<std::size_t>(edge.v1)]);
    EXPECT_LT(og::norm(start - v0), 1e-15);
    EXPECT_LT(og::norm(end - v1), 1e-15);
    for (int k = 0; k <= 10; ++k) {
      EXPECT_LT(off_part(curve.eval(k / 10.0)), 1e-9) << edge.type;
    }
    // Arc length of the converted curve equals OCCT's length of the edge.
    EXPECT_NEAR(curve.length(), edge.length, 1e-9 * edge.length);
  }
  for (const auto& face : model.faces) {
    const og::NurbsSurface surface = surface_of(face);
    for (int i = 0; i <= 6; ++i) {
      for (int j = 0; j <= 6; ++j) {
        const og::Vec3 p = surface.eval(i / 6.0, j / 6.0);
        if (face.type != "plane") {  // a cylinder (STEP) or surface of revolution (IGES)
          EXPECT_NEAR(std::hypot(p.x - 1.0, p.y - 1.0), kRadius, 1e-9);
        } else {
          EXPECT_LT(off_part(p), 1e-9);
        }
      }
    }
  }
}

TEST(CadImport, IgesFacesAreSewn) {
  const std::string path = write_iges();
  const og::cad::CadModel sewn = og::cad::read_cad(path);
  expect_closed_part(sewn);
  og::cad::CadOptions loose;
  loose.sew = false;
  const og::cad::CadModel unsewn = og::cad::read_cad(path, loose);
  const auto free_edges = std::count_if(unsewn.edges.begin(), unsewn.edges.end(),
                                        [](const og::cad::CadEdge& e) { return e.faces == 1; });
  EXPECT_GT(free_edges, 0);
  EXPECT_GT(unsewn.edges.size(), sewn.edges.size());
}

TEST(CadImport, ClosedShapesStayWholeWithoutSplitting) {
  og::cad::CadOptions options;
  options.split_closed = false;
  const og::cad::CadModel model = og::cad::read_cad(write_step(), options);
  // The hole is one periodic face with a seam edge used twice.
  EXPECT_EQ(model.faces.size(), 7U);
  const bool closed_edge = std::any_of(model.edges.begin(), model.edges.end(),
                                       [](const og::cad::CadEdge& e) { return e.v0 == e.v1; });
  EXPECT_TRUE(closed_edge);
}

TEST(CadImport, ScaleAndErrors) {
  const std::string path = write_step();
  og::cad::CadOptions options;
  options.scale = 0.5;
  const og::cad::CadModel model = og::cad::read_cad(path, options);
  double top = 0.0;
  for (const auto& v : model.vertices) {
    top = std::max({top, v[0], v[1], v[2]});
  }
  EXPECT_NEAR(top, 1.0, 1e-12);
  EXPECT_THROW((void)og::cad::read_cad(temp_file("part.txt")), std::invalid_argument);
  EXPECT_THROW((void)og::cad::read_cad(temp_file("missing.step")), std::runtime_error);
  options.tolerance = 0.0;
  EXPECT_THROW((void)og::cad::read_cad(path, options), std::invalid_argument);
}

}  // namespace
