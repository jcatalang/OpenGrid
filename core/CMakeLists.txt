add_library(og_core STATIC
  src/version.cpp
  src/util/gauss_kronrod.cpp
  src/util/roots.cpp
  src/util/predicates.cpp
  src/unstructured/size_field.cpp
  src/unstructured/triangulate.cpp
  src/unstructured/advancing_front.cpp
  src/dist/distribution.cpp
  src/dist/geometric.cpp
  src/dist/vinokur.cpp
  src/geom/curve.cpp
  src/geom/line.cpp
  src/geom/arc.cpp
  src/geom/polyline.cpp
  src/geom/cubic_spline.cpp
  src/geom/nurbs.cpp
  src/geom/curve_chain.cpp
  src/topo/connector_nodes.cpp
  src/structured/tfi.cpp
  src/structured/elliptic.cpp
  src/quality/quality2d.cpp
  src/quality/quality_tri.cpp
  src/io/cgns_structured.cpp)
add_library(OpenGrid::core ALIAS og_core)

target_include_directories(og_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(og_core PRIVATE og_options)
target_compile_definitions(og_core PRIVATE
  OG_VERSION_STRING="${PROJECT_VERSION}")
# The core library is linked into the Python extension module.
set_target_properties(og_core PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Optional CGNS support (zlib-like license; see docs/ARCHITECTURE.md).
if(OG_WITH_CGNS)
  # The CGNS package refers to HDF5's imported targets (hdf5-shared): load HDF5's config first.
  find_package(HDF5 CONFIG QUIET)
  find_package(CGNS CONFIG QUIET)
  if(TARGET CGNS::cgns_shared)
    target_link_libraries(og_core PRIVATE CGNS::cgns_shared)
    target_compile_definitions(og_core PRIVATE OG_WITH_CGNS)
    message(STATUS "OpenGrid: CGNS export enabled")
  else()
    message(WARNING "OpenGrid: CGNS library not found; CGNS export disabled")
  endif()
endif()
