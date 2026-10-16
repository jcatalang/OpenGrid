add_library(og_core STATIC
  src/version.cpp
  src/util/gauss_kronrod.cpp
  src/util/roots.cpp
  src/util/predicates.cpp
  src/util/parallel.cpp
  src/util/anderson.cpp
  src/unstructured/size_field.cpp
  src/unstructured/triangulate.cpp
  src/unstructured/advancing_front.cpp
  src/unstructured/recombine.cpp
  src/unstructured/smoothing.cpp
  src/unstructured/domain_mesher.cpp
  src/unstructured/inflation.cpp
  src/unstructured/inflation3d.cpp
  src/unstructured/optimize3d.cpp
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
  src/structured/multiblock_elliptic.cpp
  src/structured/multiblock_elliptic3d.cpp
  src/quality/quality2d.cpp
  src/quality/quality_tri.cpp
  src/quality/quality_cells.cpp
  src/quality/quality_cells3d.cpp
  src/quality/quality3d.cpp
  src/quality/volume_quality.cpp
  src/geom/surface.cpp
  src/structured/tfi3d.cpp
  src/structured/extrude.cpp
  src/structured/surface_tfi.cpp
  src/structured/elliptic3d.cpp
  src/structured/surface_elliptic.cpp
  src/structured/hyperbolic.cpp
  src/structured/wall_layers.cpp
  src/gpu/gpu.cpp
  src/gpu/multiblock3d_sweeps.cpp
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

# Optional OpenCL for the GPU sweeps (core/include/opengrid/core/gpu/gpu.hpp): the macOS
# framework, or an ICD loader elsewhere (conda-forge ocl-icd / khronos-opencl-icd-loader).
if(OG_WITH_OPENCL)
  find_package(OpenCL)
  if(OpenCL_FOUND)
    target_link_libraries(og_core PRIVATE OpenCL::OpenCL)
    target_compile_definitions(og_core PRIVATE OG_WITH_OPENCL CL_TARGET_OPENCL_VERSION=120
                                               CL_SILENCE_DEPRECATION)
    message(STATUS "OpenGrid: OpenCL ${OpenCL_VERSION_STRING} enabled (GPU sweeps)")
  else()
    message(STATUS "OpenGrid: OpenCL not found; GPU sweeps disabled")
  endif()
endif()

# Optional OpenMP for the hot loops (core/include/opengrid/core/util/parallel.hpp). Apple Clang
# has no OpenMP runtime of its own: the conda environment's libomp (llvm-openmp) is used.
if(OG_WITH_OPENMP)
  if(NOT DEFINED OpenMP_ROOT AND DEFINED ENV{CONDA_PREFIX})
    set(OpenMP_ROOT "$ENV{CONDA_PREFIX}")
  endif()
  find_package(OpenMP COMPONENTS CXX)
  if(OpenMP_CXX_FOUND)
    # Public: parallel_for is a header template, so its users compile the pragmas too.
    target_link_libraries(og_core PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(og_core PUBLIC OG_WITH_OPENMP)
    message(STATUS "OpenGrid: OpenMP ${OpenMP_CXX_VERSION} enabled")
  else()
    message(WARNING "OpenGrid: OpenMP not found; parallel loops run serially")
  endif()
endif()
