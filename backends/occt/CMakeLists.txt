# Optional OpenCASCADE backend (STEP / IGES import). Built only with OG_WITH_OCCT and when
# OpenCASCADE is found; og_core never depends on it.
set(OG_OCCT_LIBRARIES
  TKernel TKMath TKG2d TKG3d TKGeomBase TKBRep TKGeomAlgo TKTopAlgo TKShHealing
  TKXSBase TKDE TKDESTEP TKDEIGES)

add_library(og_occt STATIC src/cad_import.cpp)
target_include_directories(og_occt PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
# OCCT headers are treated as system headers so our warning flags do not apply to them.
target_include_directories(og_occt SYSTEM PRIVATE ${OpenCASCADE_INCLUDE_DIR})
target_link_libraries(og_occt PRIVATE og_options ${OG_OCCT_LIBRARIES})
set_target_properties(og_occt PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(OG_BUILD_TESTS)
  find_package(GTest CONFIG REQUIRED)
  include(GoogleTest)
  add_executable(og_occt_tests tests/test_cad_import.cpp)
  target_include_directories(og_occt_tests SYSTEM PRIVATE ${OpenCASCADE_INCLUDE_DIR})
  target_link_libraries(og_occt_tests PRIVATE og_occt og_core og_options GTest::gtest_main
    ${OG_OCCT_LIBRARIES} TKPrim TKBO TKBool)
  gtest_discover_tests(og_occt_tests)
endif()
