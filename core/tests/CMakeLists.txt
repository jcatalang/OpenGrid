find_package(GTest CONFIG REQUIRED)
include(GoogleTest)

add_executable(og_core_tests
  test_version.cpp
  test_quadrature.cpp
  test_curves.cpp
  test_arclength.cpp
  test_projection.cpp
  test_distribution.cpp
  test_connector.cpp)
target_link_libraries(og_core_tests PRIVATE og_core og_options GTest::gtest_main)
gtest_discover_tests(og_core_tests)
