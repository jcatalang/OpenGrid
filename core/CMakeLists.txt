add_library(og_core STATIC
  src/version.cpp)
add_library(OpenGrid::core ALIAS og_core)

target_include_directories(og_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(og_core PRIVATE og_options)
target_compile_definitions(og_core PRIVATE
  OG_VERSION_STRING="${PROJECT_VERSION}")
# The core library is linked into the Python extension module.
set_target_properties(og_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
