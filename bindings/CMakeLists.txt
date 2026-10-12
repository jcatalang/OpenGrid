find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development.Module)

# Locate nanobind through the active Python interpreter.
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(_core NB_STATIC og_core.cpp og_3d.cpp)
target_link_libraries(_core PRIVATE og_core og_options)

install(TARGETS _core LIBRARY DESTINATION opengrid)

# Optional CAD import module (backends/occt).
if(TARGET og_occt)
  nanobind_add_module(_occt NB_STATIC og_occt.cpp)
  target_link_libraries(_occt PRIVATE og_occt og_options)
  install(TARGETS _occt LIBRARY DESTINATION opengrid)
endif()
