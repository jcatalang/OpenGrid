cmake_minimum_required(VERSION 3.28)

project(OpenGrid VERSION 0.0.1 LANGUAGES CXX)

option(OG_BUILD_TESTS "Build the C++ GoogleTest suite" ON)
option(OG_BUILD_PYTHON "Build the nanobind Python extension" ON)
option(OG_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Apple Clang < 15 ships an incomplete C++20 standard library (see docs/ARCHITECTURE.md).
if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15)
  message(FATAL_ERROR
    "Apple Clang ${CMAKE_CXX_COMPILER_VERSION} is too old (need >= 15). "
    "Select the Command Line Tools, e.g. export DEVELOPER_DIR=/Library/Developer/CommandLineTools "
    "(see README.md).")
endif()
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11)
  message(FATAL_ERROR "GCC ${CMAKE_CXX_COMPILER_VERSION} is too old (need >= 11).")
endif()

# Common compile options for all OpenGrid targets.
add_library(og_options INTERFACE)
target_compile_options(og_options INTERFACE
  -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion
  # Determinism: never fuse multiply-add implicitly (see ARCHITECTURE.md §10).
  -ffp-contract=off)
if(OG_WARNINGS_AS_ERRORS)
  target_compile_options(og_options INTERFACE -Werror)
endif()

add_subdirectory(core)

if(OG_BUILD_PYTHON)
  add_subdirectory(bindings)
endif()

if(OG_BUILD_TESTS)
  enable_testing()
  add_subdirectory(core/tests)
endif()
