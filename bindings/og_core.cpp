// nanobind module `opengrid._core`: Python view of the C++ core library.
#include <nanobind/nanobind.h>
#include <nanobind/stl/string_view.h>

#include "opengrid/core/version.hpp"

namespace nb = nanobind;

NB_MODULE(_core, m) {
  m.doc() = "OpenGrid C++ core";
  m.def("version", &og::version, "Version string of the C++ core library.");
  m.def("add", &og::add, nb::arg("a"), nb::arg("b"), "Add two 64-bit integers (build smoke test).");
}
