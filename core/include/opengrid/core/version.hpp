// OpenGrid core — library version and a trivial smoke-test function.
#pragma once

#include <cstdint>
#include <string_view>

namespace og {

// Semantic version of the core library, e.g. "0.0.1".
std::string_view version();

// Trivial function used to verify the build and the Python bindings (Phase 0).
std::int64_t add(std::int64_t a, std::int64_t b);

}  // namespace og
