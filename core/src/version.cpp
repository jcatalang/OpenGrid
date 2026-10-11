#include "opengrid/core/version.hpp"

namespace og {

std::string_view version() {
  return OG_VERSION_STRING;
}

std::int64_t add(std::int64_t a, std::int64_t b) {
  return a + b;
}

}  // namespace og
