// Wall-clock accumulators for the phases of an iteration, enabled by the OG_PROFILE
// environment variable (otherwise every call is a branch on a flag); the report goes with the
// result's warnings so Python logs it.
#pragma once

#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

namespace og {

class PhaseTimer {
 public:
  explicit PhaseTimer(std::vector<std::string> names)
      : names_(std::move(names)),
        seconds_(names_.size(), 0.0),
        enabled_(std::getenv("OG_PROFILE") != nullptr) {}

  [[nodiscard]] bool enabled() const {
    return enabled_;
  }

  void start() {
    if (enabled_) {
      mark_ = std::chrono::steady_clock::now();
    }
  }

  void lap(std::size_t phase) {
    if (!enabled_) {
      return;
    }
    const auto now = std::chrono::steady_clock::now();
    seconds_[phase] += std::chrono::duration<double>(now - mark_).count();
    mark_ = now;
  }

  [[nodiscard]] std::string report() const {
    std::string out;
    double total = 0.0;
    for (const double s : seconds_) {
      total += s;
    }
    for (std::size_t k = 0; k < names_.size(); ++k) {
      if (!out.empty()) {
        out += ", ";
      }
      std::array<char, 64> buffer{};
      std::snprintf(buffer.data(), buffer.size(), "%s %.2f s (%.0f%%)", names_[k].c_str(),
                    seconds_[k], 100.0 * seconds_[k] / (total + 1e-300));
      out += buffer.data();
    }
    return out;
  }

 private:
  std::vector<std::string> names_;
  std::vector<double> seconds_;
  bool enabled_;
  std::chrono::steady_clock::time_point mark_;
};

}  // namespace og
