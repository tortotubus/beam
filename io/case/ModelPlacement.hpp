#pragma once
#include "elff/config/config.hpp"
#include "elff/general/random.hpp"
#include <array>
#include <optional>

namespace ELFF::IO::Case {
/** @brief Owned fixed line placement in global coordinates. */
struct InitialLineDefinition {
  std::array<real_t, 3> anchor{}, direction{1, 0, 0};
  bool anchor_at_right = false;
};
/** @brief Owned axis-aligned placement bounds and their source location. */
struct PlacementBox {
  std::array<real_t, 3> min{}, max{};
  long source_line = 0;
};
/** @brief Line template; random quantities are resolved once per attempt. */
struct InitialLineTemplate {
  InitialLineDefinition fixed;
  std::optional<PlacementBox> anchors;
  bool random_direction = false;
  long source_line = 0;
  InitialLineDefinition sample(const CounterRandom&, std::uint64_t attempt,
                               unsigned dimensions) const;
};
} // namespace ELFF::IO::Case
