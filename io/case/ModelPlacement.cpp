#include "elff/io/case/ModelPlacement.hpp"
#include "elff/general/error.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace ELFF::IO::Case {
InitialLineDefinition InitialLineTemplate::sample(const CounterRandom& random,
    std::uint64_t attempt, unsigned dimensions) const {
  ELFF_ASSERT(attempt <= (std::numeric_limits<std::uint64_t>::max() - 4) / 5,
              "Placement counter overflow.");
  auto result = fixed;
  const auto base = 5 * attempt;
  if (anchors)
    for (unsigned d = 0; d < dimensions; ++d)
      result.anchor[d] = anchors->min[d] + (anchors->max[d] - anchors->min[d]) * random.uniform01(base + d);
  if (random_direction) {
    constexpr double pi = 3.14159265358979323846;
    const double angle = 2 * pi * random.uniform01(base + 3);
    const double z = dimensions == 3 ? 2 * random.uniform01(base + 4) - 1 : 0;
    const double radius = std::sqrt(std::max(0., 1 - z*z));
    result.direction = {real_t(radius * std::cos(angle)), real_t(radius * std::sin(angle)), real_t(z)};
  }
  return result;
}
} // namespace ELFF::IO::Case
