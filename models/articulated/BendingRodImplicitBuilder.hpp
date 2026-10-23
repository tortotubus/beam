#pragma once
#include "elff/models/articulated/BendingRodImplicitDefinition.hpp"
#include "elff/models/articulated/IBBendingRodImplicit.hpp"
#include <memory>

namespace ELFF::Models {
inline std::unique_ptr<IBBendingRodImplicit> build_bending_rod_implicit(
    const BendingRodImplicitDefinition &definition)
{
  auto config = definition.config;
  const bool pin_at_right = config.boundary == BendingRod::Boundary::PinnedRight;
  if (config.boundary != BendingRod::Boundary::Free &&
      pin_at_right != definition.anchor_at_right)
    config.pin += config.length * definition.direction;
  auto rod = std::make_unique<IBBendingRodImplicit>(config);
  BendingRodMesh initial(definition.config.links + 1, definition.config.length);
  for (size_t i = 0; i < initial.get_nodes(); ++i)
    for (size_t d = 0; d < 3; ++d)
      initial.get_centerline(i)[d] = definition.config.pin[d] +
        (definition.anchor_at_right ? definition.config.length - initial.get_curvilinear_axis(i)
                                    : initial.get_curvilinear_axis(i)) * definition.direction[d];
  rod->apply_initial_condition(initial);
  return rod;
}
} // namespace ELFF::Models
