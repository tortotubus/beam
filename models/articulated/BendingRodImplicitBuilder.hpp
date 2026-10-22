#pragma once
#include "elff/models/articulated/BendingRodImplicitDefinition.hpp"
#include "elff/models/articulated/IBBendingRodImplicit.hpp"
#include <memory>

namespace ELFF::Models {
inline std::unique_ptr<IBBendingRodImplicit> build_bending_rod_implicit(
    const BendingRodImplicitDefinition &definition)
{
  auto rod = std::make_unique<IBBendingRodImplicit>(definition.config);
  BendingRodMesh initial(definition.config.links + 1, definition.config.length);
  for (size_t i = 0; i < initial.get_nodes(); ++i)
    for (size_t d = 0; d < 3; ++d)
      initial.get_centerline(i)[d] = definition.config.pin[d] +
        initial.get_curvilinear_axis(i) * definition.direction[d];
  rod->apply_initial_condition(initial);
  return rod;
}
} // namespace ELFF::Models
