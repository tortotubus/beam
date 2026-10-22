#pragma once
#include "elff/models/articulated/BendingRodImplicit.hpp"

namespace ELFF::Models {
/** @brief Owned implicit rod inputs, independent of XML and runtime ownership. */
struct BendingRodImplicitDefinition {
  long source_line = 0;
  BendingRodImplicitConfig config;
  Vec3 direction = Vec3::UnitX();
};
} // namespace ELFF::Models
