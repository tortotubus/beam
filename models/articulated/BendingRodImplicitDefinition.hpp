#pragma once
#include "elff/models/articulated/BendingRodImplicit.hpp"

namespace ELFF::Models {
/** @brief Owned implicit rod inputs, independent of XML and runtime ownership. */
struct BendingRodImplicitDefinition {
  long source_line = 0;
  BendingRodImplicitConfig config;
  /** @brief config.pin stores the initial-line anchor; the builder resolves the physical pin. */
  bool anchor_at_right = false;
  Vec3 direction = Vec3::UnitX();
};
} // namespace ELFF::Models
