#pragma once

#include "elff/models/articulated/BendingRodImplicit.hpp"
#include "elff/io/C/xmlCase.h"

namespace ELFF::Models {
/** Owned inputs: no XML handles or runtime ownership. */
struct BendingRodImplicitDefinition {
  long source_line = 0;
  BendingRodImplicitConfig config;
  Vec3 direction = Vec3::UnitX();
};
BendingRodImplicitDefinition parse_bending_rod_implicit_xml(elff_xml_node_t node);
} // namespace ELFF::Models
