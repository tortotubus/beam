#pragma once
#include "elff/models/articulated/BendingRodImplicitDefinition.hpp"
#include "elff/io/case/ModelPlacement.hpp"
#include "elff/io/C/xmlCase.h"

namespace ELFF::IO::Case {
/** @brief Parse structural settings with an already resolved line placement. */
Models::BendingRodImplicitDefinition parse_bending_rod_implicit_xml(
    elff_xml_node_t model, const InitialLineDefinition& placement);
/** @brief Parse a fixed rod; unresolved random placement is a semantic error. */
Models::BendingRodImplicitDefinition parse_bending_rod_implicit_xml(elff_xml_node_t node);
} // namespace ELFF::IO::Case
