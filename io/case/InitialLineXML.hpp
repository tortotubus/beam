#pragma once
#include "elff/io/CXX/xmlCase.hpp"
#include "elff/io/case/ModelPlacement.hpp"

namespace ELFF::IO::Case {
/** @brief Optional child read, preserving malformed/duplicate-node diagnostics. */
std::optional<CXX::XMLNode> optional_xml_child(const CXX::XMLNode&, const char*);
/** @brief Locate the initial line in a v1 or v2 model definition. */
CXX::XMLNode initial_line_xml(const CXX::XMLNode& model);
/** @brief Read finite sampling bounds or a domain with positive active extents. */
PlacementBox parse_placement_box(const CXX::XMLNode&, unsigned dimensions, bool domain);
/** @brief Read owned fixed/random placement, applying beam or rod anchor rules. */
InitialLineTemplate parse_initial_line_template(
    const CXX::XMLNode& model, bool huang, unsigned dimensions, bool allow_random);
/** @brief Read fixed placement, rejecting unresolved random quantities. */
InitialLineDefinition parse_fixed_initial_line(const CXX::XMLNode& model, bool huang);
} // namespace ELFF::IO::Case
