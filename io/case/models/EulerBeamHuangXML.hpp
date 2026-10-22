#pragma once
#include "elff/io/CXX/xmlCase.hpp"
#include "elff/models/beam/EulerBeamHuangDefinition.hpp"
#include "elff/io/case/ModelPlacement.hpp"

namespace ELFF::IO::Case {
/** @brief Read Huang XML into owned settings without constructing or adopting models. */
class EulerBeamHuangXML {
public:
  static constexpr const char* tag = "EulerBeamHuang";
  /** @brief Read and validate material, discretization and boundary settings. */
  static Models::EulerBeamHuangDefinition parse(const CXX::XMLNode&);
  /** @brief Read placement; standalone reads reject unresolved random quantities.
   * Huang requires an explicit left or right anchorEnd. Source diagnostics are
   * preserved, and returned values do not borrow the XML document. */
  static InitialLineTemplate parse_placement(
      const CXX::XMLNode&, const Models::EulerBeamHuangDefinition&,
      unsigned dimensions = 3, bool allow_random = false);
};
} // namespace ELFF::IO::Case
