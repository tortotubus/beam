#pragma once

#include "elff/models/beam/EulerBeam.hpp"
#include "elff/io/C/xmlCase.h"
#include <array>

namespace ELFF::Models {

/** Owned model inputs, independent of XML document and runtime lifetime. */
struct EulerBeamHuangDefinition
{
  long source_line;
  real_t length;
  real_t bending_rigidity;
  real_t linear_density;
  size_t nodes;
  std::array<real_t, 3> anchor;
  std::array<real_t, 3> direction;
  bool anchor_at_right;
  EulerBeam::EulerBeamBCs boundary_conditions;
};

/** Read one already-rendered EulerBeamHuang node and validate model semantics.
 * Throws IO::CXX::XMLError for XML/semantic failures. Does not construct a model,
 * assign ranks, validate a DTD, or retain XML handles. */
EulerBeamHuangDefinition parse_euler_beam_huang_xml(elff_xml_node_t model);

} // namespace ELFF::Models
