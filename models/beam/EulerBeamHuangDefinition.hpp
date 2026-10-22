#pragma once
#include "elff/models/beam/EulerBeam.hpp"

namespace ELFF::Models {
/** @brief Owned Huang material, discretization and boundary settings, independent of XML. */
struct EulerBeamHuangDefinition {
  long source_line = 0;
  real_t length = 0;
  real_t bending_rigidity = 0;
  real_t linear_density = 0;
  size_t nodes = 0;
  /** @brief Boundary types; the builder derives endpoint positions and slopes from placement. */
  EulerBeam::EulerBeamBCs boundary_conditions{};
};
} // namespace ELFF::Models
