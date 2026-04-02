#pragma once

#include "elff/models/beam/IBEulerBeamPenalty.hpp"

namespace ELFF {
namespace Models {

/**
 * @brief Backward-compatible alias wrapper for the penalty IB Euler beam.
 *
 * This preserves the legacy `IBEulerBeam` name while delegating all behavior
 * to @ref IBEulerBeamPenalty.
 */
class IBEulerBeam
  : public IBEulerBeamPenalty
{
public:
  IBEulerBeam(real_t length,
              real_t EI,
              real_t mu,
              size_t nodes,
              EulerBeamBCs bcs,
              real_t r_penalty);
};

} // namespace Models 
} // namespace ELFF 
