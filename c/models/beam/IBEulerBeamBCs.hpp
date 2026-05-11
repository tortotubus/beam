#pragma once

#include "elff/c/models/beam/IBEulerBeamBCs.h"
#include "elff/models/beam/EulerBeam.hpp"

#include <cstddef>

namespace ELFF {
namespace C {

inline Models::EulerBeam::EulerBeamBCVals
to_cpp_beam_bc_vals(const ib_euler_beam_bc_vals_t& vals)
{
  Models::EulerBeam::EulerBeamBCVals cpp_vals = {};
  cpp_vals.position[0] = vals.position.x;
  cpp_vals.position[1] = vals.position.y;
  cpp_vals.position[2] = vals.position.z;
  cpp_vals.slope[0] = vals.slope.x;
  cpp_vals.slope[1] = vals.slope.y;
  cpp_vals.slope[2] = vals.slope.z;
  cpp_vals.force[0] = vals.force.x;
  cpp_vals.force[1] = vals.force.y;
  cpp_vals.force[2] = vals.force.z;
  cpp_vals.torque[0] = vals.torque.x;
  cpp_vals.torque[1] = vals.torque.y;
  cpp_vals.torque[2] = vals.torque.z;
  return cpp_vals;
}

inline Models::EulerBeam::EulerBeamBCs
to_cpp_beam_bcs(const ib_euler_beam_bcs_t& bcs)
{
  Models::EulerBeam::EulerBeamBCs cpp_bcs = {};
  for (std::size_t bi = 0; bi < 2; ++bi) {
    cpp_bcs.end[bi] =
      static_cast<Models::EulerBeam::EulerBeamBCEnd>(bcs.end[bi]);
    cpp_bcs.type[bi] =
      static_cast<Models::EulerBeam::EulerBeamBCType>(bcs.type[bi]);
    cpp_bcs.vals[bi] = to_cpp_beam_bc_vals(bcs.vals[bi]);
  }
  return cpp_bcs;
}

} // namespace C
} // namespace ELFF
