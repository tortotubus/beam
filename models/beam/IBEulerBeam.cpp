#include "elff/models/beam/IBEulerBeam.hpp"

namespace ELFF {
namespace Models {

IBEulerBeam::IBEulerBeam(real_t length,
                         real_t EI,
                         real_t mu,
                         size_t nodes,
                         EulerBeamBCs bcs,
                         real_t r_penalty)
  : IBEulerBeamPenalty(length, EI, mu, nodes, bcs, r_penalty)
{
}

} // namespace Models
} // namespace ELFF
