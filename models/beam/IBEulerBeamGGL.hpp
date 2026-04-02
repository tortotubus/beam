#pragma once

#include "elff/models/beam/EulerBeamInextensibleGGL.hpp"
#include "elff/models/ibm/IBForceCoupled.hpp"

namespace ELFF {
namespace Models {

/**
 * @brief Immersed-boundary wrapper for the full GGL inextensible Euler beam
 * solver.
 *
 * This class couples @ref EulerBeamInextensibleGGL to the immersed-boundary
 * force interface by translating between IB mesh data and the underlying beam
 * state while preserving restart support through @ref IBModel.
 */
class IBEulerBeamGGL
  : public IBForceCoupled
  , public EulerBeamInextensibleGGL
{
private:
  void EBMeshToIBMeshNext();
  void EBMeshToIBMeshCurrent();

public:
  IBEulerBeamGGL(real_t length,
                 real_t EI,
                 real_t mu,
                 size_t nodes,
                 EulerBeamBCs bcs,
                 real_t r_penalty = 0.0);

  void apply_initial_condition(EulerBeamMesh& mesh) override;

  void ComputeNextPoints(std::vector<IBMesh::IBVertex> force,
                         real_t dt) override;

  void pack_state(IBModelState& s) const override;
  void unpack_state(const IBModelState& s) override;
};

} // namespace Models
} // namespace ELFF
