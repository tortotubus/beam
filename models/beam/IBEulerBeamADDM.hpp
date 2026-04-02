#pragma once

#include "elff/models/beam/EulerBeamInextensibleADDM.hpp"
#include "elff/models/ibm/IBForceCoupled.hpp"

namespace ELFF {
namespace Models {

/**
 * @brief Immersed-boundary wrapper for the ADDM inextensible Euler beam
 * solver.
 *
 * This class couples @ref EulerBeamInextensibleADDM with the force-coupled IB
 * interface by mirroring the Euler-beam centerline state onto IB meshes and by
 * exposing restart support through @ref IBModel.
 */
class IBEulerBeamADDM
  : public IBForceCoupled
  , public EulerBeamInextensibleADDM
{
private:
  void EBMeshToIBMeshNext();
  void EBMeshToIBMeshCurrent();

public:
  IBEulerBeamADDM(real_t length,
                  real_t EI,
                  real_t mu,
                  size_t nodes,
                  EulerBeamBCs bcs,
                  real_t r_penalty);

  void apply_initial_condition(EulerBeamMesh& mesh) override;

  void ComputeNextPoints(std::vector<IBMesh::IBVertex> force,
                         real_t dt) override;

  void pack_state(IBModelState& s) const override;
  void unpack_state(const IBModelState& s) override;
};

} // namespace Models
} // namespace ELFF
