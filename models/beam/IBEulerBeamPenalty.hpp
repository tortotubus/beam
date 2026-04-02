#pragma once

#include "elff/models/beam/EulerBeamInextensiblePenalty.hpp"
#include "elff/models/ibm/IBForceCoupled.hpp"

namespace ELFF {
namespace Models {

/**
 * @brief Immersed-boundary wrapper for the penalty-based inextensible Euler
 * beam solver.
 *
 * This class combines @ref EulerBeamInextensiblePenalty with the
 * force-coupled IB interface so the beam can advance from immersed-boundary
 * loads while keeping its solver state available for restart operations.
 */
class IBEulerBeamPenalty
  : public IBForceCoupled
  , public EulerBeamInextensiblePenalty
{
private:
  void EBMeshToIBMeshNext();
  void EBMeshToIBMeshCurrent();

public:
  IBEulerBeamPenalty(real_t length,
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
