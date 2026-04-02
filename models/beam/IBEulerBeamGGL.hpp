#pragma once

#include "elff/models/beam/EulerBeamInextensibleGGL.hpp"
#include "elff/models/ibm/IBForceCoupled.hpp"

namespace ELFF {
namespace Models {

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
