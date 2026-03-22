#pragma once

#include "models/ibm/IBForceCoupled.hpp"
#include "models/beam/EulerBeamDynamicInextensibleMoMSparse.hpp"


namespace ELFF {
namespace Models {

/**
 * @brief The immersed boundary force-coupling class for @ref EulerBeam
 */
class IBEulerBeam
  : public IBForceCoupled
  , public EulerBeamDynamicInextensibleMoMSparse
{
private:
  /**
   * @brief Copies the updated beam state into the next immersed-boundary mesh.
   */
  void EBMeshToIBMeshNext();

  /**
   * @brief Copies the updated beam state into the current immersed-boundary
   * mesh.
   */
  void EBMeshToIBMeshCurrent();

public:
  /**
   * @brief Constructs an immersed-boundary-coupled dynamic beam model.
   *
   * @param length Beam length
   * @param EI Flexural rigidity
   * @param mu Mass per unit length
   * @param nodes Number of beam and immersed-boundary nodes
   * @param bcs Boundary conditions applied to the beam
   * @param r_penalty Penalty parameter used in the inextensibility constraint
   */
  IBEulerBeam(real_t length,
              real_t EI,
              real_t mu,
              size_t nodes,
              EulerBeamBCs bcs,
              real_t r_penalty);

  /**
   * @brief Applies initial beam conditions and synchronizes the IB mesh.
   *
   * @param mesh Beam mesh containing the initial geometry
   */
  void apply_initial_condition(EulerBeamMesh& mesh) override;

  /**
   * @brief Advances the coupled beam model using immersed-boundary forces.
   *
   * @param force Nodal forces from the immersed-boundary solver
   * @param dt Time-step size
   */
  void ComputeNextPoints(std::vector<IBMesh::IBVertex> force,
                         real_t dt) override;

  /**
   * @brief Serializes the beam state into an immersed-boundary state buffer.
   *
   * @param s Output state buffer
   */
  void pack_state(IBModelState &s) const override;

  /**
   * @brief Restores the beam state from an immersed-boundary state buffer.
   *
   * @param s Input state buffer
   */
  void unpack_state(const IBModelState& s) override;
};

} // namespace Models 
} // namespace ELFF 
