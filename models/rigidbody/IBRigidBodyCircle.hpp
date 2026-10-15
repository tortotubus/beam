#pragma once

#include "elff/models/ibm/IBForceCoupled.hpp"
#include "elff/models/rigidbody/RigidBodyCircle.hpp"

namespace ELFF {
namespace Models {

/**
 * @brief Force-coupled immersed-boundary wrapper for @ref RigidBodyCircle.
 *
 * Basilisk supplies per-marker force densities through @ref IBForceCoupled;
 * this wrapper advances the rigid-body state and mirrors the resulting marker
 * positions, velocities, and quadrature measures into the IB mesh interface.
 */
class IBRigidBodyCircle
  : public IBForceCoupled
  , public RigidBodyCircle
{
public:
  using Vec3 = RigidBodyModel::Vec3;
  using RigidBodyCircle::mesh;

  IBRigidBodyCircle(real_t radius,
                    size_t point_count,
                    real_t density = 1.0,
                    const Vec3& center = { 0., 0., 0. },
                    const Vec3& velocity = { 0., 0., 0. },
                    real_t angle = 0.,
                    real_t angular_velocity_z = 0.);

  void ComputeNextPoints(std::vector<IBMesh::IBVertex> force,
                         real_t dt) override;

  void pack_state(IBModelState& s) const override;
  void unpack_state(const IBModelState& s) override;

private:
  void RigidBodyToIBMesh(IBMesh& ib_mesh) const;
  void RigidBodyToIBMeshes();
};

} // namespace Models
} // namespace ELFF
