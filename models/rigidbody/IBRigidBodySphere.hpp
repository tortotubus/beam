#pragma once

#include "elff/models/ibm/IBForceCoupled.hpp"
#include "elff/models/rigidbody/RigidBodySphere.hpp"

namespace ELFF {
namespace Models {

/**
 * @brief Force-coupled immersed-boundary wrapper for @ref RigidBodySphere.
 *
 * Basilisk supplies per-marker force densities through @ref IBForceCoupled;
 * this wrapper advances the rigid-body state and mirrors the resulting marker
 * positions, velocities, and quadrature measures into the IB mesh interface.
 */
class IBRigidBodySphere
  : public IBForceCoupled
  , public RigidBodySphere
{
public:
  using Vec3 = RigidBodyModel::Vec3;
  using Quat = std::array<real_t, 4>;
  using RigidBodySphere::mesh;

  IBRigidBodySphere(real_t radius,
                    size_t point_count,
                    real_t density = 1.0,
                    const Vec3& center = { 0., 0., 0. },
                    const Vec3& velocity = { 0., 0., 0. },
                    const Quat& q_bw = { 1., 0., 0., 0. },
                    const Vec3& angular_velocity_body = { 0., 0., 0. });

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
