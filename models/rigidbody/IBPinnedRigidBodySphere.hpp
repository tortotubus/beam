#pragma once

#include "elff/models/rigidbody/IBRigidBodySphere.hpp"

namespace ELFF {
namespace Models {

/**
 * @brief Immersed-boundary spherical rigid body with fixed marker positions.
 *
 * This is useful for fixed immersed obstacles that should expose the same
 * marker geometry as @ref IBRigidBodySphere while ignoring fluid forces.
 */
class IBPinnedRigidBodySphere : public IBRigidBodySphere
{
public:
  using Vec3 = RigidBodyModel::Vec3;
  using Quat = IBRigidBodySphere::Quat;

  IBPinnedRigidBodySphere(real_t radius,
                          size_t point_count,
                          real_t density = 1.0,
                          const Vec3& center = { 0., 0., 0. },
                          const Quat& q_bw = { 1., 0., 0., 0. });

  void ComputeNextPoints(std::vector<IBMesh::IBVertex> force,
                         real_t dt) override;
};

} // namespace Models
} // namespace ELFF
