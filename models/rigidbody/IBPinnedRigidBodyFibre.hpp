#pragma once

#include "elff/models/rigidbody/IBRigidBodyFibre.hpp"

namespace ELFF {
namespace Models {

/**
 * @brief Immersed-boundary rigid fibre with fixed marker positions.
 *
 * This is useful for fixed immersed obstacles that should expose the same
 * centerline marker geometry as @ref IBRigidBodyFibre while ignoring fluid
 * forces.
 */
class IBPinnedRigidBodyFibre : public IBRigidBodyFibre
{
public:
  using Vec3 = RigidBodyModel::Vec3;
  using Quat = IBRigidBodyFibre::Quat;

  IBPinnedRigidBodyFibre(real_t length,
                         real_t diameter,
                         size_t nodes,
                         real_t linear_density,
                         const Vec3& center = { 0., 0., 0. },
                         const Quat& q_bw = { 1., 0., 0., 0. });

  void ComputeNextPoints(std::vector<IBMesh::IBVertex> force,
                         real_t dt) override;
};

} // namespace Models
} // namespace ELFF
