#pragma once

#include "elff/models/rigidbody/IBRigidBodyCircle.hpp"

namespace ELFF {
namespace Models {

/**
 * @brief Immersed-boundary circular rigid body with fixed marker positions.
 *
 * This is useful for fixed immersed obstacles that should expose the same
 * marker geometry as @ref IBRigidBodyCircle while ignoring fluid forces.
 */
class IBPinnedRigidBodyCircle : public IBRigidBodyCircle
{
public:
  using Vec3 = RigidBodyModel::Vec3;

  IBPinnedRigidBodyCircle(real_t radius,
                          size_t point_count,
                          real_t density = 1.0,
                          const Vec3& center = { 0., 0., 0. },
                          real_t angle = 0.);

  void ComputeNextPoints(std::vector<IBMesh::IBVertex> force,
                         real_t dt) override;
};

} // namespace Models
} // namespace ELFF
