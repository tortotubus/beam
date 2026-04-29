#pragma once

#include "elff/models/rigidbody/RigidBodyMesh.hpp"

namespace ELFF {
namespace Models {

/**
 * @brief Placeholder 2D rigid-body mesh API.
 *
 * TODO: Define 2D marker storage, quadrature weights, and world/reference
 * transform utilities.
 */
class RigidBody2DMesh : public RigidBodyMesh
{
public:
  explicit RigidBody2DMesh(size_t number_of_points = 0)
    : RigidBodyMesh(number_of_points)
  {
  }

  void update_from_pose(const Vec3& x_com_world,
                        const std::array<real_t, 4>& q_bw,
                        const Vec3& v_com_world,
                        const Vec3& omega_world) override;

  Vec3 integrate_force(const std::vector<Vec3>& traction_world) const override;
  Vec3 integrate_torque_about_com(const std::vector<Vec3>& traction_world,
                                  const Vec3& x_com_world) const override;
};

} // namespace Models
} // namespace ELFF
