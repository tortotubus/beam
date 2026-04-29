#pragma once

#include "elff/config/config.hpp"
#include "elff/general/error.hpp"
#include "elff/models/rigidbody/RigidBodyMesh.hpp"

#include <array>
#include <vector>

namespace ELFF {
namespace Models {

/**
 * @brief Surface computation mesh for rigid-body traction integration.
 *
 * The reference surface is stored in body coordinates and paired with
 * per-node quadrature weights DS. The mesh can be transformed into world
 * coordinates from a rigid pose and exposes weighted integration helpers.
 */
class RigidBody3DMesh : public RigidBodyMesh
{
public:
  using Vec3 = RigidBodyMesh::Vec3;

  explicit RigidBody3DMesh(size_t number_of_points = 0);

  void update_from_pose(const Vec3& x_com_world,
                        const std::array<real_t, 4>& q_bw,
                        const Vec3& v_com_world,
                        const Vec3& omega_world) override;

  Vec3 integrate_force(const std::vector<Vec3>& traction_world) const override;

  Vec3 integrate_torque_about_com(const std::vector<Vec3>& traction_world,
                                  const Vec3& x_com_world) const override;

private:
  static std::array<real_t, 4> normalized_quaternion(
    const std::array<real_t, 4>& q);
  static std::array<std::array<real_t, 3>, 3> quaternion_to_rotation_bw(
    const std::array<real_t, 4>& q_bw);
  static Vec3 matvec(const std::array<std::array<real_t, 3>, 3>& A,
                     const Vec3& x);
};

} // namespace Models
} // namespace ELFF
