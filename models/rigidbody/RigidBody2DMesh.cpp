#include "elff/models/rigidbody/RigidBody2DMesh.hpp"

namespace ELFF {
namespace Models {

void
RigidBody2DMesh::update_from_pose(const Vec3& x_com_world,
                                  const std::array<real_t, 4>& q_bw,
                                  const Vec3& v_com_world,
                                  const Vec3& omega_world)
{
  ELFF_VERIFY(
    points_ref_.size() == ds_.size(),
    "RigidBody2DMesh::update_from_pose(): mesh is not initialized.\n");

  const real_t w = q_bw[0];
  const real_t z = q_bw[3];
  const real_t n2 = w * w + z * z;
  ELFF_VERIFY(
    n2 > 0.,
    "RigidBody2DMesh::update_from_pose(): invalid planar quaternion.\n");

  // Planar rotation about z extracted from quaternion (w, z) pair.
  const real_t c = (w * w - z * z) / n2;
  const real_t s = (2. * w * z) / n2;
  const real_t omega_z = omega_world[2];

  for (size_t i = 0; i < points_ref_.size(); ++i) {
    const Vec3 r_body = sub(points_ref_[i], cog_ref_);
    const Vec3 r_world = { c * r_body[0] - s * r_body[1],
                           s * r_body[0] + c * r_body[1],
                           r_body[2] };

    points_world_[i] = add(x_com_world, r_world);
    velocity_world_[i] =
      add(v_com_world, Vec3{ -omega_z * r_world[1], omega_z * r_world[0], 0. });
  }
}

RigidBody2DMesh::Vec3
RigidBody2DMesh::integrate_force(const std::vector<Vec3>& traction_world) const
{
  verify_traction_size(traction_world,
                       "RigidBody2DMesh::integrate_force(): traction size must "
                       "match node count.\n");

  Vec3 F = { 0., 0., 0. };
  for (size_t i = 0; i < traction_world.size(); ++i) {
    // In planar mode, the third component is redundant and ignored.
    F[0] += ds_[i] * traction_world[i][0];
    F[1] += ds_[i] * traction_world[i][1];
  }
  return F;
}

RigidBody2DMesh::Vec3
RigidBody2DMesh::integrate_torque_about_com(
  const std::vector<Vec3>& traction_world,
  const Vec3& x_com_world) const
{
  verify_traction_size(
    traction_world,
    "RigidBody2DMesh::integrate_torque_about_com(): traction "
    "size must match node count.\n");
  verify_world_points_initialized(
    "RigidBody2DMesh::integrate_torque_about_com(): "
    "world points not initialized.\n");

  real_t tau_z = 0.;
  for (size_t i = 0; i < traction_world.size(); ++i) {
    const Vec3 r = sub(points_world_[i], x_com_world);
    const real_t fx = ds_[i] * traction_world[i][0];
    const real_t fy = ds_[i] * traction_world[i][1];
    tau_z += r[0] * fy - r[1] * fx;
  }

  return { 0., 0., tau_z };
}

} // namespace Models
} // namespace ELFF
