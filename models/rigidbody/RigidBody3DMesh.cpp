#include "elff/models/rigidbody/RigidBody3DMesh.hpp"

#include <cmath>

namespace ELFF {
namespace Models {

RigidBody3DMesh::RigidBody3DMesh(size_t number_of_points)
  : points_ref_(number_of_points, {0., 0., 0.})
  , ds_(number_of_points, 1.)
  , points_world_(number_of_points, {0., 0., 0.})
  , velocity_world_(number_of_points, {0., 0., 0.})
{
}

void
RigidBody3DMesh::set_reference_surface(const std::vector<Vec3>& points_ref,
                                     const std::vector<real_t>& ds,
                                     const std::vector<Vec3>& normals_ref)
{
  ELFF_VERIFY(!points_ref.empty(),
              "RigidBody3DMesh::set_reference_surface(): points cannot be empty.\n");
  ELFF_VERIFY(points_ref.size() == ds.size(),
              "RigidBody3DMesh::set_reference_surface(): points and DS must have same size.\n");
  ELFF_VERIFY(normals_ref.empty() || normals_ref.size() == points_ref.size(),
              "RigidBody3DMesh::set_reference_surface(): normals must be empty or match points size.\n");

  for (size_t i = 0; i < ds.size(); ++i) {
    ELFF_VERIFY(ds[i] > 0.,
                "RigidBody3DMesh::set_reference_surface(): DS entries must be positive.\n");
  }

  points_ref_ = points_ref;
  ds_ = ds;
  normals_ref_ = normals_ref;

  points_world_.assign(points_ref_.size(), {0., 0., 0.});
  velocity_world_.assign(points_ref_.size(), {0., 0., 0.});
}

void
RigidBody3DMesh::set_reference_center_of_gravity(const Vec3& cog_ref)
{
  cog_ref_ = cog_ref;
}

void
RigidBody3DMesh::update_from_pose(const Vec3& x_com_world,
                                const std::array<real_t, 4>& q_bw,
                                const Vec3& v_com_world,
                                const Vec3& omega_world)
{
  ELFF_VERIFY(points_ref_.size() == ds_.size(),
              "RigidBody3DMesh::update_from_pose(): mesh is not initialized.\n");

  const auto R_bw = quaternion_to_rotation_bw(q_bw);

  for (size_t i = 0; i < points_ref_.size(); ++i) {
    const Vec3 r_body = sub(points_ref_[i], cog_ref_);
    const Vec3 r_world = matvec(R_bw, r_body);

    points_world_[i] = add(x_com_world, r_world);
    velocity_world_[i] = add(v_com_world, cross(omega_world, r_world));
  }
}

RigidBody3DMesh::Vec3
RigidBody3DMesh::integrate_force(const std::vector<Vec3>& traction_world) const
{
  ELFF_VERIFY(traction_world.size() == points_ref_.size(),
              "RigidBody3DMesh::integrate_force(): traction size must match node count.\n");

  Vec3 F = {0., 0., 0.};
  for (size_t i = 0; i < traction_world.size(); ++i) {
    F = add(F, mul(ds_[i], traction_world[i]));
  }
  return F;
}

RigidBody3DMesh::Vec3
RigidBody3DMesh::integrate_torque_about_com(const std::vector<Vec3>& traction_world,
                                          const Vec3& x_com_world) const
{
  ELFF_VERIFY(traction_world.size() == points_ref_.size(),
              "RigidBody3DMesh::integrate_torque_about_com(): traction size must match node count.\n");
  ELFF_VERIFY(points_world_.size() == points_ref_.size(),
              "RigidBody3DMesh::integrate_torque_about_com(): world points not initialized.\n");

  Vec3 tau = {0., 0., 0.};
  for (size_t i = 0; i < traction_world.size(); ++i) {
    const Vec3 r = sub(points_world_[i], x_com_world);
    const Vec3 fi = mul(ds_[i], traction_world[i]);
    tau = add(tau, cross(r, fi));
  }
  return tau;
}

real_t
RigidBody3DMesh::total_measure() const
{
  real_t sum = 0.;
  for (real_t q : ds_)
    sum += q;
  return sum;
}

RigidBody3DMesh::Vec3
RigidBody3DMesh::add(const Vec3& a, const Vec3& b)
{
  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

RigidBody3DMesh::Vec3
RigidBody3DMesh::sub(const Vec3& a, const Vec3& b)
{
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

RigidBody3DMesh::Vec3
RigidBody3DMesh::mul(real_t s, const Vec3& a)
{
  return {s * a[0], s * a[1], s * a[2]};
}

real_t
RigidBody3DMesh::dot(const Vec3& a, const Vec3& b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

RigidBody3DMesh::Vec3
RigidBody3DMesh::cross(const Vec3& a, const Vec3& b)
{
  return {a[1] * b[2] - a[2] * b[1],
          a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}

std::array<real_t, 4>
RigidBody3DMesh::normalized_quaternion(const std::array<real_t, 4>& q)
{
  const real_t n2 = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
  ELFF_VERIFY(n2 > 0.,
              "RigidBody3DMesh::normalized_quaternion(): invalid quaternion norm.\n");
  const real_t invn = 1. / std::sqrt(n2);
  return {q[0] * invn, q[1] * invn, q[2] * invn, q[3] * invn};
}

std::array<std::array<real_t, 3>, 3>
RigidBody3DMesh::quaternion_to_rotation_bw(const std::array<real_t, 4>& q_bw)
{
  const auto q = normalized_quaternion(q_bw);
  const real_t w = q[0], x = q[1], y = q[2], z = q[3];

  return {{{1. - 2. * (y * y + z * z), 2. * (x * y - z * w),
            2. * (x * z + y * w)},
           {2. * (x * y + z * w), 1. - 2. * (x * x + z * z),
            2. * (y * z - x * w)},
           {2. * (x * z - y * w), 2. * (y * z + x * w),
            1. - 2. * (x * x + y * y)}}};
}

RigidBody3DMesh::Vec3
RigidBody3DMesh::matvec(const std::array<std::array<real_t, 3>, 3>& A,
                      const Vec3& x)
{
  return {dot(A[0], x), dot(A[1], x), dot(A[2], x)};
}

} // namespace Models
} // namespace ELFF
