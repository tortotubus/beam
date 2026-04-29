#include "elff/models/rigidbody/RigidBody3DModel.hpp"

#include <cmath>

namespace ELFF {
namespace Models {

void
RigidBody3DModel::initialize_model()
{
  std::vector<Vec3> points_ref;
  std::vector<real_t> ds;
  Vec3 cog_ref = { 0., 0., 0. };
  std::vector<Vec3> normals_ref;
  real_t mass = 0.;
  Mat3 I_body = { { { 0., 0., 0. }, { 0., 0., 0. }, { 0., 0., 0. } } };

  define_reference_configuration(points_ref, ds, cog_ref, normals_ref);
  define_mass_properties(mass, I_body);
  initialize_from_reference(points_ref, ds, cog_ref, mass, I_body, normals_ref);
}

void
RigidBody3DModel::initialize_from_reference(
  const std::vector<Vec3>& points_ref,
  const std::vector<real_t>& ds,
  const Vec3& cog_ref,
  real_t mass,
  const Mat3& I_body,
  const std::vector<Vec3>& normals_ref)
{
  mesh_.set_reference_surface(points_ref, ds, normals_ref);
  mesh_.set_reference_center_of_gravity(cog_ref);

  set_mass(mass);
  I_body_ = I_body;
  I_body_inv_ = inv3(I_body_);

  set_center_of_mass({ 0., 0., 0. });
  set_com_velocity({ 0., 0., 0. });
  q_bw_ = { 1., 0., 0., 0. };
  omega_body_ = { 0., 0., 0. };
  L_body_ = { 0., 0., 0. };
}

void
RigidBody3DModel::set_initial_state(const Vec3& x_com_world,
                                    const Vec3& v_com_world,
                                    const std::array<real_t, 4>& q_bw,
                                    const Vec3& omega_body)
{
  ELFF_VERIFY(is_initialized(),
              "RigidBody3DModel::set_initial_state(): model must be "
              "initialized first.\n");

  set_center_of_mass(x_com_world);
  set_com_velocity(v_com_world);
  q_bw_ = normalize_q(q_bw);
  omega_body_ = omega_body;
  L_body_ = matvec(I_body_, omega_body_);
  reset_time_integration_state();
  update_derived_state();
}

void
RigidBody3DModel::integrate_rotation(real_t dt, const Vec3& tau_world)
{
  const Mat3 R_bw = q_to_Rbw(q_bw_);
  const Mat3 R_wb = transpose(R_bw);
  const Vec3 tau_body = matvec(R_wb, tau_world);
  const Vec3 Ldot_body =
    angular_momentum_rhs_body(tau_body, omega_body_, L_body_);

  L_body_ = add(L_body_, mul(dt, Ldot_body));
  omega_body_ = matvec(I_body_inv_, L_body_);

  const std::array<real_t, 4> omega_q = {
    0., omega_body_[0], omega_body_[1], omega_body_[2]
  };
  const auto qdot = qmul(q_bw_, omega_q);
  q_bw_[0] += 0.5 * dt * qdot[0];
  q_bw_[1] += 0.5 * dt * qdot[1];
  q_bw_[2] += 0.5 * dt * qdot[2];
  q_bw_[3] += 0.5 * dt * qdot[3];
  q_bw_ = normalize_q(q_bw_);
}

RigidBody3DModel::Vec3
RigidBody3DModel::angular_momentum_rhs_body(const Vec3& tau_body,
                                            const Vec3& omega_body,
                                            const Vec3& L_body) const
{
  const Vec3 gyroscopic = cross(omega_body, L_body);
  return sub(tau_body, gyroscopic);
}

void
RigidBody3DModel::update_derived_state()
{
  const Mat3 R_bw = q_to_Rbw(q_bw_);
  const Vec3 omega_world = matvec(R_bw, omega_body_);
  mesh_.update_from_pose(center_of_mass(), q_bw_, com_velocity(), omega_world);
}

size_t
RigidBody3DModel::expected_traction_count() const
{
  return mesh_.get_number_of_points();
}

RigidBody3DModel::Vec3
RigidBody3DModel::integrate_force_world(
  const std::vector<Vec3>& traction_world) const
{
  return mesh_.integrate_force(traction_world);
}

RigidBody3DModel::Vec3
RigidBody3DModel::integrate_torque_world_about_com(
  const std::vector<Vec3>& traction_world,
  const Vec3& x_com_world) const
{
  return mesh_.integrate_torque_about_com(traction_world, x_com_world);
}

RigidBody3DModel::Mat3
RigidBody3DModel::transpose(const Mat3& A)
{
  return { { { A[0][0], A[1][0], A[2][0] },
             { A[0][1], A[1][1], A[2][1] },
             { A[0][2], A[1][2], A[2][2] } } };
}

RigidBody3DModel::Vec3
RigidBody3DModel::matvec(const Mat3& A, const Vec3& x)
{
  return { dot(A[0], x), dot(A[1], x), dot(A[2], x) };
}

real_t
RigidBody3DModel::det3(const Mat3& A)
{
  return A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
         A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
         A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
}

RigidBody3DModel::Mat3
RigidBody3DModel::inv3(const Mat3& A)
{
  const real_t d = det3(A);
  ELFF_VERIFY(std::abs(d) > 1e-15,
              "RigidBody3DModel::inv3(): inertia tensor is singular.\n");

  const real_t invd = 1. / d;

  Mat3 B = { { { 0., 0., 0. }, { 0., 0., 0. }, { 0., 0., 0. } } };

  B[0][0] = (A[1][1] * A[2][2] - A[1][2] * A[2][1]) * invd;
  B[0][1] = (A[0][2] * A[2][1] - A[0][1] * A[2][2]) * invd;
  B[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) * invd;

  B[1][0] = (A[1][2] * A[2][0] - A[1][0] * A[2][2]) * invd;
  B[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) * invd;
  B[1][2] = (A[0][2] * A[1][0] - A[0][0] * A[1][2]) * invd;

  B[2][0] = (A[1][0] * A[2][1] - A[1][1] * A[2][0]) * invd;
  B[2][1] = (A[0][1] * A[2][0] - A[0][0] * A[2][1]) * invd;
  B[2][2] = (A[0][0] * A[1][1] - A[0][1] * A[1][0]) * invd;

  return B;
}

std::array<real_t, 4>
RigidBody3DModel::normalize_q(const std::array<real_t, 4>& q)
{
  const real_t n2 = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
  ELFF_VERIFY(n2 > 0.,
              "RigidBody3DModel::normalize_q(): invalid quaternion.\n");
  const real_t invn = 1. / std::sqrt(n2);
  return { q[0] * invn, q[1] * invn, q[2] * invn, q[3] * invn };
}

std::array<real_t, 4>
RigidBody3DModel::qmul(const std::array<real_t, 4>& a,
                       const std::array<real_t, 4>& b)
{
  return { a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
           a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
           a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
           a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0] };
}

RigidBody3DModel::Mat3
RigidBody3DModel::q_to_Rbw(const std::array<real_t, 4>& q_bw)
{
  const auto q = normalize_q(q_bw);
  const real_t w = q[0], x = q[1], y = q[2], z = q[3];

  return {
    { { 1. - 2. * (y * y + z * z), 2. * (x * y - z * w), 2. * (x * z + y * w) },
      { 2. * (x * y + z * w), 1. - 2. * (x * x + z * z), 2. * (y * z - x * w) },
      { 2. * (x * z - y * w),
        2. * (y * z + x * w),
        1. - 2. * (x * x + y * y) } }
  };
}

} // namespace Models
} // namespace ELFF
