#include "elff/models/rigidbody/RigidBody3DModel.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

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

  real_t norm = 0;
  for (auto q : q_bw)
  {
    ELFF_VERIFY(std::isfinite(q), "Invalid quaternion.");
    norm += q * q;
  }
  ELFF_VERIFY(std::isfinite(norm) && norm > 0, "Invalid quaternion norm.");
  const auto momentum = matvec(I_body_, omega_body);
  for (size_t d = 0; d < 3; ++d)
    ELFF_VERIFY(std::isfinite(x_com_world[d]) && std::isfinite(v_com_world[d]) &&
                    std::isfinite(omega_body[d]) && std::isfinite(momentum[d]),
                "Invalid rigid-body state.");
  set_center_of_mass(x_com_world);
  set_com_velocity(v_com_world);
  q_bw_ = normalize_q(q_bw);
  omega_body_ = omega_body;
  L_body_ = momentum;
  reset_time_integration_state();
  update_derived_state();
}

RigidBody3DModel::Vec3 RigidBody3DModel::world_point(const Vec3& point) const
{
  return add(center_of_mass(), matvec(q_to_Rbw(q_bw_), point));
}

RigidBody3DModel::Vec3 RigidBody3DModel::angular_velocity_world() const
{
  return matvec(q_to_Rbw(q_bw_), omega_body_);
}

RigidBody3DModel::Vec3 RigidBody3DModel::velocity_at_world(const Vec3& point) const
{
  return add(com_velocity(), cross(angular_velocity_world(), sub(point, center_of_mass())));
}

void RigidBody3DModel::verify_axisymmetric_inertia() const
{
  const real_t scale = std::max({I_body_[0][0], I_body_[1][1], I_body_[2][2]});
  const real_t tol = 32 * std::numeric_limits<real_t>::epsilon() * scale;
  ELFF_VERIFY(I_body_[0][0] > 0 && I_body_[1][1] > 0 && I_body_[2][2] > 0 && std::isfinite(scale) &&
                  std::abs(I_body_[1][1] - I_body_[2][2]) <= tol,
              "Contact stepping requires axisymmetric inertia about body x.");
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      if (i != j)
        ELFF_VERIFY(std::isfinite(I_body_[i][j]) && std::abs(I_body_[i][j]) <= tol,
                    "Contact stepping requires diagonal inertia.");
}

void RigidBody3DModel::kick(real_t dt, const Contact::Wrench& wrench)
{
  verify_axisymmetric_inertia();
  ELFF_VERIFY(std::isfinite(dt) && dt > 0, "Kick dt must be finite and positive.");
  for (size_t d = 0; d < 3; ++d)
    ELFF_VERIFY(std::isfinite(wrench.force[d]) && std::isfinite(wrench.torque[d]),
                "Invalid wrench.");
  const auto velocity = add(com_velocity(), mul(dt / mass(), wrench.force));
  const auto torque_body = matvec(transpose(q_to_Rbw(q_bw_)), wrench.torque);
  const auto momentum = add(L_body_, mul(dt, torque_body));
  const auto omega = matvec(I_body_inv_, momentum);
  for (size_t d = 0; d < 3; ++d)
    ELFF_VERIFY(std::isfinite(velocity[d]) && std::isfinite(momentum[d]) && std::isfinite(omega[d]),
                "Kick state overflow.");
  set_com_velocity(velocity);
  L_body_ = momentum;
  omega_body_ = omega;
  update_derived_state();
}

void RigidBody3DModel::drift(real_t dt)
{
  verify_axisymmetric_inertia();
  ELFF_VERIFY(std::isfinite(dt) && dt > 0 && std::isfinite(time() + dt) && time() + dt > time() &&
                  time_iteration() < std::numeric_limits<size_t>::max(),
              "Invalid drift interval.");
  const auto position = add(center_of_mass(), mul(dt, com_velocity()));
  for (auto x : position)
    ELFF_VERIFY(std::isfinite(x), "Drift position overflow.");
  // Exact free symmetric-top motion: world precession followed by body-x spin.
  const auto worldL = matvec(q_to_Rbw(q_bw_), L_body_);
  auto rotation = [&](const Vec3& angular)
  {
    const real_t speed = std::hypot(angular[0], angular[1], angular[2]);
    if (speed == 0)
      return std::array<real_t, 4>{1, 0, 0, 0};
    const real_t angle = speed * dt / 2;
    ELFF_VERIFY(std::isfinite(angle), "Rotation overflow.");
    const real_t factor = std::sin(angle) / speed;
    return std::array<real_t, 4>{std::cos(angle), factor * angular[0], factor * angular[1],
                                 factor * angular[2]};
  };
  const auto precession = rotation(mul(real_t{1} / I_body_[1][1], worldL));
  const auto spin =
      rotation({(real_t{1} / I_body_[0][0] - real_t{1} / I_body_[1][1]) * L_body_[0], 0, 0});
  q_bw_ = normalize_q(qmul(qmul(precession, q_bw_), spin));
  L_body_ = matvec(transpose(q_to_Rbw(q_bw_)), worldL);
  omega_body_ = matvec(I_body_inv_, L_body_);
  for (size_t d = 0; d < 3; ++d)
    ELFF_VERIFY(std::isfinite(L_body_[d]) && std::isfinite(omega_body_[d]),
                "Drift angular state overflow.");
  set_center_of_mass(position);
  set_time_integration_state(time() + dt, time_iteration() + 1);
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
  // Normalize before taking a determinant: dimensional absolute tolerances
  // incorrectly reject small spheres and overflow for large bodies.
  real_t scale = 0;
  for (const auto& row : A)
    for (auto value : row) {
      ELFF_VERIFY(std::isfinite(value), "Inertia tensor must be finite.");
      scale = std::max(scale, std::abs(value));
    }
  ELFF_VERIFY(scale > 0, "Rigid-body inertia tensor is singular.");
  Mat3 scaled = A;
  for (auto& row : scaled)
    for (auto& value : row) value /= scale;
  const real_t d = det3(scaled);
  ELFF_VERIFY(std::isfinite(d) && d != 0, "Rigid-body inertia tensor is singular.");
  const real_t invd = (real_t{1} / scale) / d;

  Mat3 B = { { { 0., 0., 0. }, { 0., 0., 0. }, { 0., 0., 0. } } };

  B[0][0] = (scaled[1][1] * scaled[2][2] - scaled[1][2] * scaled[2][1]) * invd;
  B[0][1] = (scaled[0][2] * scaled[2][1] - scaled[0][1] * scaled[2][2]) * invd;
  B[0][2] = (scaled[0][1] * scaled[1][2] - scaled[0][2] * scaled[1][1]) * invd;

  B[1][0] = (scaled[1][2] * scaled[2][0] - scaled[1][0] * scaled[2][2]) * invd;
  B[1][1] = (scaled[0][0] * scaled[2][2] - scaled[0][2] * scaled[2][0]) * invd;
  B[1][2] = (scaled[0][2] * scaled[1][0] - scaled[0][0] * scaled[1][2]) * invd;

  B[2][0] = (scaled[1][0] * scaled[2][1] - scaled[1][1] * scaled[2][0]) * invd;
  B[2][1] = (scaled[0][1] * scaled[2][0] - scaled[0][0] * scaled[2][1]) * invd;
  B[2][2] = (scaled[0][0] * scaled[1][1] - scaled[0][1] * scaled[1][0]) * invd;

  for (const auto& row : B)
    for (auto value : row)
      ELFF_VERIFY(std::isfinite(value), "Inverse inertia is unrepresentable.");

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
