#include "elff/models/rigidbody/RigidBodyModel.hpp"

namespace ELFF {
namespace Models {

void
RigidBodyModel::initialize()
{
  initialize_model();
  reset_time_integration_state();
  mark_initialized();
  update_derived_state();
}

void
RigidBodyModel::step(real_t dt, const std::vector<Vec3>& traction_world)
{
  verify_initialized_for_step(
    "RigidBodyModel::step(): model must be initialized before stepping.\n");
  verify_positive_dt(dt, "RigidBodyModel::step(): dt must be positive.\n");
  ELFF_VERIFY(
    traction_world.size() == expected_traction_count(),
    "RigidBodyModel::step(): traction size must match mesh point count.\n");

  const Vec3 F_world = integrate_force_world(traction_world);
  const Vec3 tau_world =
    integrate_torque_world_about_com(traction_world, x_com_world_);

  const Vec3 a_world = mul(1. / mass_, F_world);
  v_com_world_ = add(v_com_world_, mul(dt, a_world));
  x_com_world_ = add(x_com_world_, mul(dt, v_com_world_));

  integrate_rotation(dt, tau_world);
  update_derived_state();

  advance_time(dt);
}

void
RigidBodyModel::set_mass(real_t mass)
{
  ELFF_VERIFY(mass > 0.,
              "RigidBodyModel::set_mass(): mass must be positive.\n");
  mass_ = mass;
}

RigidBodyModel::Vec3
RigidBodyModel::add(const Vec3& a, const Vec3& b)
{
  return { a[0] + b[0], a[1] + b[1], a[2] + b[2] };
}

RigidBodyModel::Vec3
RigidBodyModel::sub(const Vec3& a, const Vec3& b)
{
  return { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
}

RigidBodyModel::Vec3
RigidBodyModel::mul(real_t s, const Vec3& a)
{
  return { s * a[0], s * a[1], s * a[2] };
}

real_t
RigidBodyModel::dot(const Vec3& a, const Vec3& b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

RigidBodyModel::Vec3
RigidBodyModel::cross(const Vec3& a, const Vec3& b)
{
  return { a[1] * b[2] - a[2] * b[1],
           a[2] * b[0] - a[0] * b[2],
           a[0] * b[1] - a[1] * b[0] };
}

void
RigidBodyModel::reset_time_integration_state()
{
  t_ = 0.;
  time_iter_ = 0;
}

void
RigidBodyModel::set_time_integration_state(real_t t, size_t time_iter)
{
  ELFF_VERIFY(t >= 0.,
              "RigidBodyModel::set_time_integration_state(): time must be "
              "nonnegative.\n");
  t_ = t;
  time_iter_ = time_iter;
}

void
RigidBodyModel::mark_initialized()
{
  initialized_ = true;
}

void
RigidBodyModel::verify_initialized_for_step(const char* context) const
{
  ELFF_VERIFY(initialized_, context);
}

void
RigidBodyModel::verify_positive_dt(real_t dt, const char* context) const
{
  ELFF_VERIFY(dt > 0., context);
}

void
RigidBodyModel::advance_time(real_t dt)
{
  t_ += dt;
  ++time_iter_;
}

} // namespace Models
} // namespace ELFF
