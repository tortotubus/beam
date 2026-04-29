#include "elff/models/rigidbody/RigidBody2DModel.hpp"

#include <cmath>

namespace ELFF {
namespace Models {

void
RigidBody2DModel::initialize_model()
{
  std::vector<Vec3> points_ref;
  std::vector<real_t> ds;
  Vec3 cog_ref = { 0., 0., 0. };
  std::vector<Vec3> normals_ref;
  real_t mass = 0.;
  real_t Izz = 0.;

  define_reference_configuration(points_ref, ds, cog_ref, normals_ref);
  define_mass_properties(mass, Izz);
  initialize_from_reference(points_ref, ds, cog_ref, mass, Izz, normals_ref);
}

void
RigidBody2DModel::initialize_from_reference(
  const std::vector<Vec3>& points_ref,
  const std::vector<real_t>& ds,
  const Vec3& cog_ref,
  real_t mass,
  real_t Izz,
  const std::vector<Vec3>& normals_ref)
{
  ELFF_VERIFY(
    Izz > 0.,
    "RigidBody2DModel::initialize_from_reference(): Izz must be positive.\n");

  mesh_.set_reference_surface(points_ref, ds, normals_ref);
  mesh_.set_reference_center_of_gravity(cog_ref);

  set_mass(mass);
  Izz_ = Izz;
  Izz_inv_ = 1. / Izz_;

  set_center_of_mass({ 0., 0., 0. });
  set_com_velocity({ 0., 0., 0. });
  theta_ = 0.;
  omega_z_ = 0.;
  Lz_ = 0.;
  q_bw_ = { 1., 0., 0., 0. };
}

void
RigidBody2DModel::set_initial_state(const Vec3& x_com_world,
                                    const Vec3& v_com_world,
                                    real_t theta,
                                    real_t omega_z)
{
  ELFF_VERIFY(is_initialized(),
              "RigidBody2DModel::set_initial_state(): model must be "
              "initialized first.\n");

  set_center_of_mass({ x_com_world[0], x_com_world[1], 0. });
  set_com_velocity({ v_com_world[0], v_com_world[1], 0. });
  theta_ = theta;
  omega_z_ = omega_z;
  Lz_ = Izz_ * omega_z_;
  reset_time_integration_state();
  update_derived_state();
}

size_t
RigidBody2DModel::expected_traction_count() const
{
  return mesh_.get_number_of_points();
}

RigidBody2DModel::Vec3
RigidBody2DModel::integrate_force_world(
  const std::vector<Vec3>& traction_world) const
{
  return mesh_.integrate_force(traction_world);
}

RigidBody2DModel::Vec3
RigidBody2DModel::integrate_torque_world_about_com(
  const std::vector<Vec3>& traction_world,
  const Vec3& x_com_world) const
{
  return mesh_.integrate_torque_about_com(traction_world, x_com_world);
}

void
RigidBody2DModel::integrate_rotation(real_t dt, const Vec3& tau_world)
{
  const real_t alpha_z = Izz_inv_ * tau_world[2];
  omega_z_ += dt * alpha_z;
  Lz_ = Izz_ * omega_z_;
  theta_ += dt * omega_z_;
}

void
RigidBody2DModel::update_derived_state()
{
  const real_t half = 0.5 * theta_;
  const real_t c = std::cos(half);
  const real_t s = std::sin(half);
  q_bw_ = { c, 0., 0., s };
  mesh_.update_from_pose(
    center_of_mass(), q_bw_, com_velocity(), { 0., 0., omega_z_ });
}

} // namespace Models
} // namespace ELFF
