#pragma once

#include "elff/config/config.hpp"
#include "elff/general/error.hpp"
#include "elff/models/rigidbody/RigidBody2DMesh.hpp"
#include "elff/models/rigidbody/RigidBodyModel.hpp"

#include <array>
#include <vector>

namespace ELFF {
namespace Models {

/**
 * @brief Planar rigid-body model (stored in R3 with redundant z components).
 */
class RigidBody2DModel : public RigidBodyModel
{
public:
  virtual ~RigidBody2DModel() = default;

  const RigidBody2DMesh& mesh() const { return mesh_; }
  RigidBody2DMesh& mesh() { return mesh_; }

  real_t angle() const { return theta_; }
  real_t angular_velocity_z() const { return omega_z_; }
  real_t angular_momentum_z() const { return Lz_; }
  const std::array<real_t, 4>& pose() const { return q_bw_; }

protected:
  using Vec3 = RigidBodyModel::Vec3;

  RigidBody2DModel() = default;

  virtual void define_reference_configuration(
    std::vector<Vec3>& points_ref,
    std::vector<real_t>& ds,
    Vec3& cog_ref,
    std::vector<Vec3>& normals_ref) const = 0;

  virtual void define_mass_properties(real_t& mass, real_t& Izz) const = 0;

  void initialize_from_reference(const std::vector<Vec3>& points_ref,
                                 const std::vector<real_t>& ds,
                                 const Vec3& cog_ref,
                                 real_t mass,
                                 real_t Izz,
                                 const std::vector<Vec3>& normals_ref = {});

  void set_initial_state(const Vec3& x_com_world,
                         const Vec3& v_com_world,
                         real_t theta,
                         real_t omega_z);

  void initialize_model() override;
  size_t expected_traction_count() const override;
  Vec3 integrate_force_world(
    const std::vector<Vec3>& traction_world) const override;
  Vec3 integrate_torque_world_about_com(const std::vector<Vec3>& traction_world,
                                        const Vec3& x_com_world) const override;
  void integrate_rotation(real_t dt, const Vec3& tau_world) override;
  void update_derived_state() override;

  RigidBody2DMesh mesh_;

  real_t Izz_ = 1.;
  real_t Izz_inv_ = 1.;
  real_t theta_ = 0.;
  real_t omega_z_ = 0.;
  real_t Lz_ = 0.;
  std::array<real_t, 4> q_bw_ = { 1., 0., 0., 0. };
};

} // namespace Models
} // namespace ELFF
