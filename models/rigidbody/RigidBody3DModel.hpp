#pragma once

#include "elff/config/config.hpp"
#include "elff/general/error.hpp"
#include "elff/models/rigidbody/RigidBody3DMesh.hpp"
#include "elff/models/rigidbody/RigidBodyModel.hpp"

#include <array>
#include <vector>

namespace ELFF {
namespace Models {

class RigidBody3DModel : public RigidBodyModel
{
public:
  using Vec3 = RigidBodyModel::Vec3;
  using Mat3 = std::array<std::array<real_t, 3>, 3>;

  virtual ~RigidBody3DModel() = default;

  const RigidBody3DMesh& mesh() const { return mesh_; }
  RigidBody3DMesh& mesh() { return mesh_; }
  const std::array<real_t, 4>& pose() const { return q_bw_; }
  const Vec3& angular_velocity_body() const { return omega_body_; }
  const Vec3& angular_momentum_body() const { return L_body_; }

protected:
  RigidBody3DModel() = default;

  virtual void define_reference_configuration(
    std::vector<Vec3>& points_ref,
    std::vector<real_t>& ds,
    Vec3& cog_ref,
    std::vector<Vec3>& normals_ref) const = 0;

  virtual void define_mass_properties(real_t& mass, Mat3& I_body) const = 0;

  void initialize_from_reference(const std::vector<Vec3>& points_ref,
                                 const std::vector<real_t>& ds,
                                 const Vec3& cog_ref,
                                 real_t mass,
                                 const Mat3& I_body,
                                 const std::vector<Vec3>& normals_ref = {});

  void set_initial_state(const Vec3& x_com_world,
                         const Vec3& v_com_world,
                         const std::array<real_t, 4>& q_bw,
                         const Vec3& omega_body);

  void initialize_model() override;
  size_t expected_traction_count() const override;
  Vec3 integrate_force_world(
    const std::vector<Vec3>& traction_world) const override;
  Vec3 integrate_torque_world_about_com(const std::vector<Vec3>& traction_world,
                                        const Vec3& x_com_world) const override;
  void integrate_rotation(real_t dt, const Vec3& tau_world) override;
  void update_derived_state() override;
  virtual Vec3 angular_momentum_rhs_body(const Vec3& tau_body,
                                         const Vec3& omega_body,
                                         const Vec3& L_body) const;

private:
  RigidBody3DMesh mesh_;
  Mat3 I_body_ = { { { 1., 0., 0. }, { 0., 1., 0. }, { 0., 0., 1. } } };
  Mat3 I_body_inv_ = { { { 1., 0., 0. }, { 0., 1., 0. }, { 0., 0., 1. } } };
  std::array<real_t, 4> q_bw_ = { 1., 0., 0., 0. };

  Vec3 omega_body_ = { 0., 0., 0. };
  Vec3 L_body_ = { 0., 0., 0. };

  static Mat3 transpose(const Mat3& A);
  static Vec3 matvec(const Mat3& A, const Vec3& x);
  static real_t det3(const Mat3& A);
  static Mat3 inv3(const Mat3& A);

  static std::array<real_t, 4> normalize_q(const std::array<real_t, 4>& q);
  static std::array<real_t, 4> qmul(const std::array<real_t, 4>& a,
                                    const std::array<real_t, 4>& b);
  static Mat3 q_to_Rbw(const std::array<real_t, 4>& q_bw);
};

} // namespace Models
} // namespace ELFF
