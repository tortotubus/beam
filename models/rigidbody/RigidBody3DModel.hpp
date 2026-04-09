#pragma once

#include "elff/config/config.hpp"
#include "elff/general/error.hpp"
#include "elff/models/rigidbody/RigidBody3DMesh.hpp"

#include <array>
#include <vector>

namespace ELFF {
namespace Models {

class RigidBody3DModel
{
public:
  using Vec3 = std::array<real_t, 3>;
  using Mat3 = std::array<std::array<real_t, 3>, 3>;

  virtual ~RigidBody3DModel() = default;

  void initialize();
  void step(real_t dt, const std::vector<Vec3>& traction_world);

  const RigidBody3DMesh& mesh() const { return mesh_; }
  RigidBody3DMesh& mesh() { return mesh_; }

  const Vec3& center_of_mass() const { return x_com_world_; }
  const Vec3& com_velocity() const { return v_com_world_; }
  const std::array<real_t, 4>& pose() const { return q_bw_; }
  const Vec3& angular_velocity_body() const { return omega_body_; }
  const Vec3& angular_momentum_body() const { return L_body_; }

  real_t time() const { return t_; }
  size_t time_iteration() const { return time_iter_; }

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

private:
  RigidBody3DMesh mesh_;

  real_t mass_ = 1.;
  Mat3 I_body_ = { { { 1., 0., 0. }, { 0., 1., 0. }, { 0., 0., 1. } } };
  Mat3 I_body_inv_ = { { { 1., 0., 0. }, { 0., 1., 0. }, { 0., 0., 1. } } };

  Vec3 x_com_world_ = { 0., 0., 0. };
  Vec3 v_com_world_ = { 0., 0., 0. };
  std::array<real_t, 4> q_bw_ = { 1., 0., 0., 0. };

  Vec3 omega_body_ = { 0., 0., 0. };
  Vec3 L_body_ = { 0., 0., 0. };

  real_t t_ = 0.;
  size_t time_iter_ = 0;
  bool initialized_ = false;

  static Vec3 add(const Vec3& a, const Vec3& b);
  static Vec3 sub(const Vec3& a, const Vec3& b);
  static Vec3 mul(real_t s, const Vec3& a);
  static real_t dot(const Vec3& a, const Vec3& b);
  static Vec3 cross(const Vec3& a, const Vec3& b);

  static Mat3 transpose(const Mat3& A);
  static Vec3 matvec(const Mat3& A, const Vec3& x);
  static real_t det3(const Mat3& A);
  static Mat3 inv3(const Mat3& A);

  static std::array<real_t, 4> normalize_q(const std::array<real_t, 4>& q);
  static std::array<real_t, 4> qmul(const std::array<real_t, 4>& a,
                                    const std::array<real_t, 4>& b);
  static Mat3 q_to_Rbw(const std::array<real_t, 4>& q_bw);

  void update_derived_state();
};

} // namespace Models
} // namespace ELFF
