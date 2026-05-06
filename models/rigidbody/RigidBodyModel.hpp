#pragma once

#include "elff/config/config.hpp"
#include "elff/general/error.hpp"

#include <array>
#include <vector>

namespace ELFF {
namespace Models {

/**
 * @brief Dimension-agnostic rigid-body integration shell.
 *
 * Derived classes provide geometry-specific force/torque integration and
 * rotational dynamics closures.
 */
class RigidBodyModel
{
public:
  using Vec3 = std::array<real_t, 3>;

  virtual ~RigidBodyModel() = default;

  void initialize();
  void step(real_t dt, const std::vector<Vec3>& traction_world);

  const Vec3& center_of_mass() const { return x_com_world_; }
  const Vec3& com_velocity() const { return v_com_world_; }
  real_t time() const { return t_; }
  size_t time_iteration() const { return time_iter_; }
  bool is_initialized() const { return initialized_; }

protected:
  RigidBodyModel() = default;

  virtual void initialize_model() = 0;
  virtual size_t expected_traction_count() const = 0;
  virtual Vec3 integrate_force_world(
    const std::vector<Vec3>& traction_world) const = 0;
  virtual Vec3 integrate_torque_world_about_com(
    const std::vector<Vec3>& traction_world,
    const Vec3& x_com_world) const = 0;
  virtual void integrate_rotation(real_t dt, const Vec3& tau_world) = 0;
  virtual void update_derived_state() = 0;

  void reset_time_integration_state();
  void set_time_integration_state(real_t t, size_t time_iter);
  void set_mass(real_t mass);
  real_t mass() const { return mass_; }

  void set_center_of_mass(const Vec3& x_com_world)
  {
    x_com_world_ = x_com_world;
  }
  void set_com_velocity(const Vec3& v_com_world) { v_com_world_ = v_com_world; }

  static Vec3 add(const Vec3& a, const Vec3& b);
  static Vec3 sub(const Vec3& a, const Vec3& b);
  static Vec3 mul(real_t s, const Vec3& a);
  static real_t dot(const Vec3& a, const Vec3& b);
  static Vec3 cross(const Vec3& a, const Vec3& b);

private:
  void mark_initialized();
  void verify_initialized_for_step(const char* context) const;
  void verify_positive_dt(real_t dt, const char* context) const;
  void advance_time(real_t dt);

  real_t mass_ = 1.;
  Vec3 x_com_world_ = { 0., 0., 0. };
  Vec3 v_com_world_ = { 0., 0., 0. };

  real_t t_ = 0.;
  size_t time_iter_ = 0;
  bool initialized_ = false;
};

} // namespace Models
} // namespace ELFF
