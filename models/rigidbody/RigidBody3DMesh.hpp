#pragma once

#include "elff/config/config.hpp"
#include "elff/general/error.hpp"

#include <array>
#include <vector>

namespace ELFF {
namespace Models {

/**
 * @brief Surface computation mesh for rigid-body traction integration.
 *
 * The reference surface is stored in body coordinates and paired with
 * per-node quadrature weights DS. The mesh can be transformed into world
 * coordinates from a rigid pose and exposes weighted integration helpers.
 */
class RigidBody3DMesh
{
public:
  using Vec3 = std::array<real_t, 3>;

  explicit RigidBody3DMesh(size_t number_of_points = 0);

  void set_reference_surface(const std::vector<Vec3>& points_ref,
                             const std::vector<real_t>& ds,
                             const std::vector<Vec3>& normals_ref = {});

  void set_reference_center_of_gravity(const Vec3& cog_ref);

  void update_from_pose(const Vec3& x_com_world,
                        const std::array<real_t, 4>& q_bw,
                        const Vec3& v_com_world,
                        const Vec3& omega_world);

  Vec3 integrate_force(const std::vector<Vec3>& traction_world) const;

  Vec3 integrate_torque_about_com(const std::vector<Vec3>& traction_world,
                                  const Vec3& x_com_world) const;

  size_t get_number_of_points() const { return points_ref_.size(); }

  const std::vector<Vec3>& reference_points() const { return points_ref_; }
  const std::vector<real_t>& ds() const { return ds_; }
  const std::vector<Vec3>& reference_normals() const { return normals_ref_; }
  const Vec3& reference_center_of_gravity() const { return cog_ref_; }

  const std::vector<Vec3>& world_points() const { return points_world_; }
  const std::vector<Vec3>& world_velocities() const { return velocity_world_; }

  real_t total_measure() const;

private:
  std::vector<Vec3> points_ref_;
  std::vector<real_t> ds_;
  std::vector<Vec3> normals_ref_;

  Vec3 cog_ref_ = {0., 0., 0.};

  std::vector<Vec3> points_world_;
  std::vector<Vec3> velocity_world_;

  static Vec3 add(const Vec3& a, const Vec3& b);
  static Vec3 sub(const Vec3& a, const Vec3& b);
  static Vec3 mul(real_t s, const Vec3& a);
  static real_t dot(const Vec3& a, const Vec3& b);
  static Vec3 cross(const Vec3& a, const Vec3& b);

  static std::array<real_t, 4> normalized_quaternion(
    const std::array<real_t, 4>& q);
  static std::array<std::array<real_t, 3>, 3> quaternion_to_rotation_bw(
    const std::array<real_t, 4>& q_bw);
  static Vec3 matvec(const std::array<std::array<real_t, 3>, 3>& A,
                     const Vec3& x);
};

} // namespace Models
} // namespace ELFF
