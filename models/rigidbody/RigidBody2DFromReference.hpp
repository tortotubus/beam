#pragma once

#include "elff/models/rigidbody/RigidBody2DModel.hpp"

namespace ELFF {
namespace Models {

/**
 * @brief Convenience planar rigid body backed by supplied reference data.
 */
class RigidBody2DFromReference : public RigidBody2DModel
{
public:
  RigidBody2DFromReference(const std::vector<Vec3>& points_ref,
                           const std::vector<real_t>& ds,
                           const Vec3& cog_ref,
                           real_t mass,
                           real_t Izz,
                           const std::vector<Vec3>& normals_ref = {});

protected:
  void define_reference_configuration(
    std::vector<Vec3>& points_ref,
    std::vector<real_t>& ds,
    Vec3& cog_ref,
    std::vector<Vec3>& normals_ref) const override;
  void define_mass_properties(real_t& mass, real_t& Izz) const override;

private:
  std::vector<Vec3> points_ref_;
  std::vector<real_t> ds_;
  Vec3 cog_ref_;
  real_t mass_;
  real_t Izz_;
  std::vector<Vec3> normals_ref_;
};

} // namespace Models
} // namespace ELFF
