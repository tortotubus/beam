#pragma once

#include "elff/models/rigidbody/RigidBody3DModel.hpp"

namespace ELFF {
namespace Models {

/**
 * @brief Convenience rigid body backed directly by supplied reference data.
 */
class RigidBody3DFromReference : public RigidBody3DModel
{
public:
  RigidBody3DFromReference(const std::vector<Vec3>& points_ref,
                           const std::vector<real_t>& ds,
                           const Vec3& cog_ref,
                           real_t mass,
                           const Mat3& I_body,
                           const std::vector<Vec3>& normals_ref = {});

protected:
  void define_reference_configuration(
    std::vector<Vec3>& points_ref,
    std::vector<real_t>& ds,
    Vec3& cog_ref,
    std::vector<Vec3>& normals_ref) const override;
  void define_mass_properties(real_t& mass, Mat3& I_body) const override;

private:
  std::vector<Vec3> points_ref_;
  std::vector<real_t> ds_;
  Vec3 cog_ref_;
  real_t mass_;
  Mat3 I_body_;
  std::vector<Vec3> normals_ref_;
};

} // namespace Models
} // namespace ELFF
