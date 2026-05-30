#pragma once

#include "elff/models/rigidbody/RigidBody3DModel.hpp"

namespace ELFF {
namespace Models {

/**
 * @brief Centerline-only rigid fibre with solid-cylinder mass properties.
 *
 * The immersed-boundary markers lie on the body-frame x-axis and use
 * trapezoidal line quadrature. The supplied diameter is used only for mass
 * properties so that the 3D inertia tensor is nonsingular.
 */
class RigidBodyFibre : public RigidBody3DModel
{
public:
  RigidBodyFibre(real_t length,
                 real_t diameter,
                 size_t nodes,
                 real_t linear_density);

  real_t length() const { return length_; }
  real_t diameter() const { return diameter_; }
  size_t nodes() const { return nodes_; }
  real_t linear_density() const { return linear_density_; }
  real_t total_mass() const { return total_mass_; }
  const Mat3& inertia_body() const { return inertia_body_; }

protected:
  void define_reference_configuration(
    std::vector<Vec3>& points_ref,
    std::vector<real_t>& ds,
    Vec3& cog_ref,
    std::vector<Vec3>& normals_ref) const override;

  void define_mass_properties(real_t& mass, Mat3& I_body) const override;

private:
  real_t length_;
  real_t diameter_;
  size_t nodes_;
  real_t linear_density_;
  real_t total_mass_;
  Mat3 inertia_body_;
};

} // namespace Models
} // namespace ELFF
