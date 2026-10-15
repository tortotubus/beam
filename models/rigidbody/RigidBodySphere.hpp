#pragma once

#include "elff/models/rigidbody/RigidBody3DModel.hpp"

namespace ELFF {
namespace Models {

/**
 * @brief Solid sphere rigid body with a Fibonacci-sphere surface
 * discretization.
 *
 * Reference markers are generated on a sphere of radius \f( R \f) centered at
 * the origin \f( \mathbf{x}_{\mathrm{cog,ref}} = \mathbf{0} \f). The nodal
 * quadrature measure is uniform,
 * \f( \Delta S_i = 4\pi R^2/N \f), and normals are radial/outward.
 *
 * Mass properties correspond to a uniform-density solid sphere:
 * \f( m = \rho\,\frac{4}{3}\pi R^3 \f),
 * \f( \mathbf{I} = \frac{2}{5}mR^2\mathbf{I}_3 \f).
 */
class RigidBodySphere : public RigidBody3DModel
{
public:
  /**
   * @brief Construct a sphere model.
   * @param radius Sphere radius \f( R > 0 \f)
   * @param point_count Number of surface markers \f( N > 0 \f)
   * @param density Solid density \f( \rho > 0 \f)
   */
  RigidBodySphere(real_t radius, size_t point_count, real_t density = 1.0);

  real_t radius() const { return radius_; }
  size_t point_count() const { return point_count_; }
  real_t density() const { return density_; }

protected:
  void define_reference_configuration(
    std::vector<Vec3>& points_ref,
    std::vector<real_t>& ds,
    Vec3& cog_ref,
    std::vector<Vec3>& normals_ref) const override;

  void define_mass_properties(real_t& mass, Mat3& I_body) const override;
  Vec3 angular_momentum_rhs_body(const Vec3& tau_body,
                                 const Vec3& omega_body,
                                 const Vec3& L_body) const override;

private:
  real_t radius_;
  size_t point_count_;
  real_t density_;
};

} // namespace Models
} // namespace ELFF
