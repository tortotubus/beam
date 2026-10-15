#pragma once

#include "elff/models/rigidbody/RigidBody2DModel.hpp"

namespace ELFF {
namespace Models {

/**
 * @brief Solid circular rigid body with uniform perimeter discretization.
 *
 * Reference markers are generated on a circle of radius \f( R \f) in the
 * \f$(x,y)\f$ plane and centered at the origin. Quadrature uses uniform arc
 * measure \f( \Delta s_i = 2\pi R/N \f), with outward radial normals.
 *
 * Mass properties correspond to a uniform-density solid disk of unit thickness:
 * \f( m = \rho\,\pi R^2 \f), \f( I_{zz} = \frac{1}{2} m R^2 \f).
 */
class RigidBodyCircle : public RigidBody2DModel
{
public:
  /**
   * @brief Construct a circular rigid body model.
   * @param radius Circle radius \f( R > 0 \f)
   * @param point_count Number of perimeter markers \f( N > 0 \f)
   * @param density Disk density \f( \rho > 0 \f)
   */
  RigidBodyCircle(real_t radius, size_t point_count, real_t density = 1.0);

  real_t radius() const { return radius_; }
  size_t point_count() const { return point_count_; }
  real_t density() const { return density_; }

protected:
  void define_reference_configuration(
    std::vector<Vec3>& points_ref,
    std::vector<real_t>& ds,
    Vec3& cog_ref,
    std::vector<Vec3>& normals_ref) const override;

  void define_mass_properties(real_t& mass, real_t& Izz) const override;

private:
  real_t radius_;
  size_t point_count_;
  real_t density_;
};

} // namespace Models
} // namespace ELFF
