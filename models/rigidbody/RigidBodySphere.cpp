#include "elff/models/rigidbody/RigidBodySphere.hpp"

#include "elff/general/error.hpp"

#include <algorithm>
#include <cmath>

namespace ELFF {
namespace Models {

namespace {
constexpr real_t kPi =
  static_cast<real_t>(3.141592653589793238462643383279502884L);

inline real_t
square(real_t x)
{
  return x * x;
}
} // namespace

RigidBodySphere::RigidBodySphere(real_t radius,
                                 size_t point_count,
                                 real_t density)
  : radius_(radius)
  , point_count_(point_count)
  , density_(density)
{
  ELFF_VERIFY(radius_ > 0., "RigidBodySphere: radius must be positive.\n");
  ELFF_VERIFY(point_count_ > 0,
              "RigidBodySphere: point_count must be positive.\n");
  ELFF_VERIFY(density_ > 0., "RigidBodySphere: density must be positive.\n");

  initialize();
}

void
RigidBodySphere::define_reference_configuration(
  std::vector<Vec3>& points_ref,
  std::vector<real_t>& ds,
  Vec3& cog_ref,
  std::vector<Vec3>& normals_ref) const
{
  points_ref.resize(point_count_);
  ds.assign(point_count_,
            4. * kPi * square(radius_) / static_cast<real_t>(point_count_));
  normals_ref.resize(point_count_);
  cog_ref = { 0., 0., 0. };

  // Golden angle in radians.
  const real_t golden_angle = static_cast<real_t>(kPi * (3. - std::sqrt(5.0)));

  for (size_t i = 0; i < point_count_; ++i) {
    const real_t y =
      1. - 2. * (static_cast<real_t>(i) + static_cast<real_t>(0.5)) /
             static_cast<real_t>(point_count_);
    const real_t r_xy =
      std::sqrt(std::max(static_cast<real_t>(0.), 1. - y * y));
    const real_t theta = golden_angle * static_cast<real_t>(i);

    const real_t x = r_xy * std::cos(theta);
    const real_t z = r_xy * std::sin(theta);

    normals_ref[i] = { x, y, z };
    points_ref[i] = { radius_ * x, radius_ * y, radius_ * z };
  }
}

void
RigidBodySphere::define_mass_properties(real_t& mass, Mat3& I_body) const
{
  mass = density_ * (4. / 3.) * kPi * radius_ * radius_ * radius_;
  const real_t inertia = (2. / 5.) * mass * radius_ * radius_;

  I_body = {
    { { inertia, 0., 0. }, { 0., inertia, 0. }, { 0., 0., inertia } }
  };
}

RigidBodySphere::Vec3
RigidBodySphere::angular_momentum_rhs_body(const Vec3& tau_body,
                                           const Vec3& /*omega_body*/,
                                           const Vec3& /*L_body*/) const
{
  // Sphere specialization: disable rotational coupling term (omega x L).
  return tau_body;
}

} // namespace Models
} // namespace ELFF
