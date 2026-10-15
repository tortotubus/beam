#include "elff/models/rigidbody/RigidBodyCircle.hpp"

#include "elff/general/error.hpp"

#include <cmath>

namespace ELFF {
namespace Models {

namespace {
constexpr real_t kPi =
  static_cast<real_t>(3.141592653589793238462643383279502884L);
}

RigidBodyCircle::RigidBodyCircle(real_t radius,
                                 size_t point_count,
                                 real_t density)
  : radius_(radius)
  , point_count_(point_count)
  , density_(density)
{
  ELFF_VERIFY(radius_ > 0., "RigidBodyCircle: radius must be positive.\n");
  ELFF_VERIFY(point_count_ > 0,
              "RigidBodyCircle: point_count must be positive.\n");
  ELFF_VERIFY(density_ > 0., "RigidBodyCircle: density must be positive.\n");

  initialize();
}

void
RigidBodyCircle::define_reference_configuration(
  std::vector<Vec3>& points_ref,
  std::vector<real_t>& ds,
  Vec3& cog_ref,
  std::vector<Vec3>& normals_ref) const
{
  points_ref.resize(point_count_);
  ds.assign(point_count_,
            2. * kPi * radius_ / static_cast<real_t>(point_count_));
  normals_ref.resize(point_count_);
  cog_ref = { 0., 0., 0. };

  for (size_t i = 0; i < point_count_; ++i) {
    const real_t theta =
      2. * kPi * static_cast<real_t>(i) / static_cast<real_t>(point_count_);
    const real_t c = std::cos(theta);
    const real_t s = std::sin(theta);

    normals_ref[i] = { c, s, 0. };
    points_ref[i] = { radius_ * c, radius_ * s, 0. };
  }
}

void
RigidBodyCircle::define_mass_properties(real_t& mass, real_t& Izz) const
{
  mass = density_ * kPi * radius_ * radius_;
  Izz = 0.5 * mass * radius_ * radius_;
}

} // namespace Models
} // namespace ELFF
