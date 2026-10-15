#include "elff/models/rigidbody/RigidBodyFibre.hpp"

#include "elff/general/error.hpp"

namespace ELFF {
namespace Models {

namespace {
inline real_t
square(real_t x)
{
  return x * x;
}

RigidBodyFibre::Mat3
cylinder_inertia(real_t length, real_t diameter, real_t mass)
{
  const real_t radius = static_cast<real_t>(0.5) * diameter;
  const real_t Ixx = static_cast<real_t>(0.5) * mass * square(radius);
  const real_t Iyy =
    mass * (static_cast<real_t>(3.) * square(radius) + square(length)) /
    static_cast<real_t>(12.);

  return { { { Ixx, 0., 0. }, { 0., Iyy, 0. }, { 0., 0., Iyy } } };
}
} // namespace

RigidBodyFibre::RigidBodyFibre(real_t length,
                               real_t diameter,
                               size_t nodes,
                               real_t linear_density)
  : length_(length)
  , diameter_(diameter)
  , nodes_(nodes)
  , linear_density_(linear_density)
  , total_mass_(linear_density * length)
  , inertia_body_(cylinder_inertia(length, diameter, total_mass_))
{
  ELFF_VERIFY(length_ > 0., "RigidBodyFibre: length must be positive.\n");
  ELFF_VERIFY(diameter_ > 0., "RigidBodyFibre: diameter must be positive.\n");
  ELFF_VERIFY(nodes_ >= 2, "RigidBodyFibre: nodes must be at least 2.\n");
  ELFF_VERIFY(linear_density_ > 0.,
              "RigidBodyFibre: linear_density must be positive.\n");

  initialize();
}

void
RigidBodyFibre::define_reference_configuration(
  std::vector<Vec3>& points_ref,
  std::vector<real_t>& ds,
  Vec3& cog_ref,
  std::vector<Vec3>& normals_ref) const
{
  points_ref.resize(nodes_);
  ds.assign(nodes_, length_ / static_cast<real_t>(nodes_ - 1));
  ds.front() *= static_cast<real_t>(0.5);
  ds.back() *= static_cast<real_t>(0.5);
  cog_ref = { 0., 0., 0. };
  normals_ref.clear();

  for (size_t i = 0; i < nodes_; ++i) {
    const real_t s = static_cast<real_t>(i) /
                     static_cast<real_t>(nodes_ - 1);
    points_ref[i] = { -static_cast<real_t>(0.5) * length_ + s * length_,
                      0.,
                      0. };
  }
}

void
RigidBodyFibre::define_mass_properties(real_t& mass, Mat3& I_body) const
{
  mass = total_mass_;
  I_body = inertia_body_;
}

} // namespace Models
} // namespace ELFF
