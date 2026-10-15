#include "elff/models/rigidbody/RigidBody2DFromReference.hpp"

namespace ELFF {
namespace Models {

RigidBody2DFromReference::RigidBody2DFromReference(
  const std::vector<Vec3>& points_ref,
  const std::vector<real_t>& ds,
  const Vec3& cog_ref,
  real_t mass,
  real_t Izz,
  const std::vector<Vec3>& normals_ref)
  : points_ref_(points_ref)
  , ds_(ds)
  , cog_ref_(cog_ref)
  , mass_(mass)
  , Izz_(Izz)
  , normals_ref_(normals_ref)
{
  initialize();
}

void
RigidBody2DFromReference::define_reference_configuration(
  std::vector<Vec3>& points_ref,
  std::vector<real_t>& ds,
  Vec3& cog_ref,
  std::vector<Vec3>& normals_ref) const
{
  points_ref = points_ref_;
  ds = ds_;
  cog_ref = cog_ref_;
  normals_ref = normals_ref_;
}

void
RigidBody2DFromReference::define_mass_properties(real_t& mass, real_t& Izz) const
{
  mass = mass_;
  Izz = Izz_;
}

} // namespace Models
} // namespace ELFF
