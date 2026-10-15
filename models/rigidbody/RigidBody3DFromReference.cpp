#include "elff/models/rigidbody/RigidBody3DFromReference.hpp"

namespace ELFF {
namespace Models {

RigidBody3DFromReference::RigidBody3DFromReference(
  const std::vector<Vec3>& points_ref,
  const std::vector<real_t>& ds,
  const Vec3& cog_ref,
  real_t mass,
  const Mat3& I_body,
  const std::vector<Vec3>& normals_ref)
  : points_ref_(points_ref)
  , ds_(ds)
  , cog_ref_(cog_ref)
  , mass_(mass)
  , I_body_(I_body)
  , normals_ref_(normals_ref)
{
  initialize();
}

void
RigidBody3DFromReference::define_reference_configuration(
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
RigidBody3DFromReference::define_mass_properties(real_t& mass,
                                                 Mat3& I_body) const
{
  mass = mass_;
  I_body = I_body_;
}

} // namespace Models
} // namespace ELFF
