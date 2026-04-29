#include "elff/models/rigidbody/RigidBodyMesh.hpp"

namespace ELFF {
namespace Models {

RigidBodyMesh::RigidBodyMesh(size_t number_of_points)
  : points_ref_(number_of_points, { 0., 0., 0. })
  , ds_(number_of_points, 1.)
  , points_world_(number_of_points, { 0., 0., 0. })
  , velocity_world_(number_of_points, { 0., 0., 0. })
{
}

void
RigidBodyMesh::set_reference_surface(const std::vector<Vec3>& points_ref,
                                     const std::vector<real_t>& ds,
                                     const std::vector<Vec3>& normals_ref)
{
  ELFF_VERIFY(
    !points_ref.empty(),
    "RigidBodyMesh::set_reference_surface(): points cannot be empty.\n");
  ELFF_VERIFY(points_ref.size() == ds.size(),
              "RigidBodyMesh::set_reference_surface(): points and DS must have "
              "same size.\n");
  ELFF_VERIFY(
    normals_ref.empty() || normals_ref.size() == points_ref.size(),
    "RigidBodyMesh::set_reference_surface(): normals must be empty or "
    "match points size.\n");

  for (size_t i = 0; i < ds.size(); ++i) {
    ELFF_VERIFY(
      ds[i] > 0.,
      "RigidBodyMesh::set_reference_surface(): DS entries must be positive.\n");
  }

  points_ref_ = points_ref;
  ds_ = ds;
  normals_ref_ = normals_ref;
  points_world_.assign(points_ref_.size(), { 0., 0., 0. });
  velocity_world_.assign(points_ref_.size(), { 0., 0., 0. });
}

void
RigidBodyMesh::set_reference_center_of_gravity(const Vec3& cog_ref)
{
  cog_ref_ = cog_ref;
}

real_t
RigidBodyMesh::total_measure() const
{
  real_t sum = 0.;
  for (real_t q : ds_) {
    sum += q;
  }
  return sum;
}

void
RigidBodyMesh::verify_traction_size(const std::vector<Vec3>& traction_world,
                                    const char* context) const
{
  ELFF_VERIFY(traction_world.size() == points_ref_.size(), context);
}

void
RigidBodyMesh::verify_world_points_initialized(const char* context) const
{
  ELFF_VERIFY(points_world_.size() == points_ref_.size(), context);
}

RigidBodyMesh::Vec3
RigidBodyMesh::add(const Vec3& a, const Vec3& b)
{
  return { a[0] + b[0], a[1] + b[1], a[2] + b[2] };
}

RigidBodyMesh::Vec3
RigidBodyMesh::sub(const Vec3& a, const Vec3& b)
{
  return { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
}

RigidBodyMesh::Vec3
RigidBodyMesh::mul(real_t s, const Vec3& a)
{
  return { s * a[0], s * a[1], s * a[2] };
}

real_t
RigidBodyMesh::dot(const Vec3& a, const Vec3& b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

RigidBodyMesh::Vec3
RigidBodyMesh::cross(const Vec3& a, const Vec3& b)
{
  return { a[1] * b[2] - a[2] * b[1],
           a[2] * b[0] - a[0] * b[2],
           a[0] * b[1] - a[1] * b[0] };
}

} // namespace Models
} // namespace ELFF
