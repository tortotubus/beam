#include "elff/models/rigidbody/IBPinnedRigidBodySphere.hpp"

namespace ELFF {
namespace Models {

IBPinnedRigidBodySphere::IBPinnedRigidBodySphere(real_t radius,
                                                 size_t point_count,
                                                 real_t density,
                                                 const Vec3& center,
                                                 const Quat& q_bw)
  : IBRigidBodySphere(radius,
                      point_count,
                      density,
                      center,
                      { 0., 0., 0. },
                      q_bw,
                      { 0., 0., 0. })
{}

void
IBPinnedRigidBodySphere::ComputeNextPoints(
  std::vector<IBMesh::IBVertex> force,
  real_t dt)
{
  static_cast<void>(dt);
  ELFF_ASSERT(force.size() == GetNumberOfPoints(),
              "IBPinnedRigidBodySphere::ComputeNextPoints(): force array must "
              "match point count.\n");
  CopyCurrentToNext();
}

} // namespace Models
} // namespace ELFF
