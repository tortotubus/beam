#include "elff/models/rigidbody/IBPinnedRigidBodyCircle.hpp"

namespace ELFF {
namespace Models {

IBPinnedRigidBodyCircle::IBPinnedRigidBodyCircle(real_t radius,
                                                 size_t point_count,
                                                 real_t density,
                                                 const Vec3& center,
                                                 real_t angle)
  : IBRigidBodyCircle(radius,
                      point_count,
                      density,
                      center,
                      { 0., 0., 0. },
                      angle,
                      0.)
{}

void
IBPinnedRigidBodyCircle::ComputeNextPoints(
  std::vector<IBMesh::IBVertex> force,
  real_t dt)
{
  static_cast<void>(dt);
  ELFF_ASSERT(force.size() == GetNumberOfPoints(),
              "IBPinnedRigidBodyCircle::ComputeNextPoints(): force array must "
              "match point count.\n");
  CopyCurrentToNext();
}

} // namespace Models
} // namespace ELFF
