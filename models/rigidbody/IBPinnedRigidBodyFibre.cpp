#include "elff/models/rigidbody/IBPinnedRigidBodyFibre.hpp"

namespace ELFF {
namespace Models {

IBPinnedRigidBodyFibre::IBPinnedRigidBodyFibre(real_t length,
                                               real_t diameter,
                                               size_t nodes,
                                               real_t linear_density,
                                               const Vec3& center,
                                               const Quat& q_bw)
  : IBRigidBodyFibre(length,
                     diameter,
                     nodes,
                     linear_density,
                     center,
                     { 0., 0., 0. },
                     q_bw,
                     { 0., 0., 0. })
{}

void
IBPinnedRigidBodyFibre::ComputeNextPoints(std::vector<IBMesh::IBVertex> force,
                                          real_t dt)
{
  static_cast<void>(dt);
  ELFF_ASSERT(force.size() == GetNumberOfPoints(),
              "IBPinnedRigidBodyFibre::ComputeNextPoints(): force array must "
              "match node count.\n");
  CopyCurrentToNext();
}

} // namespace Models
} // namespace ELFF
