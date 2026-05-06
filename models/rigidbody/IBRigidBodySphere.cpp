#include "elff/models/rigidbody/IBRigidBodySphere.hpp"

#include <cmath>

namespace ELFF {
namespace Models {

namespace {
constexpr int64_t state_version = 1;
constexpr real_t state_tol = static_cast<real_t>(1e-12);

inline bool
same_real(real_t a, real_t b)
{
  return std::abs(a - b) <= state_tol * (1. + std::abs(a) + std::abs(b));
}
} // namespace

IBRigidBodySphere::IBRigidBodySphere(real_t radius,
                                     size_t point_count,
                                     real_t density,
                                     const Vec3& center,
                                     const Vec3& velocity,
                                     const Quat& q_bw,
                                     const Vec3& angular_velocity_body)
  : IBForceCoupled(point_count)
  , RigidBodySphere(radius, point_count, density)
{
  SetNodalMeasures(RigidBodySphere::mesh().ds());
  set_initial_state(center, velocity, q_bw, angular_velocity_body);
  RigidBodyToIBMeshes();
}

void
IBRigidBodySphere::ComputeNextPoints(std::vector<IBMesh::IBVertex> force,
                                     real_t dt)
{
  ELFF_ASSERT(force.size() == point_count(),
              "IBRigidBodySphere::ComputeNextPoints(): force array must match "
              "point count.\n");

  std::vector<Vec3> traction(point_count());
  for (size_t i = 0; i < force.size(); ++i) {
    traction[i] = { force[i].x, force[i].y, force[i].z };
  }

  step(dt, traction);
  RigidBodyToIBMeshes();
}

void
IBRigidBodySphere::RigidBodyToIBMesh(IBMesh& ib_mesh) const
{
  const auto& points = RigidBodySphere::mesh().world_points();
  const auto& velocities = RigidBodySphere::mesh().world_velocities();
  const auto& measures = RigidBodySphere::mesh().ds();

  ELFF_ASSERT(points.size() == ib_mesh.GetNumberOfPoints(),
              "IBRigidBodySphere::RigidBodyToIBMesh(): point count mismatch.\n");

  auto& ib_points = ib_mesh.GetPoints();
  auto& ib_velocities = ib_mesh.GetVelocity();
  auto& ib_measures = ib_mesh.GetMeasures();

  for (size_t i = 0; i < points.size(); ++i) {
    ib_points[i] = { points[i][0], points[i][1], points[i][2] };
    ib_velocities[i] = { velocities[i][0], velocities[i][1], velocities[i][2] };
    ib_measures[i] = measures[i];
  }
}

void
IBRigidBodySphere::RigidBodyToIBMeshes()
{
  RigidBodyToIBMesh(IBForceCoupled::mesh);
  RigidBodyToIBMesh(IBForceCoupled::mesh_next);
}

void
IBRigidBodySphere::pack_state(IBModelState& s) const
{
  s.ints.clear();
  s.reals.clear();
  s.bytes.clear();

  s.ints.push_back(state_version);
  s.ints.push_back(static_cast<int64_t>(point_count()));
  s.ints.push_back(static_cast<int64_t>(time_iteration()));

  const Vec3& x = center_of_mass();
  const Vec3& v = com_velocity();
  const Quat& q = pose();
  const Vec3& omega = angular_velocity_body();

  s.reals.reserve(16);
  s.reals.push_back(radius());
  s.reals.push_back(density());
  s.reals.push_back(time());
  s.reals.push_back(x[0]);
  s.reals.push_back(x[1]);
  s.reals.push_back(x[2]);
  s.reals.push_back(v[0]);
  s.reals.push_back(v[1]);
  s.reals.push_back(v[2]);
  s.reals.push_back(q[0]);
  s.reals.push_back(q[1]);
  s.reals.push_back(q[2]);
  s.reals.push_back(q[3]);
  s.reals.push_back(omega[0]);
  s.reals.push_back(omega[1]);
  s.reals.push_back(omega[2]);
}

void
IBRigidBodySphere::unpack_state(const IBModelState& s)
{
  ELFF_ASSERT(s.ints.size() == 3,
              "IBRigidBodySphere::unpack_state(): invalid integer metadata.\n");
  ELFF_ASSERT(s.ints[0] == state_version,
              "IBRigidBodySphere::unpack_state(): unsupported state version.\n");
  ELFF_ASSERT(static_cast<size_t>(s.ints[1]) == point_count(),
              "IBRigidBodySphere::unpack_state(): point count mismatch.\n");
  ELFF_ASSERT(s.reals.size() == 16,
              "IBRigidBodySphere::unpack_state(): invalid real buffer size.\n");
  ELFF_ASSERT(same_real(s.reals[0], radius()),
              "IBRigidBodySphere::unpack_state(): radius mismatch.\n");
  ELFF_ASSERT(same_real(s.reals[1], density()),
              "IBRigidBodySphere::unpack_state(): density mismatch.\n");

  Vec3 center = { s.reals[3], s.reals[4], s.reals[5] };
  Vec3 velocity = { s.reals[6], s.reals[7], s.reals[8] };
  Quat q = { s.reals[9], s.reals[10], s.reals[11], s.reals[12] };
  Vec3 omega = { s.reals[13], s.reals[14], s.reals[15] };

  set_initial_state(center, velocity, q, omega);
  set_time_integration_state(s.reals[2], static_cast<size_t>(s.ints[2]));
  RigidBodyToIBMeshes();
}

} // namespace Models
} // namespace ELFF
