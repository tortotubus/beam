#include "elff/models/rigidbody/RigidBody2DMesh.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

namespace ELFF {
namespace Models {
namespace {

using Vec3 = RigidBody2DMesh::Vec3;

TEST(RigidBody2DMeshTest, UpdateFromPoseAppliesPlanarRotationAndVelocity)
{
  RigidBody2DMesh mesh;
  mesh.set_reference_surface(std::vector<Vec3>{ { 1., 0., 0. } },
                             std::vector<real_t>{ 1. });
  mesh.set_reference_center_of_gravity({ 0., 0., 0. });

  const real_t c = std::sqrt(0.5);
  mesh.update_from_pose(
    { 1., 2., 0. }, { c, 0., 0., c }, { 0.5, 0., 0. }, { 0., 0., 2. });

  const Vec3& p = mesh.world_points()[0];
  const Vec3& v = mesh.world_velocities()[0];

  EXPECT_NEAR(p[0], 1., 1e-12);
  EXPECT_NEAR(p[1], 3., 1e-12);
  EXPECT_NEAR(p[2], 0., 1e-12);

  EXPECT_NEAR(v[0], -1.5, 1e-12);
  EXPECT_NEAR(v[1], 0., 1e-12);
  EXPECT_NEAR(v[2], 0., 1e-12);
}

TEST(RigidBody2DMeshTest, IntegratesPlanarForceAndZTorque)
{
  RigidBody2DMesh mesh;
  mesh.set_reference_surface(
    std::vector<Vec3>{ { 1., 0., 0. }, { -1., 0., 0. } },
    std::vector<real_t>{ 1., 1. });
  mesh.set_reference_center_of_gravity({ 0., 0., 0. });
  mesh.update_from_pose(
    { 0., 0., 0. }, { 1., 0., 0., 0. }, { 0., 0., 0. }, { 0., 0., 0. });

  const std::vector<Vec3> traction = { { 0., 1., 5. }, { 0., -1., -7. } };

  const Vec3 F = mesh.integrate_force(traction);
  const Vec3 T = mesh.integrate_torque_about_com(traction, { 0., 0., 0. });

  EXPECT_NEAR(F[0], 0., 1e-12);
  EXPECT_NEAR(F[1], 0., 1e-12);
  EXPECT_NEAR(F[2], 0., 1e-12);

  EXPECT_NEAR(T[0], 0., 1e-12);
  EXPECT_NEAR(T[1], 0., 1e-12);
  EXPECT_NEAR(T[2], 2., 1e-12);
}

} // namespace
} // namespace Models
} // namespace ELFF
