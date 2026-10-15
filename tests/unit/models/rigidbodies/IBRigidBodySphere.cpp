#include "elff/models/rigidbody/IBRigidBodySphere.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

namespace ELFF {
namespace Models {
namespace {

using Vec3 = IBRigidBodySphere::Vec3;

TEST(IBRigidBodySphereTest, CurrentMeshMirrorsRigidBodyState)
{
  const real_t R = 1.5;
  const size_t N = 128;
  const Vec3 center = { 0.5, -0.25, 0.75 };
  const Vec3 velocity = { 0.1, -0.2, 0.3 };
  const Vec3 omega = { 0.05, -0.1, 0.2 };

  IBRigidBodySphere body(
    R, N, 1.2, center, velocity, { 1., 0., 0., 0. }, omega);

  IBMesh& ib_mesh = body.GetCurrent();
  ASSERT_EQ(ib_mesh.GetNumberOfPoints(), N);

  const auto& rb_points = body.mesh().world_points();
  const auto& rb_velocity = body.mesh().world_velocities();
  const auto& rb_ds = body.mesh().ds();
  const auto& ib_points = ib_mesh.GetPoints();
  const auto& ib_velocity = ib_mesh.GetVelocity();
  const auto& ib_measure = ib_mesh.GetMeasures();

  for (size_t i = 0; i < N; ++i) {
    EXPECT_NEAR(ib_points[i].x, rb_points[i][0], 1e-12);
    EXPECT_NEAR(ib_points[i].y, rb_points[i][1], 1e-12);
    EXPECT_NEAR(ib_points[i].z, rb_points[i][2], 1e-12);
    EXPECT_NEAR(ib_velocity[i].x, rb_velocity[i][0], 1e-12);
    EXPECT_NEAR(ib_velocity[i].y, rb_velocity[i][1], 1e-12);
    EXPECT_NEAR(ib_velocity[i].z, rb_velocity[i][2], 1e-12);
    EXPECT_NEAR(ib_measure[i], rb_ds[i], 1e-12);
  }
}

TEST(IBRigidBodySphereTest, GetNextAdvancesRigidBodyFromForceDensity)
{
  const real_t R = 1.0;
  const size_t N = 256;
  const real_t density = 2.0;
  IBRigidBodySphere body(R, N, density);

  std::vector<IBMesh::IBVertex> force(N, { 1.0, 0.0, 0.0 });
  const real_t dt = 0.05;
  IBMesh& next = body.GetNext(force, dt);

  const real_t pi = std::acos(-1.);
  const real_t mass = density * (4. / 3.) * pi * R * R * R;
  const real_t total_force = 4. * pi * R * R;
  const real_t expected_vx = dt * total_force / mass;

  EXPECT_NEAR(body.com_velocity()[0], expected_vx, 1e-10);
  EXPECT_NEAR(body.com_velocity()[1], 0.0, 1e-12);
  EXPECT_NEAR(body.com_velocity()[2], 0.0, 1e-12);
  EXPECT_NEAR(body.time(), dt, 1e-12);
  EXPECT_EQ(body.time_iteration(), 1u);

  ASSERT_EQ(next.GetNumberOfPoints(), N);
  const auto& rb_points = body.mesh().world_points();
  const auto& ib_points = next.GetPoints();
  for (size_t i = 0; i < N; ++i) {
    EXPECT_NEAR(ib_points[i].x, rb_points[i][0], 1e-12);
    EXPECT_NEAR(ib_points[i].y, rb_points[i][1], 1e-12);
    EXPECT_NEAR(ib_points[i].z, rb_points[i][2], 1e-12);
  }
}

TEST(IBRigidBodySphereTest, PackUnpackRestoresDynamicStateAndIBMesh)
{
  const real_t R = 0.75;
  const size_t N = 192;
  const real_t density = 1.4;

  IBRigidBodySphere body(R,
                         N,
                         density,
                         { 0.25, -0.1, 0.2 },
                         { 0.2, 0.1, -0.05 },
                         { 0.98, 0.1, -0.05, 0.15 },
                         { 0.05, 0.03, -0.02 });

  std::vector<IBMesh::IBVertex> force(N, { 0.5, 0.25, -0.1 });
  body.GetNext(force, 0.05);

  IBModelState state;
  body.pack_state(state);

  IBRigidBodySphere restored(R, N, density);
  restored.unpack_state(state);

  for (size_t i = 0; i < 3; ++i) {
    EXPECT_NEAR(restored.center_of_mass()[i], body.center_of_mass()[i], 1e-12);
    EXPECT_NEAR(restored.com_velocity()[i], body.com_velocity()[i], 1e-12);
    EXPECT_NEAR(restored.angular_velocity_body()[i],
                body.angular_velocity_body()[i],
                1e-12);
  }
  for (size_t i = 0; i < 4; ++i) {
    EXPECT_NEAR(restored.pose()[i], body.pose()[i], 1e-12);
  }
  EXPECT_NEAR(restored.time(), body.time(), 1e-12);
  EXPECT_EQ(restored.time_iteration(), body.time_iteration());

  const auto& restored_points = restored.GetCurrent().GetPoints();
  const auto& body_points = body.GetCurrent().GetPoints();
  ASSERT_EQ(restored_points.size(), body_points.size());
  for (size_t i = 0; i < restored_points.size(); ++i) {
    EXPECT_NEAR(restored_points[i].x, body_points[i].x, 1e-12);
    EXPECT_NEAR(restored_points[i].y, body_points[i].y, 1e-12);
    EXPECT_NEAR(restored_points[i].z, body_points[i].z, 1e-12);
  }
}

} // namespace
} // namespace Models
} // namespace ELFF
