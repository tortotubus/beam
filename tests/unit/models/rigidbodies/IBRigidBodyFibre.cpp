#include "elff/models/rigidbody/IBRigidBodyFibre.hpp"

#include <gtest/gtest.h>

#include <vector>

namespace ELFF {
namespace Models {
namespace {

using Vec3 = IBRigidBodyFibre::Vec3;

TEST(IBRigidBodyFibreTest, CurrentMeshMirrorsRigidBodyState)
{
  const real_t L = 2.0;
  const real_t d = 0.1;
  const size_t N = 7;
  const Vec3 center = { 0.5, -0.25, 0.75 };
  const Vec3 velocity = { 0.1, -0.2, 0.3 };
  const Vec3 omega = { 0.05, -0.1, 0.2 };

  IBRigidBodyFibre fibre(
    L, d, N, 1.2, center, velocity, { 1., 0., 0., 0. }, omega);

  IBMesh& ib_mesh = fibre.GetCurrent();
  ASSERT_EQ(ib_mesh.GetNumberOfPoints(), N);

  const auto& rb_points = fibre.mesh().world_points();
  const auto& rb_velocity = fibre.mesh().world_velocities();
  const auto& rb_ds = fibre.mesh().ds();
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

TEST(IBRigidBodyFibreTest, GetNextAdvancesRigidBodyFromForceDensity)
{
  const real_t L = 1.5;
  const size_t N = 9;
  const real_t linear_density = 2.0;
  IBRigidBodyFibre fibre(L, 0.1, N, linear_density);

  std::vector<IBMesh::IBVertex> force(N, { 1.0, 0.0, 0.0 });
  const real_t dt = 0.05;
  IBMesh& next = fibre.GetNext(force, dt);

  const real_t expected_vx = dt * L / (linear_density * L);
  EXPECT_NEAR(fibre.com_velocity()[0], expected_vx, 1e-12);
  EXPECT_NEAR(fibre.com_velocity()[1], 0.0, 1e-12);
  EXPECT_NEAR(fibre.com_velocity()[2], 0.0, 1e-12);
  EXPECT_NEAR(fibre.time(), dt, 1e-12);
  EXPECT_EQ(fibre.time_iteration(), 1u);

  ASSERT_EQ(next.GetNumberOfPoints(), N);
  const auto& rb_points = fibre.mesh().world_points();
  const auto& ib_points = next.GetPoints();
  for (size_t i = 0; i < N; ++i) {
    EXPECT_NEAR(ib_points[i].x, rb_points[i][0], 1e-12);
    EXPECT_NEAR(ib_points[i].y, rb_points[i][1], 1e-12);
    EXPECT_NEAR(ib_points[i].z, rb_points[i][2], 1e-12);
  }
}

TEST(IBRigidBodyFibreTest, PackUnpackRestoresDynamicStateAndIBMesh)
{
  const real_t L = 1.25;
  const real_t d = 0.08;
  const size_t N = 8;
  const real_t linear_density = 1.4;

  IBRigidBodyFibre fibre(L,
                         d,
                         N,
                         linear_density,
                         { 0.25, -0.1, 0.2 },
                         { 0.2, 0.1, -0.05 },
                         { 0.98, 0.1, -0.05, 0.15 },
                         { 0.05, 0.03, -0.02 });

  std::vector<IBMesh::IBVertex> force(N, { 0.5, 0.25, -0.1 });
  fibre.GetNext(force, 0.05);

  IBModelState state;
  fibre.pack_state(state);

  IBRigidBodyFibre restored(L, d, N, linear_density);
  restored.unpack_state(state);

  for (size_t i = 0; i < 3; ++i) {
    EXPECT_NEAR(restored.center_of_mass()[i], fibre.center_of_mass()[i], 1e-12);
    EXPECT_NEAR(restored.com_velocity()[i], fibre.com_velocity()[i], 1e-12);
    EXPECT_NEAR(restored.angular_velocity_body()[i],
                fibre.angular_velocity_body()[i],
                1e-12);
  }
  for (size_t i = 0; i < 4; ++i) {
    EXPECT_NEAR(restored.pose()[i], fibre.pose()[i], 1e-12);
  }
  EXPECT_NEAR(restored.time(), fibre.time(), 1e-12);
  EXPECT_EQ(restored.time_iteration(), fibre.time_iteration());

  const auto& restored_points = restored.GetCurrent().GetPoints();
  const auto& fibre_points = fibre.GetCurrent().GetPoints();
  ASSERT_EQ(restored_points.size(), fibre_points.size());
  for (size_t i = 0; i < restored_points.size(); ++i) {
    EXPECT_NEAR(restored_points[i].x, fibre_points[i].x, 1e-12);
    EXPECT_NEAR(restored_points[i].y, fibre_points[i].y, 1e-12);
    EXPECT_NEAR(restored_points[i].z, fibre_points[i].z, 1e-12);
  }
}

} // namespace
} // namespace Models
} // namespace ELFF
