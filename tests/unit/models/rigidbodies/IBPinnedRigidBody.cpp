#include "elff/models/rigidbody/IBPinnedRigidBodyCircle.hpp"
#include "elff/models/rigidbody/IBPinnedRigidBodySphere.hpp"

#include <gtest/gtest.h>

#include <vector>

namespace ELFF {
namespace Models {
namespace {

TEST(IBPinnedRigidBodyTest, CircleIgnoresForcesAndDoesNotAdvanceTime)
{
  const real_t R = 1.0;
  const size_t N = 64;
  IBPinnedRigidBodyCircle body(R, N, 2.0, { 0.5, -0.25, 0.0 }, 0.2);

  const auto initial_points = body.GetCurrent().GetPoints();
  const auto initial_velocity = body.GetCurrent().GetVelocity();

  std::vector<IBMesh::IBVertex> force(N, { 10.0, -4.0, 0.0 });
  IBMesh& next = body.GetNext(force, 0.1);

  EXPECT_NEAR(body.time(), 0.0, 1e-12);
  EXPECT_EQ(body.time_iteration(), 0u);
  ASSERT_EQ(next.GetNumberOfPoints(), N);

  const auto& next_points = next.GetPoints();
  const auto& next_velocity = next.GetVelocity();
  for (size_t i = 0; i < N; ++i) {
    EXPECT_NEAR(next_points[i].x, initial_points[i].x, 1e-12);
    EXPECT_NEAR(next_points[i].y, initial_points[i].y, 1e-12);
    EXPECT_NEAR(next_points[i].z, initial_points[i].z, 1e-12);
    EXPECT_NEAR(next_velocity[i].x, initial_velocity[i].x, 1e-12);
    EXPECT_NEAR(next_velocity[i].y, initial_velocity[i].y, 1e-12);
    EXPECT_NEAR(next_velocity[i].z, initial_velocity[i].z, 1e-12);
  }
}

TEST(IBPinnedRigidBodyTest, SphereIgnoresForcesAndDoesNotAdvanceTime)
{
  const real_t R = 0.75;
  const size_t N = 96;
  IBPinnedRigidBodySphere body(
    R, N, 1.5, { -0.1, 0.2, 0.3 }, { 0.98, 0.1, -0.05, 0.15 });

  const auto initial_points = body.GetCurrent().GetPoints();
  const auto initial_velocity = body.GetCurrent().GetVelocity();

  std::vector<IBMesh::IBVertex> force(N, { 3.0, -2.0, 1.0 });
  IBMesh& next = body.GetNext(force, 0.05);

  EXPECT_NEAR(body.time(), 0.0, 1e-12);
  EXPECT_EQ(body.time_iteration(), 0u);
  ASSERT_EQ(next.GetNumberOfPoints(), N);

  const auto& next_points = next.GetPoints();
  const auto& next_velocity = next.GetVelocity();
  for (size_t i = 0; i < N; ++i) {
    EXPECT_NEAR(next_points[i].x, initial_points[i].x, 1e-12);
    EXPECT_NEAR(next_points[i].y, initial_points[i].y, 1e-12);
    EXPECT_NEAR(next_points[i].z, initial_points[i].z, 1e-12);
    EXPECT_NEAR(next_velocity[i].x, initial_velocity[i].x, 1e-12);
    EXPECT_NEAR(next_velocity[i].y, initial_velocity[i].y, 1e-12);
    EXPECT_NEAR(next_velocity[i].z, initial_velocity[i].z, 1e-12);
  }
}

} // namespace
} // namespace Models
} // namespace ELFF
