#include "elff/c/models/ibm/IBForceCoupled.h"
#include "elff/c/models/rigidbody/IBRigidBody.h"

#include <gtest/gtest.h>

#include <cmath>
#include <cstdlib>

namespace {

TEST(IBRigidBodyCAPITest, CircleConstructsAndAdvancesThroughForceCoupledAPI)
{
  const double R = 1.0;
  const int N = 96;
  const double density = 2.0;
  const vertex_t center = { 0.25, -0.5, 0.0 };
  const vertex_t velocity = { 0.0, 0.0, 0.0 };

  ib_rigid_body_circle_t circle =
    ib_rigid_body_circle_new(R, N, density, center, velocity, 0.0, 0.0);
  ASSERT_NE(circle, nullptr);

  ib_force_coupled_t model = reinterpret_cast<ib_force_coupled_t>(circle);
  EXPECT_EQ(ib_force_coupled_get_number_of_nodes(model), N);

  ib_mesh_t current = ib_force_coupled_get_current(model);
  ASSERT_EQ(current.n, N);
  EXPECT_NEAR(current.position[0].x, center.x + R, 1e-12);
  EXPECT_NEAR(current.position[0].y, center.y, 1e-12);

  vertex_t* force =
    static_cast<vertex_t*>(std::calloc(static_cast<size_t>(N), sizeof(vertex_t)));
  ASSERT_NE(force, nullptr);
  for (int i = 0; i < N; ++i) {
    force[i].x = 1.0;
  }

  const double dt = 0.1;
  ib_mesh_t next = ib_force_coupled_get_next(model, force, N, dt);
  ASSERT_EQ(next.n, N);

  const double pi = std::acos(-1.0);
  const double mass = density * pi * R * R;
  const double total_force = 2.0 * pi * R;
  const double expected_vx = dt * total_force / mass;
  EXPECT_NEAR(next.velocity[0].x, expected_vx, 1e-10);

  ib_mesh_free(&current);
  ib_mesh_free(&next);
  std::free(force);
  ib_rigid_body_circle_destroy(circle);
}

TEST(IBRigidBodyCAPITest, SphereConstructsAndAdvancesThroughForceCoupledAPI)
{
  const double R = 0.75;
  const int N = 160;
  const double density = 1.5;
  const vertex_t center = { -0.25, 0.1, 0.2 };
  const vertex_t velocity = { 0.0, 0.0, 0.0 };
  const quaternion_t q = { 1.0, 0.0, 0.0, 0.0 };
  const vertex_t omega = { 0.0, 0.0, 0.0 };

  ib_rigid_body_sphere_t sphere =
    ib_rigid_body_sphere_new(R, N, density, center, velocity, q, omega);
  ASSERT_NE(sphere, nullptr);

  ib_force_coupled_t model = reinterpret_cast<ib_force_coupled_t>(sphere);
  EXPECT_EQ(ib_force_coupled_get_number_of_nodes(model), N);

  ib_mesh_t current = ib_force_coupled_get_current(model);
  ASSERT_EQ(current.n, N);

  vertex_t* force =
    static_cast<vertex_t*>(std::calloc(static_cast<size_t>(N), sizeof(vertex_t)));
  ASSERT_NE(force, nullptr);
  for (int i = 0; i < N; ++i) {
    force[i].x = 1.0;
  }

  const double dt = 0.05;
  ib_mesh_t next = ib_force_coupled_get_next(model, force, N, dt);
  ASSERT_EQ(next.n, N);

  const double pi = std::acos(-1.0);
  const double mass = density * (4.0 / 3.0) * pi * R * R * R;
  const double total_force = 4.0 * pi * R * R;
  const double expected_vx = dt * total_force / mass;
  EXPECT_NEAR(next.velocity[0].x, expected_vx, 1e-10);

  ib_mesh_free(&current);
  ib_mesh_free(&next);
  std::free(force);
  ib_rigid_body_sphere_destroy(sphere);
}

TEST(IBRigidBodyCAPITest, FibreConstructsAndAdvancesThroughForceCoupledAPI)
{
  const double L = 1.5;
  const double d = 0.1;
  const int N = 9;
  const double linear_density = 2.0;
  const vertex_t center = { 0.25, -0.5, 0.0 };
  const vertex_t velocity = { 0.0, 0.0, 0.0 };
  const quaternion_t q = { 1.0, 0.0, 0.0, 0.0 };
  const vertex_t omega = { 0.0, 0.0, 0.0 };

  ib_rigid_body_fibre_t fibre = ib_rigid_body_fibre_new(
    L, d, N, linear_density, center, velocity, q, omega);
  ASSERT_NE(fibre, nullptr);

  ib_force_coupled_t model = reinterpret_cast<ib_force_coupled_t>(fibre);
  EXPECT_EQ(ib_force_coupled_get_number_of_nodes(model), N);

  ib_mesh_t current = ib_force_coupled_get_current(model);
  ASSERT_EQ(current.n, N);
  EXPECT_NEAR(current.position[0].x, center.x - L / 2.0, 1e-12);
  EXPECT_NEAR(current.position[N - 1].x, center.x + L / 2.0, 1e-12);

  vertex_t* force =
    static_cast<vertex_t*>(std::calloc(static_cast<size_t>(N), sizeof(vertex_t)));
  ASSERT_NE(force, nullptr);
  for (int i = 0; i < N; ++i) {
    force[i].x = 1.0;
  }

  const double dt = 0.05;
  ib_mesh_t next = ib_force_coupled_get_next(model, force, N, dt);
  ASSERT_EQ(next.n, N);

  const double expected_vx = dt * L / (linear_density * L);
  EXPECT_NEAR(next.velocity[0].x, expected_vx, 1e-12);

  ib_mesh_free(&current);
  ib_mesh_free(&next);
  std::free(force);
  ib_rigid_body_fibre_destroy(fibre);
}

TEST(IBRigidBodyCAPITest, PinnedCircleIgnoresForcesThroughForceCoupledAPI)
{
  const double R = 1.0;
  const int N = 80;
  const vertex_t center = { 0.1, -0.2, 0.0 };

  ib_pinned_rigid_body_circle_t circle =
    ib_pinned_rigid_body_circle_new(R, N, 1.0, center, 0.0);
  ASSERT_NE(circle, nullptr);

  ib_force_coupled_t model = reinterpret_cast<ib_force_coupled_t>(circle);
  ib_mesh_t current = ib_force_coupled_get_current(model);
  ASSERT_EQ(current.n, N);

  vertex_t* force =
    static_cast<vertex_t*>(std::calloc(static_cast<size_t>(N), sizeof(vertex_t)));
  ASSERT_NE(force, nullptr);
  for (int i = 0; i < N; ++i) {
    force[i].x = 100.0;
    force[i].y = -50.0;
  }

  ib_mesh_t next = ib_force_coupled_get_next(model, force, N, 0.2);
  ASSERT_EQ(next.n, N);
  for (int i = 0; i < N; ++i) {
    EXPECT_NEAR(next.position[i].x, current.position[i].x, 1e-12);
    EXPECT_NEAR(next.position[i].y, current.position[i].y, 1e-12);
    EXPECT_NEAR(next.position[i].z, current.position[i].z, 1e-12);
    EXPECT_NEAR(next.velocity[i].x, 0.0, 1e-12);
    EXPECT_NEAR(next.velocity[i].y, 0.0, 1e-12);
    EXPECT_NEAR(next.velocity[i].z, 0.0, 1e-12);
  }

  ib_mesh_free(&current);
  ib_mesh_free(&next);
  std::free(force);
  ib_pinned_rigid_body_circle_destroy(circle);
}

TEST(IBRigidBodyCAPITest, PinnedFibreIgnoresForcesThroughForceCoupledAPI)
{
  const double L = 1.0;
  const double d = 0.05;
  const int N = 7;
  const vertex_t center = { 0.1, -0.2, 0.0 };
  const quaternion_t q = { 1.0, 0.0, 0.0, 0.0 };

  ib_pinned_rigid_body_fibre_t fibre =
    ib_pinned_rigid_body_fibre_new(L, d, N, 1.0, center, q);
  ASSERT_NE(fibre, nullptr);

  ib_force_coupled_t model = reinterpret_cast<ib_force_coupled_t>(fibre);
  ib_mesh_t current = ib_force_coupled_get_current(model);
  ASSERT_EQ(current.n, N);

  vertex_t* force =
    static_cast<vertex_t*>(std::calloc(static_cast<size_t>(N), sizeof(vertex_t)));
  ASSERT_NE(force, nullptr);
  for (int i = 0; i < N; ++i) {
    force[i].x = 100.0;
    force[i].y = -50.0;
    force[i].z = 25.0;
  }

  ib_mesh_t next = ib_force_coupled_get_next(model, force, N, 0.2);
  ASSERT_EQ(next.n, N);
  for (int i = 0; i < N; ++i) {
    EXPECT_NEAR(next.position[i].x, current.position[i].x, 1e-12);
    EXPECT_NEAR(next.position[i].y, current.position[i].y, 1e-12);
    EXPECT_NEAR(next.position[i].z, current.position[i].z, 1e-12);
    EXPECT_NEAR(next.velocity[i].x, 0.0, 1e-12);
    EXPECT_NEAR(next.velocity[i].y, 0.0, 1e-12);
    EXPECT_NEAR(next.velocity[i].z, 0.0, 1e-12);
  }

  ib_mesh_free(&current);
  ib_mesh_free(&next);
  std::free(force);
  ib_pinned_rigid_body_fibre_destroy(fibre);
}

TEST(IBRigidBodyCAPITest, PinnedSphereIgnoresForcesThroughForceCoupledAPI)
{
  const double R = 0.5;
  const int N = 100;
  const vertex_t center = { 0.2, 0.1, -0.3 };
  const quaternion_t q = { 1.0, 0.0, 0.0, 0.0 };

  ib_pinned_rigid_body_sphere_t sphere =
    ib_pinned_rigid_body_sphere_new(R, N, 1.0, center, q);
  ASSERT_NE(sphere, nullptr);

  ib_force_coupled_t model = reinterpret_cast<ib_force_coupled_t>(sphere);
  ib_mesh_t current = ib_force_coupled_get_current(model);
  ASSERT_EQ(current.n, N);

  vertex_t* force =
    static_cast<vertex_t*>(std::calloc(static_cast<size_t>(N), sizeof(vertex_t)));
  ASSERT_NE(force, nullptr);
  for (int i = 0; i < N; ++i) {
    force[i].x = 100.0;
    force[i].y = -50.0;
    force[i].z = 25.0;
  }

  ib_mesh_t next = ib_force_coupled_get_next(model, force, N, 0.2);
  ASSERT_EQ(next.n, N);
  for (int i = 0; i < N; ++i) {
    EXPECT_NEAR(next.position[i].x, current.position[i].x, 1e-12);
    EXPECT_NEAR(next.position[i].y, current.position[i].y, 1e-12);
    EXPECT_NEAR(next.position[i].z, current.position[i].z, 1e-12);
    EXPECT_NEAR(next.velocity[i].x, 0.0, 1e-12);
    EXPECT_NEAR(next.velocity[i].y, 0.0, 1e-12);
    EXPECT_NEAR(next.velocity[i].z, 0.0, 1e-12);
  }

  ib_mesh_free(&current);
  ib_mesh_free(&next);
  std::free(force);
  ib_pinned_rigid_body_sphere_destroy(sphere);
}

} // namespace
