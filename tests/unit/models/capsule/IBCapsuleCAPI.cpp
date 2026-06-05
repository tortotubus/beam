#include <gtest/gtest.h>

#include <elff/c/models/capsule/IBCapsule.h>
#include <elff/c/models/ibm/IBMesh.h>
#include <elff/c/models/ibm/IBVelocityCoupled.h>

#include <vector>

namespace {

vertex_t vertex(double x, double y, double z)
{
  return { x, y, z };
}

void expect_valid_current_mesh(ib_capsule_t capsule, int expected_nodes)
{
  ib_velocity_coupled_t vc = ib_capsule_as_velocity_coupled(capsule);
  ib_mesh_t mesh = ib_velocity_coupled_get_current(vc);

  EXPECT_EQ(mesh.n, expected_nodes);
  ASSERT_NE(mesh.position, nullptr);
  ASSERT_NE(mesh.velocity, nullptr);
  ASSERT_NE(mesh.forces, nullptr);
  ASSERT_NE(mesh.measure, nullptr);
  EXPECT_GT(mesh.measure[0], 0.0);

  ib_mesh_free(&mesh);
}

} // namespace

TEST(IBCapsuleCAPITest, CreatesBuilderCapsulesAsVelocityCoupledModels)
{
  ib_capsule_t sphere =
    ib_capsule_sphere_new(1.0, vertex(0.0, 0.0, 0.0), 1, 1.0);
  ib_capsule_t ellipsoid =
    ib_capsule_ellipsoid_new(vertex(1.0, 0.75, 0.5),
                             vertex(0.25, 0.5, 0.75),
                             1,
                             1.0);
  ib_capsule_t biconcave =
    ib_capsule_biconcave_new(1.0, vertex(0.0, 0.0, 0.0), 1, 1.0);

  expect_valid_current_mesh(sphere, 42);
  expect_valid_current_mesh(ellipsoid, 42);
  expect_valid_current_mesh(biconcave, 42);

  ib_capsule_destroy(sphere);
  ib_capsule_destroy(ellipsoid);
  ib_capsule_destroy(biconcave);
}

TEST(IBCapsuleCAPITest, AdvancesThroughVelocityCoupledInterface)
{
  ib_capsule_t capsule =
    ib_capsule_sphere_new(1.0, vertex(0.0, 0.0, 0.0), 0, 1.0);
  ib_velocity_coupled_t vc = ib_capsule_as_velocity_coupled(capsule);

  ib_mesh_t current = ib_velocity_coupled_get_current(vc);
  ASSERT_EQ(current.n, 12);
  const vertex_t x0 = current.position[0];

  std::vector<vertex_t> velocity(static_cast<size_t>(current.n),
                                 vertex(1.0, 2.0, 3.0));

  ib_mesh_t next =
    ib_velocity_coupled_get_next(vc, velocity.data(), current.n, 0.25);

  ASSERT_EQ(next.n, current.n);
  EXPECT_DOUBLE_EQ(next.position[0].x, x0.x + 0.25);
  EXPECT_DOUBLE_EQ(next.position[0].y, x0.y + 0.5);
  EXPECT_DOUBLE_EQ(next.position[0].z, x0.z + 0.75);
  EXPECT_DOUBLE_EQ(next.velocity[0].x, 1.0);
  EXPECT_DOUBLE_EQ(next.velocity[0].y, 2.0);
  EXPECT_DOUBLE_EQ(next.velocity[0].z, 3.0);

  ib_mesh_free(&current);
  ib_mesh_free(&next);
  ib_capsule_destroy(capsule);
}

TEST(IBCapsuleCAPITest, SetsMembraneAndBendingLaws)
{
  ib_capsule_t capsule =
    ib_capsule_sphere_new(1.0, vertex(0.0, 0.0, 0.0), 1, 1.0);

  ib_capsule_set_skalak_law(capsule, 2.0, 10.0);
  ib_capsule_set_neo_hookean_law(capsule, 1.5);
  ib_capsule_set_linear_bending_law(capsule, 0.01);
  ib_capsule_set_helfrich_bending_law(capsule, 0.02);
  ib_capsule_clear_bending_law(capsule);

  expect_valid_current_mesh(capsule, 42);
  ib_capsule_destroy(capsule);
}

TEST(IBCapsuleCAPITest, ExposesTriangleTopology)
{
  ib_capsule_t capsule =
    ib_capsule_sphere_new(1.0, vertex(0.0, 0.0, 0.0), 0, 1.0);

  const int triangle_count = ib_capsule_get_triangle_count(capsule);
  EXPECT_EQ(triangle_count, 20);

  int nodes[3] = { -1, -1, -1 };
  ib_capsule_get_triangle_node_ids(capsule, 0, nodes);

  EXPECT_GE(nodes[0], 0);
  EXPECT_GE(nodes[1], 0);
  EXPECT_GE(nodes[2], 0);
  EXPECT_LT(nodes[0], 12);
  EXPECT_LT(nodes[1], 12);
  EXPECT_LT(nodes[2], 12);
  EXPECT_NE(nodes[0], nodes[1]);
  EXPECT_NE(nodes[1], nodes[2]);
  EXPECT_NE(nodes[2], nodes[0]);

  ib_capsule_destroy(capsule);
}
