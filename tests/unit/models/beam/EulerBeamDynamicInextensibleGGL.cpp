#if __has_include(<elff/models/beam/EulerBeamDynamicInextensibleGGL.hpp>)

#include <gtest/gtest.h>

#include <elff/models/beam/EulerBeamDynamicInextensibleGGL.hpp>

#include <cmath>

namespace ELFF {
namespace Models {
namespace {

void
expect_mesh_state_finite(EulerBeamMesh& mesh)
{
  for (size_t i = 0; i < mesh.get_nodes(); ++i) {
    const auto& centerline = mesh.get_centerline(i);
    const auto& slope = mesh.get_slope(i);
    const auto& velocity = mesh.get_centerline_velocity(i);

    for (real_t value : centerline) {
      EXPECT_TRUE(std::isfinite(value));
    }
    for (real_t value : slope) {
      EXPECT_TRUE(std::isfinite(value));
    }
    for (real_t value : velocity) {
      EXPECT_TRUE(std::isfinite(value));
    }
  }
}

TEST(EulerBeamDynamicInextensibleGGLUnitTest,
     KeepsSimpleBoundaryFixedAndStateFinite)
{
  const real_t length = 1.;
  const real_t EI = 0.05;
  const real_t mu = 1.;
  const real_t r_penalty = 1e4;
  const size_t nodes = 12;
  const real_t dt = 1e-2;
  const std::array<real_t, 3> load = { 0., -0.1, 0. };

  EulerBeam::EulerBeamBCs boundary_conditions = {
    .end = { EulerBeam::left, EulerBeam::right },
    .type = { EulerBeam::simple_bc, EulerBeam::free_bc },
    .vals = { {
                .position = { 0., 0., 0. },
              },
              {
              } }
  };

  EulerBeamDynamicInextensibleGGL beam(
    length, EI, mu, nodes, boundary_conditions, r_penalty);
  beam.apply_initial_condition();

  for (size_t step = 0; step < 5; ++step) {
    ASSERT_NO_THROW(beam.solve(dt, load));
  }

  EulerBeamMesh& mesh = beam.get_mesh();
  expect_mesh_state_finite(mesh);

  EXPECT_NEAR(mesh.get_centerline(0)[0], 0., 1e-12);
  EXPECT_NEAR(mesh.get_centerline(0)[1], 0., 1e-12);
  EXPECT_NEAR(mesh.get_centerline(0)[2], 0., 1e-12);

  EXPECT_NEAR(mesh.get_centerline_velocity(0)[0], 0., 1e-12);
  EXPECT_NEAR(mesh.get_centerline_velocity(0)[1], 0., 1e-12);
  EXPECT_NEAR(mesh.get_centerline_velocity(0)[2], 0., 1e-12);
}

TEST(EulerBeamDynamicInextensibleGGLUnitTest,
     PreservesClampedPositionAndSlopeAfterStep)
{
  const real_t length = 1.;
  const real_t EI = 0.05;
  const real_t mu = 1.;
  const real_t r_penalty = 1e4;
  const size_t nodes = 10;
  const real_t dt = 1e-3;
  const std::array<real_t, 3> load = { 0., 0., 0. };

  EulerBeam::EulerBeamBCs boundary_conditions = {
    .end = { EulerBeam::left, EulerBeam::right },
    .type = { EulerBeam::clamped_bc, EulerBeam::free_bc },
    .vals = { {
                .position = { 0., 0., 0. },
                .slope = { 1., 0., 0. },
              },
              {
              } }
  };

  EulerBeamDynamicInextensibleGGL beam(
    length, EI, mu, nodes, boundary_conditions, r_penalty);
  beam.apply_initial_condition();

  ASSERT_NO_THROW(beam.solve(dt, load));

  EulerBeamMesh& mesh = beam.get_mesh();
  expect_mesh_state_finite(mesh);

  EXPECT_NEAR(mesh.get_centerline(0)[0], 0., 1e-12);
  EXPECT_NEAR(mesh.get_centerline(0)[1], 0., 1e-12);
  EXPECT_NEAR(mesh.get_centerline(0)[2], 0., 1e-12);

  EXPECT_NEAR(mesh.get_slope(0)[0], 1., 1e-12);
  EXPECT_NEAR(mesh.get_slope(0)[1], 0., 1e-12);
  EXPECT_NEAR(mesh.get_slope(0)[2], 0., 1e-12);

  EXPECT_NEAR(mesh.get_centerline_velocity(0)[0], 0., 1e-12);
  EXPECT_NEAR(mesh.get_centerline_velocity(0)[1], 0., 1e-12);
  EXPECT_NEAR(mesh.get_centerline_velocity(0)[2], 0., 1e-12);
}

} // namespace
} // namespace Models
} // namespace ELFF

#else

#include <gtest/gtest.h>

TEST(EulerBeamDynamicInextensibleGGLUnitTest, HeaderMissing)
{
  GTEST_SKIP()
    << "EulerBeamDynamicInextensibleGGL.hpp is not available in this workspace.";
}

#endif
