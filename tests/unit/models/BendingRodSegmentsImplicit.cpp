#include <elff/models/articulated/BendingRodSegmentsImplicit.hpp>
#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <stdexcept>

namespace ELFF {
namespace {
using namespace Models;

TEST(BendingRodSegmentsImplicitTest, StepRefinementIsSecondOrder)
{
  for (const double radius : {0.0, 0.5, 1.0})
  for (const auto boundary : {Boundary::Pinned, Boundary::Free}) {
    SCOPED_TRACE(boundary == Boundary::Pinned ? "pinned" : "free");
    BendingRodSegmentsImplicitConfig config;
    config.links = 3;
    config.spectral_radius = radius;
    SCOPED_TRACE(radius);
    config.nonlinear_tolerance = 1e-12;
    config.boundary = boundary;
    config.ei = 0.05;
    config.masses = {1.0, 1.5, 0.8};
    config.bending_damping = {0.01, 0.02};
    BendingRodSegmentsImplicit rod(config);
    const auto initial = rod.initialState(
        {Vec3(1, 0.1, 0.1), Vec3(1, -0.1, 0.2), Vec3(1, 0.2, -0.1)},
        {Vec3(0.2, -0.1, 0.3)}, Vec3(0.1, 0.2, -0.1), Vec3(0.1, 0.2, 0.3));
    std::vector<IBMesh::IBVertex> loads(4, {0.1, 0.2, -0.1});
    std::vector<Vec3> integrated_loads(4);
    for (int i = 0; i < 4; ++i)
      integrated_loads[i] = rod.GetCurrent().GetMeasures()[i] * Vec3(0.1, 0.2, -0.1);
    const auto reference = integrateRK45(rod, initial, 0.1, 1e-5,
                                         1e-6, 1e-9, false, integrated_loads);
    double errors[3];
    for (int level = 0; level < 3; ++level) {
      rod.setState(initial);
      const int steps = 50 * (1 << level);
      for (int i = 0; i < steps; ++i) rod.GetNext(loads, 0.1 / steps);
      errors[level] = (rod.state() - reference).norm();
    }
    for (int level = 1; level < 3; ++level) {
      const double order = std::log2(errors[level - 1] / errors[level]);
      EXPECT_GT(order, 1.8);
      EXPECT_LT(order, 2.2);
      GTEST_LOG_(INFO) << (boundary == Boundary::Pinned ? "pinned" : "free")
                       << " implicit observed order=" << order;
    }
    EXPECT_LT(errors[2], 1e-4);
  }
}

TEST(BendingRodSegmentsImplicitTest, StiffBendingRemainsBoundedAndPreservesLengths)
{
  BendingRodSegmentsImplicitConfig config;
  config.links = 8;
  config.ei = 100;
  config.gravity.setZero();
  BendingRodSegmentsImplicit rod(config);
  std::vector<Vec3> tangents(config.links);
  for (int i = 0; i < config.links; ++i)
    tangents[i] = Vec3(1, i % 2 == 0 ? 0.02 : -0.02, 0.01);
  rod.setState(rod.initialState(tangents));
  const double initial_energy = rod.energy(rod.state());
  std::vector<IBMesh::IBVertex> loads(9, {0, 0, 0});
  for (int i = 0; i < 100; ++i) {
    SCOPED_TRACE(i);
    ASSERT_NO_THROW(rod.GetNext(loads, 0.01));
    ASSERT_TRUE(rod.state().allFinite());
    // Generalized-alpha can overshoot during startup; it is not a
    // monotonically energy-decreasing method. Check boundedness and decay.
    EXPECT_LE(rod.energy(rod.state()), initial_energy * 2.0);
    const auto positions = rod.nodes(rod.state());
    EXPECT_NEAR(positions.row(0).norm(), 0.0, 1e-14);
    for (int j = 0; j < config.links; ++j) {
      EXPECT_NEAR((positions.row(j + 1) - positions.row(j)).norm(), 0.125, 1e-12);
      const Vec3 tangent = rod.state().segment<3>(3 * j);
      const Vec3 omega = rod.state().segment<3>(3 * config.links + 3 * j);
      EXPECT_NEAR(tangent.norm(), 1.0, 1e-12);
      EXPECT_NEAR(tangent.dot(omega), 0.0, 1e-12);
    }
  }
  EXPECT_LT(rod.energy(rod.state()), initial_energy * 0.01);
}

TEST(BendingRodSegmentsImplicitTest, ForceCouplingAndRestart)
{
  BendingRodSegmentsImplicitConfig config;
  config.links = 3;
  config.boundary = Boundary::Free;
  config.gravity.setZero();
  BendingRodSegmentsImplicit rod(config), restored(config);
  IBForceCoupled &coupled = rod;
  std::vector<IBMesh::IBVertex> force(4, {0, 2, 0});
  coupled.GetNext(force, 0.01);
  EXPECT_NEAR(rod.state()[1], 0.0001, 1e-12);
  EXPECT_NEAR(coupled.GetCurrent().GetVelocity()[0].y, 0.02, 1e-12);
  IBModelState checkpoint;
  coupled.pack_state(checkpoint);
  restored.unpack_state(checkpoint);
  EXPECT_TRUE(restored.state().isApprox(rod.state(), 1e-14));
  rod.GetNext(force, 0.01);
  restored.GetNext(force, 0.01);
  EXPECT_TRUE(restored.state().isApprox(rod.state(), 1e-14));
  IBModelState wrong_scheme = checkpoint;
  wrong_scheme.ints[4] = 2;
  EXPECT_THROW(restored.unpack_state(wrong_scheme), std::invalid_argument);
  BendingRodSegmentsConfig rk_config;
  static_cast<BendingRodSegmentsParameters &>(rk_config) = config;
  BendingRodSegments rk(rk_config);
  EXPECT_THROW(rk.unpack_state(checkpoint), std::invalid_argument);
  rk.pack_state(checkpoint);
  EXPECT_THROW(restored.unpack_state(checkpoint), std::invalid_argument);
}

TEST(BendingRodSegmentsImplicitTest, FullStepLargeSwingRemainsBounded)
{
  // Regression: implicit bending with explicit centrifugal terms blew up
  // near t=8 in this case, despite passing small-amplitude stiffness checks.
  BendingRodSegmentsImplicitConfig config;
  config.links = 60;
  config.length = 32.6;
  config.total_mass = 7.67 * config.length;
  config.ei = 700;
  BendingRodSegmentsImplicit rod(config);
  std::vector<IBMesh::IBVertex> force(61, {0, 0, 0});
  const double energy_scale = config.total_mass * config.gravity.norm() * config.length;
  for (int i = 0; i < 1000; ++i) {
    SCOPED_TRACE(i);
    ASSERT_NO_THROW(rod.GetNext(force, 0.01));
    ASSERT_TRUE(rod.state().allFinite());
    // This is a boundedness check, not an energy-conservation claim for the
    // generalized-alpha scheme. Second-order accuracy is checked separately above.
    ASSERT_LT(std::abs(rod.energy(rod.state())), 0.1 * energy_scale);
  }
}

TEST(BendingRodSegmentsImplicitTest, RotatingCheckpointPreservesHistory)
{
  BendingRodSegmentsImplicitConfig config;
  config.links = 3;
  config.ei = 0.2;
  BendingRodSegmentsImplicit rod(config), restored(config);
  std::vector<IBMesh::IBVertex> force(4, {0.1, -0.2, 0.3});
  rod.setState(rod.initialState({Vec3(1, 0.2, 0), Vec3(1, 0, 0.1), Vec3(1, -0.1, 0)},
                                {Vec3(0.1, 0.2, 0.3)}));
  IBModelState checkpoint;
  rod.pack_state(checkpoint); // Unstarted history is also valid.
  restored.unpack_state(checkpoint);
  for (int i = 0; i < 4; ++i) rod.GetNext(force, 0.01);
  rod.pack_state(checkpoint);
  restored.unpack_state(checkpoint);
  for (double dt : {0.01, 0.01, 0.005, 0.005}) {
    rod.GetNext(force, dt);
    restored.GetNext(force, dt);
    EXPECT_TRUE(restored.state().isApprox(rod.state(), 1e-13));
  }
  checkpoint.reals.back() = std::numeric_limits<double>::quiet_NaN();
  const auto before = restored.state();
  EXPECT_THROW(restored.unpack_state(checkpoint), std::invalid_argument);
  EXPECT_TRUE(restored.state().isApprox(before, 0));
}

TEST(BendingRodSegmentsImplicitTest, NonconvergenceDoesNotAdvanceOrChangeHistory)
{
  BendingRodSegmentsImplicitConfig config;
  config.links = 3;
  config.ei = 1;
  config.max_iterations = 1;
  config.nonlinear_tolerance = 1e-13;
  BendingRodSegmentsImplicit rod(config);
  rod.setState(rod.initialState({Vec3(1, 0.2, 0), Vec3(1, 0, 0.1), Vec3(1, -0.1, 0)}));
  IBModelState before, after;
  rod.pack_state(before);
  std::vector<IBMesh::IBVertex> force(4, {0, 0, 0});
  EXPECT_THROW(rod.GetNext(force, 0.1), std::runtime_error);
  rod.pack_state(after);
  EXPECT_EQ(before.ints, after.ints);
  EXPECT_EQ(before.reals, after.reals);
}

TEST(BendingRodSegmentsImplicitTest, RejectsInvalidConfigurationAndRestartsHistory)
{
  BendingRodSegmentsImplicitConfig config;
  config.links = 2;
  config.spectral_radius = -0.1;
  EXPECT_THROW(BendingRodSegmentsImplicit rod(config), std::invalid_argument);
  config.spectral_radius = 1.1;
  EXPECT_THROW(BendingRodSegmentsImplicit rod(config), std::invalid_argument);
  config.spectral_radius = 0.5;
  config.max_iterations = 0;
  EXPECT_THROW(BendingRodSegmentsImplicit rod(config), std::invalid_argument);
  config.max_iterations = 20;
  BendingRodSegmentsImplicit rod(config), fresh(config);
  std::vector<IBMesh::IBVertex> force(3, {0, 0, 0});
  rod.GetNext(force, 0.01);
  // Reset through the base interface must also discard algorithmic history.
  BendingRodSegments &base = rod;
  base.setState(fresh.state());
  rod.GetNext(force, 0.01);
  fresh.GetNext(force, 0.01);
  EXPECT_TRUE(rod.state().isApprox(fresh.state(), 1e-14));
  const auto before = rod.state();
  EXPECT_THROW(rod.GetNext(force, -0.01), std::invalid_argument);
  EXPECT_TRUE(rod.state().isApprox(before, 0));
}

} // namespace
} // namespace ELFF
