#include <elff/models/articulated/BendingRodSegmentsIMEX.hpp>
#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <stdexcept>

namespace ELFF {
namespace {
using namespace Models;

TEST(BendingRodSegmentsIMEXTest, SmallStepMatchesArticulatedDynamics)
{
  for (const auto boundary : {Boundary::Pinned, Boundary::Free}) {
    SCOPED_TRACE(boundary == Boundary::Pinned ? "pinned" : "free");
    BendingRodSegmentsIMEXConfig config;
    config.links = 3;
    config.boundary = boundary;
    config.masses = {1.0, 2.0, 3.0};
    config.ei = 0.2;
    config.bending_damping = {0.03, 0.04};
    BendingRodSegmentsIMEX rod(config);
    const auto state = rod.initialState(
        {Vec3(1, 0.1, 0), Vec3(1, 0.2, 0.1), Vec3(1, -0.1, 0.2)},
        {Vec3(0.2, 0.3, -0.1)}, Vec3(1, 2, 3), Vec3(0.1, 0.2, 0.3));
    std::vector<Vec3> forces = {Vec3(1, 2, 3), Vec3(-1, 0, 1),
                                Vec3(0, -1, 2), Vec3(0.3, 0.5, -0.2)};
    const Eigen::VectorXd exact = rod.rhs(0, state, forces);
    const double dt = 1e-7;
    const Eigen::VectorXd measured = (rod.imexStep(state, dt, forces) - state) / dt;
    EXPECT_LT((measured - exact).norm() / (1.0 + exact.norm()), 1e-4);
  }
}

TEST(BendingRodSegmentsIMEXTest, StepRefinementIsSecondOrder)
{
  for (const auto boundary : {Boundary::Pinned, Boundary::Free}) {
    SCOPED_TRACE(boundary == Boundary::Pinned ? "pinned" : "free");
    BendingRodSegmentsIMEXConfig config;
    config.links = 3;
    config.boundary = boundary;
    config.ei = 0.05;
    config.masses = {1.0, 1.5, 0.8};
    config.bending_damping = {0.01, 0.02};
    BendingRodSegmentsIMEX rod(config);
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
                       << " IMEX observed order=" << order;
    }
    EXPECT_LT(errors[2], 1e-4);
  }
}

TEST(BendingRodSegmentsIMEXTest, StiffBendingRemainsBoundedAndPreservesLengths)
{
  BendingRodSegmentsIMEXConfig config;
  config.links = 8;
  config.ei = 100;
  config.gravity.setZero();
  BendingRodSegmentsIMEX rod(config);
  std::vector<Vec3> tangents(config.links);
  for (int i = 0; i < config.links; ++i)
    tangents[i] = Vec3(1, i % 2 == 0 ? 0.02 : -0.02, 0.01);
  rod.setState(rod.initialState(tangents));
  const double initial_energy = rod.energy(rod.state());
  std::vector<IBMesh::IBVertex> loads(9, {0, 0, 0});
  for (int i = 0; i < 100; ++i) {
    rod.GetNext(loads, 0.01);
    ASSERT_TRUE(rod.state().allFinite());
    EXPECT_LE(rod.energy(rod.state()), initial_energy * (1.0 + 1e-8));
    const auto positions = rod.nodes(rod.state());
    EXPECT_NEAR(positions.row(0).norm(), 0.0, 1e-14);
    for (int j = 0; j < config.links; ++j)
      EXPECT_NEAR((positions.row(j + 1) - positions.row(j)).norm(), 0.125, 1e-12);
  }
}

TEST(BendingRodSegmentsIMEXTest, ForceCouplingAndRestart)
{
  BendingRodSegmentsIMEXConfig config;
  config.links = 3;
  config.boundary = Boundary::Free;
  config.gravity.setZero();
  BendingRodSegmentsIMEX rod(config), restored(config);
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
  IBModelState first_order = checkpoint;
  first_order.ints.back() = 1;
  EXPECT_THROW(restored.unpack_state(first_order), std::invalid_argument);
  BendingRodSegmentsConfig rk_config;
  static_cast<BendingRodSegmentsParameters &>(rk_config) = config;
  BendingRodSegments rk(rk_config);
  EXPECT_THROW(rk.unpack_state(checkpoint), std::invalid_argument);
  rk.pack_state(checkpoint);
  EXPECT_THROW(restored.unpack_state(checkpoint), std::invalid_argument);
}

TEST(BendingRodSegmentsIMEXTest, FullStepLargeSwingRemainsBounded)
{
  // Regression: implicit bending with explicit centrifugal terms blew up
  // near t=8 in this case, despite passing small-amplitude stiffness checks.
  BendingRodSegmentsIMEXConfig config;
  config.links = 60;
  config.length = 32.6;
  config.total_mass = 7.67 * config.length;
  config.ei = 700;
  BendingRodSegmentsIMEX rod(config);
  std::vector<IBMesh::IBVertex> force(61, {0, 0, 0});
  const double energy_scale = config.total_mass * config.gravity.norm() * config.length;
  for (int i = 0; i < 1000; ++i) {
    SCOPED_TRACE(i);
    rod.GetNext(force, 0.01);
    ASSERT_TRUE(rod.state().allFinite());
    // This is a boundedness check, not an energy-conservation claim for the
    // extrapolated scheme. Second-order accuracy is checked separately above.
    ASSERT_LT(std::abs(rod.energy(rod.state())), 0.1 * energy_scale);
  }
}

TEST(BendingRodSegmentsIMEXTest, RejectsInvalidStepsAndAntiparallelJoints)
{
  BendingRodSegmentsIMEXConfig config;
  config.links = 2;
  BendingRodSegmentsIMEX rod(config);
  EXPECT_THROW(rod.imexStep(rod.state(), 0), std::invalid_argument);
  EXPECT_THROW(rod.imexStep(rod.state(), -0.01), std::invalid_argument);
  EXPECT_THROW(rod.imexStep(rod.state(), std::numeric_limits<double>::infinity()),
               std::invalid_argument);
  EXPECT_THROW(rod.imexStep(rod.state(), 0.01, {Vec3::Zero()}), std::invalid_argument);
  const auto folded = rod.initialState({Vec3::UnitX(), -Vec3::UnitX()});
  EXPECT_THROW(rod.imexStep(folded, 0.01), std::runtime_error);
}

} // namespace
} // namespace ELFF
