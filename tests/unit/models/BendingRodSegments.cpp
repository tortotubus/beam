#include "elff/models/articulated/BendingRodSegments.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <stdexcept>

using namespace ELFF::Models;

TEST(BendingRodSegmentsTest, PinnedSingleSegmentHasPendulumAcceleration)
{
  BendingRodSegmentsConfig config;
  config.links = 1;
  config.length = 2.0;
  config.pin = Vec3(1.0, 2.0, 3.0);
  BendingRodSegments rod(config);
  const auto state = rod.initialState();
  const auto derivative = rod.rhs(0.0, state);
  EXPECT_NEAR(derivative[5], -3.0 * 9.81 / (2.0 * config.length), 1e-12);
  EXPECT_TRUE(rod.nodes(state).row(0).transpose().isApprox(config.pin));
}

TEST(BendingRodSegmentsTest, FreeRodFallsWithoutChangingShape)
{
  BendingRodSegmentsConfig config;
  config.links = 3;
  config.boundary = Boundary::Free;
  config.masses = {1.0, 2.0, 3.0};
  BendingRodSegments rod(config);
  const Vec3 center(1.0, 2.0, 3.0);
  const Vec3 velocity(0.4, 0.3, 0.2);
  const auto initial = rod.initialState({}, {}, center, velocity);
  const double duration = 0.1;
  const auto final = integrateRK45(rod, initial, duration, 0.01);
  const Vec3 expected = center + duration * velocity +
                        0.5 * duration * duration * config.gravity;
  EXPECT_TRUE(final.head<3>().isApprox(expected, 1e-10));
  EXPECT_NEAR(rod.energy(initial), rod.energy(final), 1e-9);
  const auto points = rod.nodes(final);
  Vec3 measured_center = Vec3::Zero();
  for (int i = 0; i < config.links; ++i) {
    EXPECT_NEAR((points.row(i + 1) - points.row(i)).norm(),
                config.length / config.links, 1e-12);
    measured_center += config.masses[i] *
                       0.5 * (points.row(i) + points.row(i + 1)).transpose();
  }
  EXPECT_TRUE((measured_center / 6.0).isApprox(expected, 1e-10));
}

TEST(BendingRodSegmentsTest, BendingRestoresAndDampingDissipates)
{
  BendingRodSegmentsConfig config;
  config.links = 2;
  config.gravity.setZero();
  config.ei = 1.0;
  config.bending_damping = {0.2};
  BendingRodSegments rod(config);
  const auto initial = rod.initialState({Vec3::UnitX(), Vec3(1.0, 0.2, 0.0)});
  const auto final = integrateRK45(rod, initial, 0.05, 0.001, 1e-9, 1e-12);
  EXPECT_LT(rod.energy(final), rod.energy(initial));
  const auto positions = rod.nodes(final);
  const Vec3 first = (positions.row(1) - positions.row(0)).normalized();
  const Vec3 second = (positions.row(2) - positions.row(1)).normalized();
  EXPECT_GT(first.dot(second), 1.0 / std::sqrt(1.04));
}

TEST(BendingRodSegmentsTest, RejectsInvalidConfigurationAndState)
{
  BendingRodSegmentsConfig config;
  config.links = 0;
  EXPECT_THROW(BendingRodSegments rod(config), std::invalid_argument);
  config.links = 2;
  config.masses = {1.0};
  EXPECT_THROW(BendingRodSegments rod(config), std::invalid_argument);
  config.masses.clear();
  BendingRodSegments rod(config);
  EXPECT_THROW(rod.initialState({Vec3::Zero()}), std::invalid_argument);
  EXPECT_THROW(rod.rhs(0.0, Eigen::VectorXd::Zero(1)), std::invalid_argument);
}

TEST(BendingRodSegmentsTest, FixedStepsAreDefaultAndIgnoreErrorTolerances)
{
  BendingRodSegmentsConfig config;
  config.links = 1;
  BendingRodSegments rod(config);
  const auto initial = rod.initialState();
  const auto result = integrateRK45(rod, initial, 0.25, 0.125);
  auto manual = integrateRK45(rod, initial, 0.125, 0.125);
  manual = integrateRK45(rod, manual, 0.125, 0.125);
  EXPECT_TRUE(result.isApprox(manual, 1e-14));
  EXPECT_TRUE(result.isApprox(
      integrateRK45(rod, initial, 0.25, 0.125, 1e-14, 1e-16), 1e-14));
  const auto shortened = integrateRK45(rod, initial, 0.1875, 0.125);
  manual = integrateRK45(rod, initial, 0.125, 0.125);
  manual = integrateRK45(rod, manual, 0.0625, 0.0625);
  EXPECT_TRUE(shortened.isApprox(manual, 1e-14));
}

TEST(BendingRodSegmentsTest, AdaptiveOptionMatchesFineStepReference)
{
  BendingRodSegmentsConfig config;
  config.links = 1;
  BendingRodSegments rod(config);
  const auto initial = rod.initialState();
  const auto reference = integrateRK45(rod, initial, 0.5, 0.001);
  const auto adaptive = integrateRK45(rod, initial, 0.5, 0.5, 1e-10, 1e-12, true);
  EXPECT_LT((adaptive - reference).norm(), 1e-8);
  const auto coarse = integrateRK45(rod, initial, 0.5, 0.25);
  const auto fine = integrateRK45(rod, initial, 0.5, 0.125);
  EXPECT_LT((fine - reference).norm(), (coarse - reference).norm());
}

TEST(BendingRodSegmentsTest, ForceCoupledInterfaceAdvancesUniformlyLoadedFreeRod)
{
  BendingRodSegmentsConfig config;
  config.links = 3;
  config.length = 2.0;
  config.total_mass = 4.0;
  config.gravity.setZero();
  config.boundary = Boundary::Free;
  BendingRodSegments rod(config);
  IBForceCoupled &model = rod;
  ASSERT_EQ(model.GetNumberOfPoints(), 4u);
  double measure = 0.0;
  for (auto weight : model.GetCurrent().GetMeasures()) measure += weight;
  EXPECT_NEAR(measure, config.length, 1e-12);
  EXPECT_NEAR(model.GetCurrent().GetMeasures().front(), config.length / 6.0, 1e-12);
  std::vector<IBMesh::IBVertex> forces(4, {0, 3, 0});
  const double dt = 0.01;
  model.GetNext(forces.data(), forces.size(), dt);
  EXPECT_NEAR(rod.state()[1], 0.5 * 1.5 * dt * dt, 1e-12);
  EXPECT_NEAR(rod.state()[3 + 3 * config.links + 1], 1.5 * dt, 1e-12);
  for (const auto &velocity : model.GetCurrent().GetVelocity()) {
    EXPECT_NEAR(velocity.y, 1.5 * dt, 1e-12);
  }
  const auto &next = model.GetNext(forces, dt);
  EXPECT_NEAR(rod.state()[1], 0.5 * 1.5 * 4.0 * dt * dt, 1e-12);
  EXPECT_EQ(next.GetPoints()[0].y, model.GetCurrent().GetPoints()[0].y);
}

TEST(BendingRodSegmentsTest, EndpointLoadsProduceCorrectTorque)
{
  BendingRodSegmentsConfig config;
  config.links = 1;
  config.length = 2.0;
  config.total_mass = 3.0;
  config.gravity.setZero();
  BendingRodSegments pinned(config);
  std::vector<Vec3> loads = {Vec3(0, 100, 0), Vec3(0, 4, 0)};
  const auto derivative = pinned.rhs(0, pinned.state(), loads);
  // Pin inertia is m L^2 / 3; the load at the pin contributes no torque.
  EXPECT_NEAR(derivative[5], 2.0, 1e-12);
  config.boundary = Boundary::Free;
  BendingRodSegments free(config);
  loads = {Vec3(0, -4, 0), Vec3(0, 4, 0)};
  const auto free_derivative = free.rhs(0, free.state(), loads);
  EXPECT_NEAR(free_derivative[11], 8.0, 1e-12);
  EXPECT_NEAR(free_derivative.segment<3>(6).norm(), 0.0, 1e-12);
}

TEST(BendingRodSegmentsTest, MeshVelocitiesMatchMovingCurvedGeometry)
{
  BendingRodSegmentsConfig config;
  config.links = 2;
  config.boundary = Boundary::Free;
  config.masses = {1, 3};
  BendingRodSegments rod(config);
  rod.setState(rod.initialState({Vec3::UnitX(), Vec3::UnitY()},
                               {Vec3(0, 0, 1), Vec3(0, 0, -2)},
                               Vec3(1, 2, 3), Vec3(0.5, -0.3, 0.1)));
  const auto derivative = rod.rhs(0, rod.state());
  const double eps = 1e-6;
  const Eigen::MatrixXd velocity =
      (rod.nodes(rod.state() + eps * derivative) -
       rod.nodes(rod.state() - eps * derivative)) / (2 * eps);
  for (int i = 0; i <= config.links; ++i) {
    const auto v = rod.GetCurrent().GetVelocity()[i];
    EXPECT_NEAR(v.x, velocity(i, 0), 1e-9);
    EXPECT_NEAR(v.y, velocity(i, 1), 1e-9);
    EXPECT_NEAR(v.z, velocity(i, 2), 1e-9);
  }
}

TEST(BendingRodSegmentsTest, CoupledSubstepsAndRestartPreserveContinuation)
{
  BendingRodSegmentsConfig config;
  config.links = 2;
  config.ei = 0.1;
  config.integration_step = 0.001;
  BendingRodSegments rod(config), manual(config), restored(config);
  std::vector<IBMesh::IBVertex> force(3, {0.1, 0.2, 0.3});
  rod.GetNext(force, 0.004);
  for (int i = 0; i < 4; ++i) manual.GetNext(force, 0.001);
  EXPECT_TRUE(rod.state().isApprox(manual.state(), 1e-12));
  IBModelState checkpoint;
  IBModel &model = rod;
  model.pack_state(checkpoint);
  restored.unpack_state(checkpoint);
  EXPECT_TRUE(rod.state().isApprox(restored.state(), 1e-14));
  EXPECT_EQ(restored.GetCurrent().GetForces()[0].z, force[0].z);
  rod.GetNext(force, 0.004);
  restored.GetNext(force, 0.004);
  EXPECT_TRUE(rod.state().isApprox(restored.state(), 1e-14));
  config.length = 2;
  BendingRodSegments different(config);
  EXPECT_THROW(different.unpack_state(checkpoint), std::invalid_argument);
  checkpoint.reals.pop_back();
  EXPECT_THROW(restored.unpack_state(checkpoint), std::invalid_argument);
}

TEST(BendingRodSegmentsTest, InvalidCouplingInputLeavesStateUnchanged)
{
  BendingRodSegmentsConfig config;
  config.links = 1;
  BendingRodSegments rod(config);
  const Eigen::VectorXd initial = rod.state();
  EXPECT_THROW(rod.ComputeNextPoints({}, 0.01), std::invalid_argument);
  EXPECT_THROW(rod.GetNext(std::vector<IBMesh::IBVertex>(2), -0.01),
               std::invalid_argument);
  EXPECT_THROW(rod.setState(Eigen::VectorXd::Zero(rod.stateSize())),
               std::invalid_argument);
  EXPECT_TRUE(rod.state().isApprox(initial));
}
