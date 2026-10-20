#include <elff/models/articulated/IBBendingRodImplicit.hpp>
#include <gtest/gtest.h>
#include <elff/c/models/ibm/IBForceCoupled.h>
#include <type_traits>

#include <cmath>
#include <limits>
#include <elff/general/error.hpp>

namespace ELFF {
namespace {
using namespace Models;

// Match the configured ELFF error policy; exception builds also verify rollback
// in-process, while abort builds exercise the same failures in a child process.
#ifdef ELFF_USE_EXCEPTIONS
#define EXPECT_ROD_ERROR(statement) EXPECT_THROW(statement, ELFF::ErrorException)
#else
#define EXPECT_ROD_ERROR(statement) EXPECT_DEATH({ statement; }, "BEAM abort:")
#endif

// Numerical diagnostics are exposed only by the test harness.
template<class Implementation>
class RodHarness : public Implementation {
public:
  using Implementation::Implementation;
  using Implementation::advance;
  using Implementation::state;
  using Implementation::setState;
  using Implementation::initialState;
  using Implementation::rhs;
  using Implementation::nodes;
  using Implementation::velocities;
  using Implementation::energy;
};
using TestRod = RodHarness<BendingRodImplicit>;
using TestIBRod = RodHarness<IBBendingRodImplicit>;


static_assert(std::is_abstract_v<BendingRod>);
static_assert(std::is_base_of_v<BendingRod, BendingRodImplicit>);

TEST(BendingRodImplicitTest, CommonMeshAndPolymorphicUniformSolve)
{
  BendingRodMesh initial(4, 2);
  EXPECT_EQ(initial.get_nodes(), 4u);
  EXPECT_DOUBLE_EQ(initial.get_curvilinear_axis(3), 2);
  EXPECT_DOUBLE_EQ(initial.get_ds(), 2.0 / 3);
  EXPECT_EQ(initial.get_centerline(3), (BendingRod::Vector{2, 0, 0}));
  EXPECT_ROD_ERROR(BendingRodMesh(1, 2));
  EXPECT_ROD_ERROR(BendingRodMesh(4, 0));
  BendingRodImplicitConfig config;
  config.links = 3;
  config.length = 2;
  config.total_mass = 4;
  config.gravity.setZero();
  config.boundary = Boundary::Free;
  BendingRodImplicit uniform(config), nonuniform(config);
  BendingRod &base = uniform;
  base.apply_initial_condition(initial);
  nonuniform.apply_initial_condition(initial);
  base.solve(0.01, BendingRod::Vector{0, 3, 0});
  nonuniform.solve(0.01, std::vector<BendingRod::Vector>(4, {0, 3, 0}));
  const BendingRod &view = uniform;
  for (size_t i = 0; i < 4; ++i) {
    EXPECT_EQ(view.get_mesh().get_centerline(i), nonuniform.get_mesh().get_centerline(i));
    EXPECT_NEAR(view.get_mesh().get_centerline(i)[1], 0.000075, 1e-12);
    EXPECT_NEAR(view.get_mesh().get_centerline_velocity(i)[1], 0.015, 1e-12);
  }
  // Reinitialize the moving mesh so the no-load step starts with zero
  // acceleration, rather than the previous load's generalized-alpha history.
  const BendingRodMesh moving = base.get_mesh();
  base.apply_initial_condition(moving);
  base.solve(0.01);
  EXPECT_NEAR(base.get_mesh().get_centerline(0)[1], 0.000225, 1e-12);
}

TEST(BendingRodImplicitTest, CurvedMovingMeshInitializesAndRestartsThroughBase)
{
  BendingRodImplicitConfig config;
  config.links = 2;
  config.length = 2;
  config.boundary = Boundary::Free;
  config.masses = {1, 3};
  IBBendingRodImplicit rod(config), restored(config);
  BendingRodMesh initial(3, 2);
  initial.get_centerline() = {{2, 3, 4}, {3, 3, 4}, {3, 4, 4}};
  initial.get_centerline_velocity() = {{0.2, 0.3, 0.4}, {0.2, 0.8, 0.4}, {0.45, 0.8, 0.4}};
  BendingRod &base = rod;
  base.apply_initial_condition(initial);
  for (size_t i = 0; i < 3; ++i) {
    for (size_t d = 0; d < 3; ++d) {
      EXPECT_NEAR(base.get_mesh().get_centerline(i)[d], initial.get_centerline(i)[d], 1e-12);
      EXPECT_NEAR(base.get_mesh().get_centerline_velocity(i)[d], initial.get_centerline_velocity(i)[d], 1e-12);
    }
    EXPECT_DOUBLE_EQ(rod.GetCurrent().GetPoints()[i].y, base.get_mesh().get_centerline(i)[1]);
  }
  base.solve(0.001);
  ModelState saved;
  rod.pack_state(saved);
  restored.unpack_state(saved);
  EXPECT_EQ(rod.get_mesh().get_centerline(), restored.get_mesh().get_centerline());
  EXPECT_EQ(rod.get_mesh().get_centerline_velocity(), restored.get_mesh().get_centerline_velocity());
  EXPECT_DOUBLE_EQ(restored.GetCurrent().GetPoints()[2].x, restored.get_mesh().get_centerline(2)[0]);
  base.apply_initial_condition(initial); // A mesh reset also discards algorithmic history.
  restored.apply_initial_condition(initial);
  base.solve(0.001);
  restored.solve(0.001);
  EXPECT_EQ(rod.get_mesh().get_centerline(), restored.get_mesh().get_centerline());
}

TEST(BendingRodImplicitTest, InvalidInitialMeshesLeaveStateAndHistoryUnchanged)
{
  BendingRodImplicitConfig config;
  config.links = 2;
  config.length = 2;
  IBBendingRodImplicit rod(config);
  rod.solve(0.001);
  ModelState before, after;
  rod.pack_state(before);
  for (int invalid = 0; invalid < 5; ++invalid) {
    BendingRodMesh initial(3, 2);
    if (invalid == 0) initial.get_centerline(0)[1] = 1; // Moves the pin.
    if (invalid == 1) initial.get_centerline(2)[0] += 0.1; // Stretches a segment.
    if (invalid == 2) initial.get_centerline_velocity(2)[0] = 1; // Axial velocity.
    if (invalid == 3) initial.get_centerline_velocity(0)[1] = 1; // Moving pin.
    if (invalid == 4) initial.get_centerline().pop_back();
    EXPECT_ROD_ERROR(rod.apply_initial_condition(initial));
    rod.pack_state(after);
    EXPECT_EQ(before.ints, after.ints);
    EXPECT_EQ(before.reals, after.reals);
  }
  EXPECT_ROD_ERROR(rod.apply_initial_condition(BendingRodMesh(4, 2)));
}

TEST(BendingRodImplicitTest, PinnedSingleSegmentHasPendulumAcceleration)
{
  BendingRodImplicitConfig config;
  config.links = 1;
  config.length = 2.0;
  config.pin = Vec3(1.0, 2.0, 3.0);
  TestRod rod(config);
  const auto state = rod.initialState();
  const auto derivative = rod.rhs(0.0, state);
  EXPECT_NEAR(derivative[5], -3.0 * 9.81 / (2.0 * config.length), 1e-12);
  EXPECT_TRUE(rod.nodes(state).row(0).transpose().isApprox(config.pin));
}

TEST(BendingRodImplicitTest, FreeRodFallsWithoutChangingShape)
{
  BendingRodImplicitConfig config;
  config.links = 3;
  config.boundary = Boundary::Free;
  config.masses = {1.0, 2.0, 3.0};
  TestRod rod(config);
  const Vec3 center(1.0, 2.0, 3.0);
  const Vec3 velocity(0.4, 0.3, 0.2);
  const auto initial = rod.initialState({}, {}, center, velocity);
  const double duration = 0.1;
  rod.setState(initial);
  for (int i = 0; i < 10; ++i) rod.advance(duration / 10);
  const auto final = rod.state();
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

TEST(BendingRodImplicitTest, BendingRestoresAndDampingDissipates)
{
  BendingRodImplicitConfig config;
  config.links = 2;
  config.gravity.setZero();
  config.ei = 1.0;
  config.bending_damping = {0.2};
  TestRod rod(config);
  const auto initial = rod.initialState({Vec3::UnitX(), Vec3(1.0, 0.2, 0.0)});
  rod.setState(initial);
  for (int i = 0; i < 50; ++i) rod.advance(0.001);
  const auto final = rod.state();
  EXPECT_LT(rod.energy(final), rod.energy(initial));
  const auto positions = rod.nodes(final);
  const Vec3 first = (positions.row(1) - positions.row(0)).normalized();
  const Vec3 second = (positions.row(2) - positions.row(1)).normalized();
  EXPECT_GT(first.dot(second), 1.0 / std::sqrt(1.04));
}

TEST(BendingRodImplicitTest, RejectsInvalidConfigurationAndState)
{
  BendingRodImplicitConfig config;
  config.links = 0;
  EXPECT_ROD_ERROR(TestRod rod(config));
  config.links = 2;
  config.masses = {1.0};
  EXPECT_ROD_ERROR(TestRod rod(config));
  config.masses.clear();
  TestRod rod(config);
  EXPECT_ROD_ERROR(rod.initialState({Vec3::Zero()}));
  EXPECT_ROD_ERROR(rod.rhs(0.0, Eigen::VectorXd::Zero(1)));
}

TEST(BendingRodImplicitTest, EndpointLoadsProduceCorrectTorque)
{
  BendingRodImplicitConfig config;
  config.links = 1;
  config.length = 2.0;
  config.total_mass = 3.0;
  config.gravity.setZero();
  TestRod pinned(config);
  std::vector<Vec3> loads = {Vec3(0, 100, 0), Vec3(0, 4, 0)};
  const auto derivative = pinned.rhs(0, pinned.state(), loads);
  // Pin inertia is m L^2 / 3; the load at the pin contributes no torque.
  EXPECT_NEAR(derivative[5], 2.0, 1e-12);
  config.boundary = Boundary::Free;
  TestRod free(config);
  loads = {Vec3(0, -4, 0), Vec3(0, 4, 0)};
  const auto free_derivative = free.rhs(0, free.state(), loads);
  EXPECT_NEAR(free_derivative[11], 8.0, 1e-12);
  EXPECT_NEAR(free_derivative.segment<3>(6).norm(), 0.0, 1e-12);
}

static_assert(!std::is_base_of_v<IBForceCoupled, BendingRodImplicit>);

TEST(BendingRodImplicitTest, IndependentSolverMatchesDensityAdapterAndVelocities)
{
  for (const auto boundary : {Boundary::Pinned, Boundary::Free}) {
    BendingRodImplicitConfig config;
    config.links = 3;
    config.length = 2;
    config.ei = 0.2;
    config.boundary = boundary;
    TestRod structural(config);
    TestIBRod coupled(config);
    const auto initial = structural.initialState(
      {Vec3(1, 0.2, 0), Vec3(1, 0, 0.1), Vec3(1, -0.1, 0)},
      {Vec3(0.1, 0.2, 0.3)}, Vec3(2, 3, 4), Vec3(0.3, 0.2, 0.1));
    structural.setState(initial);
    coupled.setState(initial);
    std::vector<IBMesh::IBVertex> density(4, {0, 0, 0});
    density.back() = {0.2, -0.1, 0.3}; // Distal load exercises torque transmission.
    std::vector<Vec3> loads(4, Vec3::Zero());
    loads.back() = config.length / (2 * config.links) * Vec3(0.2, -0.1, 0.3);
    for (double dt : {0.01, 0.01, 0.005}) {
      structural.advance(dt, loads);
      coupled.GetNext(density, dt);
      EXPECT_TRUE(structural.state().isApprox(coupled.state(), 1e-14));
      const auto positions = structural.nodes(structural.state());
      const auto velocities = structural.velocities(structural.state());
      const auto derivative = structural.rhs(0, structural.state(), loads);
      const Eigen::MatrixXd difference = (structural.nodes(structural.state() + 1e-7 * derivative) - positions) / 1e-7;
      EXPECT_LT((velocities - difference).norm(), 1e-6);
      EXPECT_NEAR(coupled.GetCurrent().GetVelocity().back().z, velocities(3, 2), 1e-13);
    }
    ModelState before, after;
    structural.pack_state(before);
    structural.advance(0, loads);
    structural.pack_state(after);
    EXPECT_EQ(before.reals, after.reals);
    EXPECT_ROD_ERROR(structural.advance(-0.1, loads));
    EXPECT_ROD_ERROR(structural.advance(0.1, {Vec3::Zero()}));
    structural.pack_state(after);
    EXPECT_EQ(before.reals, after.reals);
  }
}

TEST(BendingRodImplicitTest, CheckedCFailurePreservesStateAndMeshes)
{
  BendingRodImplicitConfig config;
  config.links = 3;
  config.ei = 1;
  config.max_iterations = 1;
  config.nonlinear_tolerance = 1e-13;
  TestIBRod rod(config);
  rod.setState(rod.initialState({Vec3(1, 0.2, 0), Vec3(1, 0, 0.1), Vec3(1, -0.1, 0)}));
  auto handle = static_cast<IBForceCoupled*>(&rod);
  std::vector<vertex_t> force(4, {0.1, 0.2, 0.3});
  char error[512];
  ib_mesh_t output{};
  ASSERT_EQ(ib_force_coupled_get_next_checked(handle, force.data(), 4, 0,
    &output, error, sizeof(error)), 0) << error;
  ib_mesh_free(&output);
  ModelState before, after;
  rod.pack_state(before);
  const auto mesh = rod.GetCurrent();
#ifdef ELFF_USE_EXCEPTIONS
  EXPECT_NE(ib_force_coupled_get_next_checked(handle, force.data(), 4, 0.1,
    &output, error, sizeof(error)), 0);
  EXPECT_NE(error[0], '\0');
  EXPECT_EQ(output.n, 0);
  EXPECT_EQ(output.position, nullptr);
#else
  EXPECT_DEATH(ib_force_coupled_get_next_checked(handle, force.data(), 4, 0.1,
    &output, error, sizeof(error)), "generalized-alpha Newton solve failed");
#endif
  rod.pack_state(after);
  EXPECT_EQ(before.ints, after.ints);
  EXPECT_EQ(before.reals, after.reals);
  EXPECT_DOUBLE_EQ(rod.GetCurrent().GetPoints().back().y, mesh.GetPoints().back().y);
  EXPECT_DOUBLE_EQ(rod.GetCurrent().GetVelocity().back().z, mesh.GetVelocity().back().z);
  EXPECT_NE(ib_force_coupled_get_next_checked(handle, force.data(), 3, 0.01,
    &output, error, sizeof(error)), 0);
  EXPECT_NE(ib_force_coupled_get_next_checked(handle, force.data(), 4, -0.01,
    &output, error, sizeof(error)), 0);
  force[0].x = std::numeric_limits<double>::quiet_NaN();
  EXPECT_NE(ib_force_coupled_get_next_checked(handle, force.data(), 4, 0.01,
    &output, error, sizeof(error)), 0);
  const auto legacy = ib_force_coupled_get_next(handle, force.data(), 4, 0.01);
  EXPECT_EQ(legacy.n, 0);
  rod.pack_state(after);
  EXPECT_EQ(before.reals, after.reals);
}

TEST(BendingRodImplicitTest, StepRefinementIsSecondOrder)
{
  for (const double radius : {0.0, 0.5, 1.0})
  for (const auto boundary : {Boundary::Pinned, Boundary::Free}) {
    SCOPED_TRACE(boundary == Boundary::Pinned ? "pinned" : "free");
    BendingRodImplicitConfig config;
    config.links = 3;
    config.spectral_radius = radius;
    SCOPED_TRACE(radius);
    config.nonlinear_tolerance = 1e-12;
    config.boundary = boundary;
    config.ei = 0.05;
    config.masses = {1.0, 1.5, 0.8};
    config.bending_damping = {0.01, 0.02};
    TestIBRod rod(config);
    const auto initial = rod.initialState(
        {Vec3(1, 0.1, 0.1), Vec3(1, -0.1, 0.2), Vec3(1, 0.2, -0.1)},
        {Vec3(0.2, -0.1, 0.3)}, Vec3(0.1, 0.2, -0.1), Vec3(0.1, 0.2, 0.3));
    std::vector<IBMesh::IBVertex> loads(4, {0.1, 0.2, -0.1});
    std::vector<Vec3> integrated_loads(4);
    for (int i = 0; i < 4; ++i)
      integrated_loads[i] = rod.GetCurrent().GetMeasures()[i] * Vec3(0.1, 0.2, -0.1);
    // A much finer implicit solution supplies the convergence reference.
    TestRod reference_rod(config);
    reference_rod.setState(initial);
    constexpr int reference_steps = 4000;
    for (int i = 0; i < reference_steps; ++i)
      reference_rod.advance(0.1 / reference_steps, integrated_loads);
    const Eigen::VectorXd reference = reference_rod.state();
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

TEST(BendingRodImplicitTest, StiffBendingRemainsBoundedAndPreservesLengths)
{
  BendingRodImplicitConfig config;
  config.links = 8;
  config.ei = 100;
  config.gravity.setZero();
  TestIBRod rod(config);
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

TEST(BendingRodImplicitTest, ForceCouplingAndRestart)
{
  BendingRodImplicitConfig config;
  config.links = 3;
  config.boundary = Boundary::Free;
  config.gravity.setZero();
  TestIBRod rod(config), restored(config);
  IBForceCoupled &coupled = rod;
  std::vector<IBMesh::IBVertex> force(4, {0, 2, 0});
  coupled.GetNext(force, 0.01);
  EXPECT_NEAR(rod.state()[1], 0.0001, 1e-12);
  EXPECT_NEAR(coupled.GetCurrent().GetVelocity()[0].y, 0.02, 1e-12);
  ModelState checkpoint;
  coupled.pack_state(checkpoint);
  restored.unpack_state(checkpoint);
  EXPECT_TRUE(restored.state().isApprox(rod.state(), 1e-14));
  rod.GetNext(force, 0.01);
  restored.GetNext(force, 0.01);
  EXPECT_TRUE(restored.state().isApprox(rod.state(), 1e-14));
  ModelState wrong_scheme = checkpoint;
  wrong_scheme.ints[4] = 2;
  EXPECT_ROD_ERROR(restored.unpack_state(wrong_scheme));
  // Structural and coupled checkpoints have distinct payloads.
  TestRod structural(config);
  EXPECT_ROD_ERROR(structural.unpack_state(checkpoint));
  structural.pack_state(checkpoint);
  EXPECT_ROD_ERROR(restored.unpack_state(checkpoint));

}

TEST(BendingRodImplicitTest, FullStepLargeSwingRemainsBounded)
{
  // Regression: implicit bending with explicit centrifugal terms blew up
  // near t=8 in this case, despite passing small-amplitude stiffness checks.
  BendingRodImplicitConfig config;
  config.links = 60;
  config.length = 32.6;
  config.total_mass = 7.67 * config.length;
  config.ei = 700;
  TestIBRod rod(config);
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

TEST(BendingRodImplicitTest, RotatingCheckpointPreservesHistory)
{
  BendingRodImplicitConfig config;
  config.links = 3;
  config.ei = 0.2;
  TestIBRod rod(config), restored(config);
  std::vector<IBMesh::IBVertex> force(4, {0.1, -0.2, 0.3});
  rod.setState(rod.initialState({Vec3(1, 0.2, 0), Vec3(1, 0, 0.1), Vec3(1, -0.1, 0)},
                                {Vec3(0.1, 0.2, 0.3)}));
  ModelState checkpoint;
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
  EXPECT_ROD_ERROR(restored.unpack_state(checkpoint));
  EXPECT_TRUE(restored.state().isApprox(before, 0));
}

TEST(BendingRodImplicitTest, NonconvergenceDoesNotAdvanceOrChangeHistory)
{
  BendingRodImplicitConfig config;
  config.links = 3;
  config.ei = 1;
  config.max_iterations = 1;
  config.nonlinear_tolerance = 1e-13;
  TestIBRod rod(config);
  rod.setState(rod.initialState({Vec3(1, 0.2, 0), Vec3(1, 0, 0.1), Vec3(1, -0.1, 0)}));
  ModelState before, after;
  rod.pack_state(before);
  std::vector<IBMesh::IBVertex> force(4, {0, 0, 0});
  EXPECT_ROD_ERROR(rod.GetNext(force, 0.1));
  rod.pack_state(after);
  EXPECT_EQ(before.ints, after.ints);
  EXPECT_EQ(before.reals, after.reals);
}

TEST(BendingRodImplicitTest, RejectsInvalidConfigurationAndRestartsHistory)
{
  BendingRodImplicitConfig config;
  config.links = 2;
  config.spectral_radius = -0.1;
  EXPECT_ROD_ERROR(TestIBRod rod(config));
  config.spectral_radius = 1.1;
  EXPECT_ROD_ERROR(TestIBRod rod(config));
  config.spectral_radius = 0.5;
  config.max_iterations = 0;
  EXPECT_ROD_ERROR(TestIBRod rod(config));
  config.max_iterations = 20;
  TestIBRod rod(config), fresh(config);
  std::vector<IBMesh::IBVertex> force(3, {0, 0, 0});
  rod.GetNext(force, 0.01);
  // Reset through the adapter must also discard algorithmic history.
  rod.setState(fresh.state());
  rod.GetNext(force, 0.01);
  fresh.GetNext(force, 0.01);
  EXPECT_TRUE(rod.state().isApprox(fresh.state(), 1e-14));
  const auto before = rod.state();
  EXPECT_ROD_ERROR(rod.GetNext(force, -0.01));
  EXPECT_TRUE(rod.state().isApprox(before, 0));
}

#undef EXPECT_ROD_ERROR

} // namespace
} // namespace ELFF
