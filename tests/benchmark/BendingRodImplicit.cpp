#include <elff/io/CXX/vtkHDFPolyData.hpp>
#include <elff/models/articulated/IBBendingRodImplicit.hpp>
#include <gtest/gtest.h>
#include <Eigen/Dense>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <vector>

namespace ELFF {
namespace {
using namespace Models;

IO::CXX::vtkPolyData rod_polydata(const IBMesh &mesh)
{
  IO::CXX::vtkPolyData data;
  for (const auto &point : mesh.GetPoints())
    data.add_point(point.x, point.y, point.z);
  for (std::size_t i = 1; i < mesh.GetNumberOfPoints(); ++i)
    data.add_line(i - 1, i);
  const auto field = data.add_pointdata_vector("velocity", 3);
  for (std::size_t i = 0; i < mesh.GetNumberOfPoints(); ++i) {
    const auto &v = mesh.GetVelocity()[i];
    data.set_pointdata_vector3(field, i, {v.x, v.y, v.z});
  }
  return data;
}

// Maximum absolute error in the individual segment lengths.
double length_error(const IBMesh &mesh, double segment_length)
{
  double error = 0.0;
  const auto &points = mesh.GetPoints();
  for (std::size_t i = 1; i < points.size(); ++i) {
    const Vec3 delta(points[i].x - points[i - 1].x,
                     points[i].y - points[i - 1].y,
                     points[i].z - points[i - 1].z);
    error = std::max(error, std::abs(delta.norm() - segment_length));
  }
  return error;
}

// Benchmark-only planar equilibrium: minimize bending plus gravitational
// energy with both endpoints pinned. Angles enforce segment lengths exactly;
// two Lagrange multipliers enforce the horizontal and vertical endpoint span.
// This initializer is specific to uniform Glowinski material and vertical gravity.
::testing::AssertionResult static_glowinski_shape(const BendingRodImplicitConfig &config,
                           double span, BendingRodMesh &mesh)
{
  const int n = config.links;
  const double h = config.length / n;
  const double mass = config.total_mass / n;
  const double gravity = 9.81;
  const double scale = config.total_mass * gravity * config.length;
  const double stiffness = config.ei / h / scale;
  Eigen::VectorXd state = Eigen::VectorXd::Zero(n + 2);

  // A symmetric, downward circular arc provides a feasible initial guess.
  double low = 0, high = 2 * std::acos(-1.0) / n;
  for (int iteration = 0; iteration < 80; ++iteration) {
    const double delta = (low + high) / 2;
    const double chord = h * std::sin(n * delta / 2) / std::sin(delta / 2);
    if (chord > span) low = delta;
    else high = delta;
  }
  for (int i = 0; i < n; ++i)
    state[i] = (i - (n - 1) / 2.0) * (low + high) / 2;

  // Normalize energy and endpoint constraints to keep the KKT system scaled.
  auto evaluate = [&](const Eigen::VectorXd &x, Eigen::VectorXd &residual,
                      Eigen::MatrixXd &jacobian) {
    residual = Eigen::VectorXd::Zero(n + 2);
    jacobian = Eigen::MatrixXd::Zero(n + 2, n + 2);
    residual[n] = -span / config.length;
    for (int i = 0; i < n; ++i) {
      const double sine = std::sin(x[i]), cosine = std::cos(x[i]);
      const double weight = mass * gravity * h * (n - i - 0.5) / scale;
      const double dx = -sine / n, dy = cosine / n;
      residual[i] = weight * cosine + x[n] * dx + x[n + 1] * dy;
      residual[n] += cosine / n;
      residual[n + 1] += sine / n;
      jacobian(i, i) = -weight * sine - (x[n] * cosine + x[n + 1] * sine) / n;
      jacobian(i, n) = jacobian(n, i) = dx;
      jacobian(i, n + 1) = jacobian(n + 1, i) = dy;
      if (i > 0) {
        const double bend = stiffness * (x[i] - x[i - 1]);
        residual[i] += bend;
        residual[i - 1] -= bend;
        jacobian(i, i) += stiffness;
        jacobian(i - 1, i - 1) += stiffness;
        jacobian(i, i - 1) -= stiffness;
        jacobian(i - 1, i) -= stiffness;
      }
    }
  };

  Eigen::VectorXd residual;
  Eigen::MatrixXd jacobian;
  evaluate(state, residual, jacobian);
  const Eigen::MatrixXd initial_constraints = jacobian.bottomLeftCorner(2, n);
  state.tail<2>() = (initial_constraints * initial_constraints.transpose()).ldlt().solve(
      -initial_constraints * residual.head(n));
  evaluate(state, residual, jacobian);
  for (int iteration = 0; residual.lpNorm<Eigen::Infinity>() > 1e-11; ++iteration) {
    if (iteration == 100) return ::testing::AssertionFailure() << "Static iteration limit; residual=" << residual.norm();
    const Eigen::VectorXd correction = jacobian.fullPivLu().solve(-residual);
    if (!correction.allFinite()) return ::testing::AssertionFailure() << "Nonfinite static correction";
    bool accepted = false;
    for (double alpha = 1; alpha >= 1.0 / 65536; alpha *= 0.5) {
      Eigen::VectorXd trial_residual;
      Eigen::MatrixXd trial_jacobian;
      const Eigen::VectorXd trial = state + alpha * correction;
      evaluate(trial, trial_residual, trial_jacobian);
      if (trial_residual.allFinite() &&
          trial_residual.norm() < (1 - 1e-4 * alpha) * residual.norm()) {
        state = trial;
        residual = std::move(trial_residual);
        jacobian = std::move(trial_jacobian);
        accepted = true;
        break;
      }
    }
    if (!accepted) return ::testing::AssertionFailure() << "Static line search failed; residual=" << residual.norm();
  }

  // Check that this is a constrained local minimum, not merely a stationary
  // configuration. Also stay on the principal joint-angle branch used by the rod.
  const Eigen::MatrixXd constraints = jacobian.bottomLeftCorner(2, n);
  const Eigen::MatrixXd nullspace = constraints.fullPivLu().kernel();
  const Eigen::MatrixXd reduced =
      nullspace.transpose() * jacobian.topLeftCorner(n, n) * nullspace;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> stability(reduced);
  if (stability.info() != Eigen::Success || stability.eigenvalues().minCoeff() <= 0)
    return ::testing::AssertionFailure() << "Static equilibrium is not a stable minimum";
  for (int i = 1; i < n; ++i)
    if (std::abs(state[i] - state[i - 1]) >= std::acos(-1.0))
      return ::testing::AssertionFailure() << "Joint angle exceeds principal branch";

  mesh.get_centerline(0) = {config.pin.x(), config.pin.y(), config.pin.z()};
  for (int i = 0; i < n; ++i) {
    const auto &previous = mesh.get_centerline(i);
    mesh.get_centerline(i + 1) =
        {previous[0] + h * std::cos(state[i]),
         previous[1] + h * std::sin(state[i]), previous[2]};
  }
  for (auto &velocity : mesh.get_centerline_velocity()) velocity = {0, 0, 0};
  ::testing::Test::RecordProperty("static_residual", residual.lpNorm<Eigen::Infinity>());
  return ::testing::AssertionSuccess();
}

} // namespace

void run_glowinski(BendingRodImplicitConfig config, const std::string &output_name)
{
  // Match the Glowinski beam benchmarks: 61 nodes, length 32.6, EI 700,
  // linear density 7.67 and gravity 9.81; run for ten seconds.
  constexpr int nodes = 61;
  constexpr double linear_density = 7.67;
  constexpr double dt = 1e-2;
  constexpr double final_time = 10.0;
  const int steps = static_cast<int>(std::ceil(final_time / dt));
  const int output_stride = static_cast<int>(std::lround(0.01 / dt));

  config.links = nodes - 1;
  config.length = 32.6;
  config.total_mass = linear_density * config.length;
  config.ei = 700.0;
  // The first endpoint is pinned at the origin, with free rotation; the
  // last endpoint is free. No added joint damping, as in the Glowinski cases.
  config.boundary = Boundary::Pinned;
  config.pin = Vec3::Zero();
  IBBendingRodImplicit rod(config);

  // Equilibrate with both endpoints supported, then release the second support.
  // The dynamic rod retains only the pin at node zero.
  BendingRodMesh initial(nodes, config.length);
  ASSERT_TRUE(static_glowinski_shape(config, 20.0, initial))
      << "Static Glowinski equilibrium did not converge to a local minimum";
  ASSERT_NEAR(initial.get_centerline(nodes - 1)[0], 20.0, 1e-8);
  ASSERT_NEAR(initial.get_centerline(nodes - 1)[1], 0.0, 1e-8);
  ASSERT_LT(initial.get_centerline(nodes / 2)[1], -1.0);
  for (int i = 0; i < config.links; ++i) {
    const auto &a = initial.get_centerline(i);
    const auto &b = initial.get_centerline(i + 1);
    ASSERT_NEAR(Vec3(b[0] - a[0], b[1] - a[1], b[2] - a[2]).norm(),
                config.length / config.links, 1e-12);
  }
  rod.apply_initial_condition(initial);
  const auto initial_tip = initial.get_centerline(nodes - 1);
  const double geometry_tolerance = std::max(1e-12, 20.0 * std::numeric_limits<real_t>::epsilon());

  // Supply weight as force density; the IB adapter applies nodal measures.
  std::vector<IBMesh::IBVertex> traction(config.links + 1, {0, -9.81 * linear_density, 0});
  IBForceCoupled &coupled = rod;
  std::ofstream csv(output_name + ".csv");
  std::ofstream centerlines(output_name + "_centerlines_5_10.csv");
  std::ofstream gp(output_name + "_5_10.gp");
  ASSERT_TRUE(csv.good());
  ASSERT_TRUE(centerlines.good());
  ASSERT_TRUE(gp.good());
  csv << std::setprecision(17)
      << "time,tip_x,tip_y,tip_z,energy,max_length_error\n";
  centerlines << std::setprecision(17) << "# time,x,y,z\n";
  // Capture [5, 10] seconds, including both endpoints, at 0.25-second intervals.
  constexpr int snapshots = 21;
  const int first_snapshot_step = static_cast<int>(std::lround(5.0 / dt));
  int next_snapshot = 0;
  gp << "if (!exists(\"csv_file\")) csv_file = \""
     << output_name << "_centerlines_5_10.csv\"\n"
     << "if (!exists(\"plot_output\")) plot_output = \""
     << output_name << "_5_10.tex\"\n"
     << R"gnuplot(set datafile separator comma
set terminal cairolatex pdf size 4.8in,3.2in color colortext font ",10"
set output plot_output
set xlabel "$x$"
set ylabel "$y$"
set size ratio -1
set grid
unset key
set cblabel "$t$ (s)"
)gnuplot"
     << "set cbrange [5:" << final_time << "]\n"
     << "plot for [i=0:" << snapshots - 1
     << "] csv_file index i using 2:3:1 with lines lc palette\n";
  gp.close();
  ASSERT_TRUE(gp.good());

  const auto start = std::chrono::steady_clock::now();
  for (int step = 0; step <= steps; ++step) {
    SCOPED_TRACE(::testing::Message() << "step=" << step);
    if (step > 0)
      coupled.GetNext(traction, dt);
    for (const auto &point : rod.get_mesh().get_centerline())
      for (auto value : point) ASSERT_TRUE(std::isfinite(value));
    if (next_snapshot < snapshots &&
        step == first_snapshot_step + static_cast<int>(std::lround(
                    next_snapshot * static_cast<double>(steps - first_snapshot_step) / (snapshots - 1)))) {
      for (const auto &point : rod.get_mesh().get_centerline())
        centerlines << step * dt << ',' << point[0] << ',' << point[1]
                    << ',' << point[2] << '\n';
      centerlines << "\n\n";
      ++next_snapshot;
    }
    if (step % output_stride != 0 && step != steps)
      continue;

    const auto &mesh = coupled.GetCurrent();
    const auto &pin = mesh.GetPoints().front();
    const auto &pin_velocity = mesh.GetVelocity().front();
    ASSERT_NEAR(pin.x, config.pin.x(), geometry_tolerance);
    ASSERT_NEAR(pin.y, config.pin.y(), geometry_tolerance);
    ASSERT_NEAR(pin.z, config.pin.z(), geometry_tolerance);
    ASSERT_DOUBLE_EQ(pin_velocity.x, 0.0);
    ASSERT_DOUBLE_EQ(pin_velocity.y, 0.0);
    ASSERT_DOUBLE_EQ(pin_velocity.z, 0.0);
    const double error = length_error(mesh, config.length / config.links);
    ASSERT_LT(error, geometry_tolerance);

    double energy = rod.energy();
    // Include the external weight potential in this benchmark diagnostic.
    for (int i = 0; i <= config.links; ++i) {
      const double measure = config.length / config.links *
                             ((i == 0 || i == config.links) ? 0.5 : 1.0);
      energy += linear_density * 9.81 * measure * mesh.GetPoints()[i].y;
    }
    ASSERT_TRUE(std::isfinite(energy));
    const auto &tip = mesh.GetPoints().back();
    csv << step * dt << ',' << tip.x << ',' << tip.y << ',' << tip.z << ','
        << energy << ',' << error << '\n';
    {
      auto data = rod_polydata(mesh);
      IO::CXX::vtkHDFPolyData output(output_name + ".vtkhdf", data);
      if (step == 0)
        output.write_new_transient(true, 0.0);
      else
        output.append_transient(step * dt);
    }
  }
  csv.flush();
  ASSERT_TRUE(csv.good());
  centerlines.close();
  ASSERT_TRUE(centerlines.good());
  EXPECT_EQ(next_snapshot, snapshots);

  const auto final_tip = rod.get_mesh().get_centerline(nodes - 1);
  double displacement_squared = 0;
  for (size_t d = 0; d < 3; ++d)
    displacement_squared += std::pow(final_tip[d] - initial_tip[d], 2);
  EXPECT_GT(std::sqrt(displacement_squared), 0.05 * config.length);

  const double elapsed = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - start).count();
  ::testing::Test::RecordProperty("elapsed_seconds", elapsed);
  GTEST_LOG_(INFO) << output_name << " (output dt " << dt << "): "
                   << elapsed << " s";
}

TEST(BendingRodImplicitBenchmark, GlowinskiSaggedPinnedFree)
{
  // One fully implicit generalized-alpha step per output interval.
  run_glowinski(BendingRodImplicitConfig{},
      "glowinski_bending_rod_implicit_sagged");
}

} // namespace ELFF
