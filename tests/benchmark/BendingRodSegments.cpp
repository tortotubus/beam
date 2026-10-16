#include <elff/io/CXX/vtkHDFPolyData.hpp>
#include <elff/models/articulated/BendingRodSegmentsIMEX.hpp>
#include <elff/models/articulated/BendingRodSegmentsImplicit.hpp>
#include <gtest/gtest.h>

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

} // namespace

template <class Rod, class Config>
void run_glowinski_straight(Config config, const std::string &output_name)
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
  config.gravity = Vec3(0.0, -9.81, 0.0);
  Rod rod(config);

  // Start horizontal and at rest, rather than from the static sagged shape
  // used by the other Glowinski benchmarks.
  const Eigen::VectorXd initial = rod.initialState({Vec3::UnitX()});
  rod.setState(initial);
  const Vec3 initial_tip = rod.nodes(initial).row(config.links).transpose();
  const double geometry_tolerance = std::max(1e-12, 20.0 * std::numeric_limits<real_t>::epsilon());

  // Gravity is already in the model; no additional IB traction is applied.
  std::vector<IBMesh::IBVertex> traction(config.links + 1, {0, 0, 0});
  IBForceCoupled &coupled = rod;
  std::ofstream csv(output_name + ".csv");
  ASSERT_TRUE(csv.good());
  csv << std::setprecision(17)
      << "time,tip_x,tip_y,tip_z,energy,max_length_error\n";

  const auto start = std::chrono::steady_clock::now();
  for (int step = 0; step <= steps; ++step) {
    SCOPED_TRACE(::testing::Message() << "step=" << step);
    if (step > 0)
      coupled.GetNext(traction, dt);
    ASSERT_TRUE(rod.state().allFinite());
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

    const double energy = rod.energy(rod.state());
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

  const Vec3 final_tip = rod.nodes(rod.state()).row(config.links).transpose();
  EXPECT_GT((final_tip - initial_tip).norm(), 0.05 * config.length);

  const double elapsed = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - start).count();
  ::testing::Test::RecordProperty("elapsed_seconds", elapsed);
  GTEST_LOG_(INFO) << output_name << " (output dt " << dt << "): "
                   << elapsed << " s";
}

TEST(BendingRodSegmentsBenchmark, GlowinskiStraightPinnedFree)
{
  BendingRodSegmentsConfig config;
  config.integration_step = 0.01 / 8;
  run_glowinski_straight<BendingRodSegments>(config,
      "glowinski_bending_rod_segments_straight");
}

TEST(BendingRodSegmentsIMEXBenchmark, GlowinskiStraightPinnedFree)
{
  // One second-order IMEX update every 0.01 s; no adaptive step control.
  run_glowinski_straight<BendingRodSegmentsIMEX>(BendingRodSegmentsIMEXConfig{},
      "glowinski_bending_rod_segments_imex_straight");
}

TEST(BendingRodSegmentsImplicitBenchmark, GlowinskiStraightPinnedFree)
{
  // One fully implicit generalized-alpha step per output interval.
  run_glowinski_straight<BendingRodSegmentsImplicit>(BendingRodSegmentsImplicitConfig{},
      "glowinski_bending_rod_segments_implicit_straight");
}

} // namespace ELFF
