#include <elff/io/CXX/vtkHDFPolyData.hpp>
#include <elff/models/beam/EulerBeamInextensiblePenalty.hpp>
#include <gtest/gtest.h>

#include "EulerBeamStaticInextensibleReferences.hpp"

#include <cmath>
#include <string>

namespace ELFF {

using namespace IO::CXX;
using namespace Models;

TEST(EulerBeamInextensiblePenaltyTest, BisshoppAndDrucker)
{
  const real_t length = 1.0;
  const real_t EI = 1.0;
  const real_t r_penalty = 1e5;
  const size_t nodes = 512;
  const real_t tip_force_y = -1.0;
  const double comparison_tol = 1e-6;

  EulerBeam::EulerBeamBCs boundary_conditions = {
    .end = { EulerBeam::left, EulerBeam::right },
    .type = { EulerBeam::clamped_bc, EulerBeam::point_force_bc },
    .vals = { { .position = { 0, 0, 0 }, .slope = { 1, 0, 0 } },
              { .force = { 0, tip_force_y, 0 } } }
  };

  EulerBeamInextensiblePenalty static_beam(
    length, EI, nodes, boundary_conditions, r_penalty);

  static_beam.solve();
  static_beam.get_mesh().plot_gnuplot("Bisshopp and Drucker Penalty");

  EulerBeamMesh& mesh = static_beam.get_mesh();
  const auto centerline = mesh.get_centerline();
  const std::array<real_t, 3> tip = centerline[nodes - 1];

  const BisshoppAndDrucker1945Result res =
    BisshoppAndDrucker1945(length, EI, -tip_force_y);

  EXPECT_NEAR(std::abs(length - tip[0]), res.A, comparison_tol);
  EXPECT_NEAR(std::abs(tip[1]), res.delta, comparison_tol);
}

TEST(EulerBeamInextensiblePenaltyTest, Glowinski)
{
  GTEST_LOG_(INFO) << "CTEST_FULL_OUTPUT";

  const real_t length = 32.6;
  const real_t EI = 700.0;
  const real_t mu = 7.67;
  const real_t r_penalty = 1e6;
  const std::array<real_t, 3> load = { 0.0, -9.81 * mu, 0.0 };

  const real_t dt = 1e-2;
  const real_t tf = 10.0;
  const size_t Nt = static_cast<size_t>(std::ceil(tf / dt));
  const size_t nodes = 61;

  EulerBeam::EulerBeamBCs boundary_conditions = {
    .end = { EulerBeam::left, EulerBeam::right },
    .type = { EulerBeam::simple_bc, EulerBeam::simple_bc },
    .vals = { {
                .position = { 0, 0, 0 },
              },
              {
                .position = { 20, 0, 0 },
              } }
  };

  EulerBeamInextensiblePenalty static_beam(
    length, EI, nodes, boundary_conditions, r_penalty);
  static_beam.apply_initial_condition();

  ELFF_LOG("Static Solve:");
  static_beam.solve(load);

  boundary_conditions.type[1] = EulerBeam::free_bc;

  EulerBeamInextensiblePenalty dynamic_beam(
    length, EI, mu, nodes, boundary_conditions, r_penalty);
  dynamic_beam.apply_initial_condition(static_beam.get_mesh());

  ELFF_LOG("Dynamic Solve:");
  for (size_t ti = 0; ti < Nt; ++ti) {
    const std::string filename = "glowinski_penalty.vtkhdf";
    vtkPolyData pd = dynamic_beam.get_mesh().to_vtk_polydata();
    vtkHDFPolyData hdf_pd(filename, pd);

    if (ti == 0) {
      hdf_pd.write_new_transient(true, ti * dt);
    } else {
      hdf_pd.append_transient(ti * dt);
    }

    dynamic_beam.solve(dt, load);
  }
}

TEST(EulerBeamInextensiblePenaltyTest, Huang)
{
  GTEST_LOG_(INFO) << "CTEST_FULL_OUTPUT";

  const real_t length = 1.0;
  const size_t nodes = 30;
  const real_t kappa = 0.1 * M_PI;

  EulerBeamMesh ic_mesh(nodes, length);

  auto& ic_centerline = ic_mesh.get_centerline();
  auto& ic_slope = ic_mesh.get_slope();
  auto& ic_velocity = ic_mesh.get_centerline_velocity();
  auto& ic_s = ic_mesh.get_curvilinear_axis();

  for (size_t ni = 0; ni < nodes; ++ni) {
    const real_t s = ic_s[ni];

    ic_centerline[ni][0] = (length - s) * std::cos(kappa);
    ic_centerline[ni][1] = (length - s) * std::sin(kappa);
    ic_centerline[ni][2] = 0.0;

    ic_slope[ni][0] = -std::cos(kappa);
    ic_slope[ni][1] = -std::sin(kappa);
    ic_slope[ni][2] = 0.0;

    ic_velocity[ni][0] = 0.0;
    ic_velocity[ni][1] = 0.0;
    ic_velocity[ni][2] = 0.0;
  }

  EulerBeam::EulerBeamBCs boundary_conditions = {
    .end = { EulerBeam::left, EulerBeam::right },
    .type = { EulerBeam::free_bc, EulerBeam::simple_bc },
    .vals = { {
                .position = { 0, 0, 0 },
              },
              {
              } }
  };

  const real_t EI = 0.01;
  const real_t mu = 1.0;
  const real_t r_penalty = 1e5;
  const real_t dt = 0.02;
  const real_t tf = 0.8;
  const size_t Nt = static_cast<size_t>(std::ceil(tf / dt));
  const std::array<real_t, 3> load = { 10.0, 0.0, 0.0 };

  EulerBeamInextensiblePenalty dynamic_beam(
    length, EI, mu, nodes, boundary_conditions, r_penalty);
  dynamic_beam.apply_initial_condition(ic_mesh);

  for (size_t ti = 0; ti < Nt; ++ti) {
    const std::string filename = "huang_penalty.vtkhdf";
    vtkPolyData pd = dynamic_beam.get_mesh().to_vtk_polydata();
    vtkHDFPolyData hdf_pd(filename, pd);

    if (ti == 0) {
      hdf_pd.write_new_transient(true, ti * dt);
    } else {
      hdf_pd.append_transient(ti * dt);
    }

    dynamic_beam.solve(dt, load);
  }
}

} // namespace ELFF
