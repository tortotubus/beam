

#include <elff/io/CXX/vtkHDFPolyData.hpp>
#include <elff/models/beam/EulerBeamGGL.hpp>
#include <elff/models/beam/EulerBeamMoM.hpp>
#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <string>

#include "EulerBeamGlowinskiCsv.hpp"

namespace ELFF {

using namespace IO::CXX;
using namespace Models;

TEST(EulerBeamGGLTest, Glowinski)
{
  GTEST_LOG_(INFO) << "CTEST_FULL_OUTPUT";

  const real_t length = 32.6;
  const real_t EI = 700.;
  const real_t mu = 7.67;
  const real_t r_penalty = 1e7;
  const std::array<real_t, 3> load = { 0., -9.81 * mu, 0. };

  const real_t dt = 1e-2;
  const real_t tf = 10.;
  const size_t Nt = static_cast<size_t>(std::ceil(tf / dt));
  const size_t nodes = 61;

  EulerBeam::EulerBeamBCs boundary_conditions = {
    .end = { EulerBeam::left, EulerBeam::right },
    .type = { EulerBeam::simple_bc, EulerBeam::simple_bc },
    .vals = { {
                .position = { 0., 0., 0. },
              },
              {
                .position = { 20., 0., 0. },
              } }
  };

  EulerBeamMoM static_beam(
    length, EI, nodes, boundary_conditions, r_penalty);
  static_beam.apply_initial_condition();

  ELFF_LOG("Static Solve:");
  static_beam.solve(load);

  boundary_conditions.type[1] = EulerBeam::free_bc;

  EulerBeamGGL dynamic_beam(
    length, EI, mu, nodes, boundary_conditions, r_penalty);
  dynamic_beam.apply_initial_condition(static_beam.get_mesh());

  std::ofstream csv("glowinski_ggl.csv");
  ASSERT_TRUE(csv.good());
  Benchmark::write_glowinski_csv_header(csv);

  std::ofstream snapshots("glowinski_ggl_centerlines_5_10.csv");
  std::ofstream gp("glowinski_ggl_5_10.gp");
  ASSERT_TRUE(snapshots.good());
  ASSERT_TRUE(gp.good());
  snapshots << std::setprecision(17) << "# time,x,y,z\n";
  gp << R"gnuplot(if (!exists("csv_file")) csv_file = "glowinski_ggl_centerlines_5_10.csv"
if (!exists("plot_output")) plot_output = "glowinski_ggl_5_10.tex"
set datafile separator comma
set terminal cairolatex pdf size 4.8in,3.2in color colortext font ",10"
set output plot_output
set xlabel "$x$"
set ylabel "$y$"
set size ratio -1
set grid
unset key
set cblabel "$t$ (s)"
set cbrange [5:10]
plot for [i=0:20] csv_file index i using 2:3:1 with lines lc palette
)gnuplot";
  gp.close();
  ASSERT_TRUE(gp.good());
  const size_t first_snapshot = static_cast<size_t>(std::lround(5.0 / dt));
  const size_t snapshot_stride = static_cast<size_t>(std::lround(0.25 / dt));
  size_t snapshot_count = 0;

  ELFF_LOG("Dynamic Solve:");
  for (size_t ti = 0; ti <= Nt; ++ti) {
    const std::string filename = "glowinski_ggl.vtkhdf";
    vtkPolyData pd = dynamic_beam.get_mesh().to_vtk_polydata();
    vtkHDFPolyData hdf_pd(filename, pd);

    if (ti == 0) {
      hdf_pd.write_new_transient(true, ti * dt);
    } else {
      hdf_pd.append_transient(ti * dt);
    }

    Benchmark::write_glowinski_csv_frame(
      csv, ti, ti * dt, dynamic_beam.get_mesh());
    if (ti >= first_snapshot && (ti - first_snapshot) % snapshot_stride == 0) {
      for (const auto &point : dynamic_beam.get_mesh().get_centerline())
        snapshots << ti * dt << ',' << point[0] << ',' << point[1]
                  << ',' << point[2] << '\n';
      snapshots << "\n\n";
      ++snapshot_count;
    }
    if (ti < Nt) dynamic_beam.solve(dt, load);
  }
  snapshots.close();
  ASSERT_TRUE(snapshots.good());
  EXPECT_EQ(snapshot_count, 21u);
}

TEST(EulerBeamGGLTest, Huang)
{
  GTEST_LOG_(INFO) << "CTEST_FULL_OUTPUT";

  const real_t length = 1.;
  const size_t nodes = 30;
  const real_t kappa = 0.1 * M_PI;

  EulerBeamMesh ic_mesh(nodes, length);

  std::vector<std::array<real_t, 3>>& ic_centerline = ic_mesh.get_centerline();
  std::vector<std::array<real_t, 3>>& ic_slope = ic_mesh.get_slope();
  std::vector<std::array<real_t, 3>>& ic_velocity =
    ic_mesh.get_centerline_velocity();
  std::vector<real_t>& ic_s = ic_mesh.get_curvilinear_axis();

  for (size_t ni = 0; ni < nodes; ++ni) {
    const real_t s = ic_s[ni];

    ic_centerline[ni][0] = (length - s) * std::cos(kappa);
    ic_centerline[ni][1] = (length - s) * std::sin(kappa);
    ic_centerline[ni][2] = 0.;

    ic_slope[ni][0] = -std::cos(kappa);
    ic_slope[ni][1] = -std::sin(kappa);
    ic_slope[ni][2] = 0.;

    ic_velocity[ni][0] = 0.;
    ic_velocity[ni][1] = 0.;
    ic_velocity[ni][2] = 0.;
  }

  EulerBeam::EulerBeamBCs boundary_conditions = {
    .end = { EulerBeam::left, EulerBeam::right },
    .type = { EulerBeam::free_bc, EulerBeam::simple_bc },
    .vals = { {
                .position = { 0., 0., 0. },
              },
              {
              } }
  };

  const real_t EI = 0.01;
  const real_t mu = 1.;
  const real_t r_penalty = 1e4;
  const real_t dt = 0.02;
  const real_t tf = 0.8;
  const size_t Nt = static_cast<size_t>(std::ceil(tf / dt));
  const std::array<real_t, 3> load = { 10., 0., 0. };

  EulerBeamGGL dynamic_beam(
    length, EI, mu, nodes, boundary_conditions, r_penalty);
  dynamic_beam.apply_initial_condition(ic_mesh);

  for (size_t ti = 0; ti < Nt; ++ti) {
    const std::string filename = "huang_ggl.vtkhdf";
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
