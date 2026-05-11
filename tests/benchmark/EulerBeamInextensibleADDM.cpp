#include <cmath>
#include <cstdio>

#include <elff/io/CXX/vtkHDFPolyData.hpp>
#include <elff/models/beam/EulerBeamInextensibleADDM.hpp>
#include <elff/models/beam/EulerBeamInextensibleADDMTestHarness.hpp>
#include <elff/models/beam/EulerBeamInextensibleMoM.hpp>

#include <fstream>
#include <gtest/gtest.h>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "EulerBeamDynamicInextensibleReferences.hpp"
#include "EulerBeamStaticInextensibleReferences.hpp"

namespace ELFF {

using namespace IO::CXX;
using namespace Models;

std::vector<real_t> logspace(double xmin, double xmax, std::size_t N) {
  std::vector<real_t> values;
  values.reserve(N);

  if (N == 0) {
    return values;
  }

  if (N == 1) {
    values.push_back(xmin);
    return values;
  }

  const real_t log_min = std::log10(xmin);
  const real_t log_max = std::log10(xmax);

  for (std::size_t i = 0; i < N; ++i) {
    const real_t theta = static_cast<real_t>(i) / static_cast<real_t>(N - 1);
    const real_t exponent = (1.0 - theta) * log_min + theta * log_max;
    values.push_back(std::pow(10.0, exponent));
  }

  return values;
}

TEST(EulerBeamInextensibleADDMTest, NonlinearMMF1) {

  std::string csv_filename = "addm-nonlinear-convergence.csv";
  std::string tex_filename = "addm-nonlinear-convergence.tex";
  std::string gp_filename = "addm-nonlinear-convergence.gp";
  std::string filename = "addm-nonlinear-convergence.vtkhdf";

  std::ofstream gp(gp_filename);

  gp << "if (!exists(\"csv_file\")) csv_file = \"" << csv_filename << "\"\n";
  gp << "if (!exists(\"plot_output\")) plot_output = \"../vector/"
     << tex_filename << "\"\n";

  gp << R"gnuplot(
set datafile separator comma
set terminal cairolatex pdf size 4.8in,3.2in color colortext font ",10"
set output plot_output
)gnuplot";

  gp << R"gnuplot(
set xlabel "$\\log{(1/\\Delta t)}$"
set ylabel "$\\log{\\lVert \\vec x(T) - \\vec x_h(T)\\rVert}$"
set grid
set key outside right 
plot)gnuplot";

  std::ofstream csv(csv_filename);
  csv << std::setprecision(16);
  csv << "dt,epsilon,L2_error,log_inv_dt,log_error\n";

  std::vector<real_t> stopping_tolerances = {1e-7, 1e-8, 1e-9};
  // std::vector<real_t> dt_values = {0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625};
  std::vector<real_t> dt_values = logspace(0.00625, 0.05, 10);

  real_t length = M_PI_2;
  size_t nodes = 451;
  EulerBeamMesh mmf_mesh(nodes, length);
  real_t ds = mmf_mesh.get_ds();

  real_t EI = 1.;
  real_t mu = 1.;
  real_t r_penalty = 1e3;
  real_t tf = 1.0;

  int time_history_required = 3;

  for (int tol_idx = 0; tol_idx < stopping_tolerances.size(); tol_idx++) {

    std::vector<real_t> error_measure;

    for (int dt_idx = 0; dt_idx < dt_values.size(); dt_idx++) {
      real_t dt = dt_values[dt_idx];
      int Nt = static_cast<int>(std::ceil(tf / dt));

      std::vector<std::array<real_t, 3>> position_history_l;
      std::vector<std::array<real_t, 3>> slope_history_l;
      std::vector<std::array<real_t, 3>> velocity_history_l;
      std::vector<std::array<real_t, 3>> acceleration_history_l;
      std::vector<std::array<real_t, 3>> slope_velocity_history_l;
      std::vector<std::array<real_t, 3>> slope_acceleration_history_l;

      std::vector<std::array<real_t, 3>> position_history_r;
      std::vector<std::array<real_t, 3>> slope_history_r;
      std::vector<std::array<real_t, 3>> velocity_history_r;
      std::vector<std::array<real_t, 3>> acceleration_history_r;
      std::vector<std::array<real_t, 3>> slope_velocity_history_r;
      std::vector<std::array<real_t, 3>> slope_acceleration_history_r;

      // real_t epsilon = stopping_tolerances[tol_idx] * std::pow(dt,2);
      real_t epsilon = stopping_tolerances[tol_idx];

      for (int time_idx = -time_history_required + 1; time_idx <= Nt;
           time_idx++) {
        mmf_mesh = ManufacturedDynamicResult1(nodes, time_idx * dt);
        auto &current_position = mmf_mesh.get_centerline();
        auto &current_slope = mmf_mesh.get_slope();
        auto &current_velocity = mmf_mesh.get_centerline_velocity();
        auto &current_acceleration = mmf_mesh.get_centerline_acceleration();
        size_t l_idx = 0;
        size_t r_idx = mmf_mesh.get_nodes() - 1;
        const real_t s_left = mmf_mesh.get_curvilinear_axis(l_idx);
        const real_t s_right = mmf_mesh.get_curvilinear_axis(r_idx);
        position_history_l.push_back(current_position[l_idx]);
        slope_history_l.push_back(current_slope[l_idx]);
        velocity_history_l.push_back(current_velocity[l_idx]);
        acceleration_history_l.push_back(current_acceleration[l_idx]);
        slope_velocity_history_l.push_back(
            ManufacturedDynamicResult1SlopeVelocity(s_left, time_idx * dt));
        slope_acceleration_history_l.push_back(
            ManufacturedDynamicResult1SlopeAcceleration(s_left, time_idx * dt));
        position_history_r.push_back(current_position[r_idx]);
        slope_history_r.push_back(current_slope[r_idx]);
        velocity_history_r.push_back(current_velocity[r_idx]);
        acceleration_history_r.push_back(current_acceleration[r_idx]);
        slope_velocity_history_r.push_back(
            ManufacturedDynamicResult1SlopeVelocity(s_right, time_idx * dt));
        slope_acceleration_history_r.push_back(
            ManufacturedDynamicResult1SlopeAcceleration(s_right,
                                                        time_idx * dt));
      }

      EulerBeam::EulerBeamTimeDependentBCs bcs_time_dependent = {
          .end = {EulerBeam::left, EulerBeam::right},
          .type = {EulerBeam::clamped_bc, EulerBeam::clamped_bc},
          .history = {{.position_history = position_history_l,
                       .slope_history = slope_history_l,
                       .velocity_history = velocity_history_l,
                       .acceleration_history = acceleration_history_l,
                       .slope_velocity_history = slope_velocity_history_l,
                       .slope_acceleration_history =
                           slope_acceleration_history_l},
                      {.position_history = position_history_r,
                       .slope_history = slope_history_r,
                       .velocity_history = velocity_history_r,
                       .acceleration_history = acceleration_history_r,
                       .slope_velocity_history = slope_velocity_history_r,
                       .slope_acceleration_history =
                           slope_acceleration_history_r}},
          .time_zero_idx = static_cast<size_t>(time_history_required - 1),
      };

      // EulerBeam::EulerBeamBCs bcs_static = {
      //   .end = {EulerBeam::left, EulerBeam::right},
      //   .type = {EulerBeam::clamped_bc, EulerBeam::free_bc},
      //   .vals = {
      //     {
      //       .slope = slope_history_l[2],
      //       .position = position_history_l[2],
      //     } , {}
      //   }
      // };

      mmf_mesh = ManufacturedDynamicResult1(nodes, 0.0);
      EulerBeamInextensibleADDMTestHarness beam(length, EI, mu, nodes,
                                                bcs_time_dependent, r_penalty);

      beam.apply_initial_condition(mmf_mesh);

      beam.set_outer_relaxation(1.0);
      beam.set_outer_tolerance(epsilon);
      beam.set_outer_iter_max(30000);

      ELFF_LOG("Epsilon = " << epsilon << ", dt = " << dt);

      for (int time_idx = 0; time_idx < Nt; time_idx++) {
        const auto linear_old = ManufacturedDynamicResult1LinearLoad(
            static_cast<int>(nodes), time_idx * dt, EI, mu);
        const auto linear_new = ManufacturedDynamicResult1LinearLoad(
            static_cast<int>(nodes), (time_idx + 1) * dt, EI, mu);
        const auto nonlinear_new = ManufacturedDynamicResult1NonlinearLoad(
            static_cast<int>(nodes), (time_idx + 1) * dt, EI, mu);

        std::vector<std::array<real_t, 3>> desired_averaged_load(nodes);
        for (size_t i = 0; i < nodes; ++i) {
          for (size_t d = 0; d < 3; ++d) {
            const real_t internal_new = nonlinear_new[i][d] - linear_new[i][d];
            desired_averaged_load[i][d] =
                0.5 * (linear_old[i][d] + linear_new[i][d]) + internal_new;
          }
        }

        beam.solve_averaged_load(dt, desired_averaged_load);
      }

      mmf_mesh = ManufacturedDynamicResult1(nodes, Nt * dt);
      real_t err2 = 0.;
      const auto &xh = beam.get_mesh().get_centerline();
      const auto &xex = mmf_mesh.get_centerline();

      for (size_t node_idx = 0; node_idx < nodes; ++node_idx) {
        real_t w = ds;
        if (node_idx == 0 || node_idx == nodes - 1) {
          w *= .5;
        }
        const real_t dx = xh[node_idx][0] - xex[node_idx][0];
        const real_t dy = xh[node_idx][1] - xex[node_idx][1];
        err2 += w * (dx * dx + dy * dy);
      }

      csv << dt << "," << stopping_tolerances[tol_idx] << "," << std::sqrt(err2)
          << "," << std::log10(1.0 / dt) << "," << std::log10(std::sqrt(err2))
          << "\n";
    }
    gp << " \\\n\tcsv_file every ::1 using ($2 == "
       << stopping_tolerances[tol_idx]
       << " ? $4 : 1/0):5 with linespoints title \"$\\\\epsilon = "
       << stopping_tolerances[tol_idx] << "$\",";
  }

  csv.close();
  gp.close();
}

TEST(EulerBeamInextensibleADDMTest, LinearMMF1) {

  std::string csv_filename = "addm-linear-convergence.csv";
  std::string tex_filename = "addm-linear-convergence.tex";
  std::string gp_filename = "addm-linear-convergence.gp";
  std::string filename = "addm-linear-convergence.vtkhdf";

  std::ofstream gp(gp_filename);

  gp << "if (!exists(\"csv_file\")) csv_file = \"" << csv_filename << "\"\n";
  gp << "if (!exists(\"plot_output\")) plot_output = \"../vector/"
     << tex_filename << "\"\n";

  gp << R"gnuplot(
set datafile separator comma
set terminal cairolatex pdf size 4.8in,3.2in color colortext font ",10"
set output plot_output
)gnuplot";

  gp << R"gnuplot(
set xlabel "$\\log{(1/\\Delta t)}$"
set ylabel "$\\log{\\lVert \\vec x(T) - \\vec x_h(T)\\rVert}$"
set grid
set key off 
plot)gnuplot";

  std::ofstream csv(csv_filename);
  csv << std::setprecision(16);
  csv << "dt,epsilon,L2_error,log_inv_dt,log_error\n";

  real_t stopping_tolerance = 1e-1;
  // std::vector<real_t> dt_values = logspace(.2, .00625, 6);
  std::vector<real_t> dt_values = logspace(0.00625, 0.05, 10);

  real_t length = M_PI_2;
  size_t nodes = 451;
  EulerBeamMesh mmf_mesh(nodes, length);
  real_t ds = mmf_mesh.get_ds();

  real_t EI = 1.;
  real_t mu = 1.;
  real_t r_penalty = 1;
  real_t tf = 1.0;

  int time_history_required = 3;

  std::vector<real_t> error_measure;

  for (int dt_idx = 0; dt_idx < dt_values.size(); dt_idx++) {
    real_t dt = dt_values[dt_idx];
    int Nt = static_cast<int>(std::ceil(tf / dt));

    std::vector<std::array<real_t, 3>> position_history_l;
    std::vector<std::array<real_t, 3>> slope_history_l;
    std::vector<std::array<real_t, 3>> velocity_history_l;
    std::vector<std::array<real_t, 3>> acceleration_history_l;
    std::vector<std::array<real_t, 3>> slope_velocity_history_l;
    std::vector<std::array<real_t, 3>> slope_acceleration_history_l;

    std::vector<std::array<real_t, 3>> position_history_r;
    std::vector<std::array<real_t, 3>> slope_history_r;
    std::vector<std::array<real_t, 3>> velocity_history_r;
    std::vector<std::array<real_t, 3>> acceleration_history_r;
    std::vector<std::array<real_t, 3>> slope_velocity_history_r;
    std::vector<std::array<real_t, 3>> slope_acceleration_history_r;

    // real_t epsilon = stopping_tolerances[tol_idx] * std::pow(dt,2);
    real_t epsilon = stopping_tolerance;

    for (int time_idx = -time_history_required + 1; time_idx <= Nt;
         time_idx++) {
      mmf_mesh = ManufacturedDynamicResult1(nodes, time_idx * dt);
      auto &current_position = mmf_mesh.get_centerline();
      auto &current_slope = mmf_mesh.get_slope();
      auto &current_velocity = mmf_mesh.get_centerline_velocity();
      auto &current_acceleration = mmf_mesh.get_centerline_acceleration();
      size_t l_idx = 0;
      size_t r_idx = mmf_mesh.get_nodes() - 1;
      const real_t s_left = mmf_mesh.get_curvilinear_axis(l_idx);
      const real_t s_right = mmf_mesh.get_curvilinear_axis(r_idx);
      position_history_l.push_back(current_position[l_idx]);
      slope_history_l.push_back(current_slope[l_idx]);
      velocity_history_l.push_back(current_velocity[l_idx]);
      acceleration_history_l.push_back(current_acceleration[l_idx]);
      slope_velocity_history_l.push_back(
          ManufacturedDynamicResult1SlopeVelocity(s_left, time_idx * dt));
      slope_acceleration_history_l.push_back(
          ManufacturedDynamicResult1SlopeAcceleration(s_left, time_idx * dt));
      position_history_r.push_back(current_position[r_idx]);
      slope_history_r.push_back(current_slope[r_idx]);
      velocity_history_r.push_back(current_velocity[r_idx]);
      acceleration_history_r.push_back(current_acceleration[r_idx]);
      slope_velocity_history_r.push_back(
          ManufacturedDynamicResult1SlopeVelocity(s_right, time_idx * dt));
      slope_acceleration_history_r.push_back(
          ManufacturedDynamicResult1SlopeAcceleration(s_right, time_idx * dt));
    }

    EulerBeam::EulerBeamTimeDependentBCs bcs_time_dependent = {
        .end = {EulerBeam::left, EulerBeam::right},
        .type = {EulerBeam::clamped_bc, EulerBeam::clamped_bc},
        .history = {{.position_history = position_history_l,
                     .slope_history = slope_history_l,
                     .velocity_history = velocity_history_l,
                     .acceleration_history = acceleration_history_l,
                     .slope_velocity_history = slope_velocity_history_l,
                     .slope_acceleration_history =
                         slope_acceleration_history_l},
                    {.position_history = position_history_r,
                     .slope_history = slope_history_r,
                     .velocity_history = velocity_history_r,
                     .acceleration_history = acceleration_history_r,
                     .slope_velocity_history = slope_velocity_history_r,
                     .slope_acceleration_history =
                         slope_acceleration_history_r}},
        .time_zero_idx = static_cast<size_t>(time_history_required - 1),
    };
 

    mmf_mesh = ManufacturedDynamicResult1(nodes, 0.0);
    EulerBeamInextensibleADDMTestHarness beam(length, EI, mu, nodes,
                                              bcs_time_dependent, r_penalty);

    beam.apply_initial_condition(mmf_mesh);

    beam.set_outer_relaxation(1.);
    beam.set_outer_tolerance(epsilon);
    beam.set_outer_iter_max(30000);

    ELFF_LOG("Epsilon = " << epsilon << ", dt = " << dt);

    for (int time_idx = 0; time_idx < Nt; time_idx++) {
      const auto linear_old = ManufacturedDynamicResult1LinearLoad(
          static_cast<int>(nodes), time_idx * dt, EI, mu);
      const auto linear_new = ManufacturedDynamicResult1LinearLoad(
          static_cast<int>(nodes), (time_idx + 1) * dt, EI, mu);

      std::vector<std::array<real_t, 3>> averaged_load(nodes);
      for (size_t i = 0; i < nodes; ++i) {
        for (size_t d = 0; d < 3; ++d) {
          averaged_load[i][d] = 0.5 * (linear_old[i][d] + linear_new[i][d]);
        }
      }

      beam.solve_linear(dt, averaged_load);
    }

    mmf_mesh = ManufacturedDynamicResult1(nodes, Nt * dt);
    real_t err2 = 0.;
    const auto &xh = beam.get_mesh().get_centerline();
    const auto &xex = mmf_mesh.get_centerline();

    for (size_t node_idx = 0; node_idx < nodes; ++node_idx) {
      real_t w = ds;
      if (node_idx == 0 || node_idx == nodes - 1) {
        w *= .5;
      }
      const real_t dx = xh[node_idx][0] - xex[node_idx][0];
      const real_t dy = xh[node_idx][1] - xex[node_idx][1];
      err2 += w * (dx * dx + dy * dy);
    }

    csv << dt << "," << stopping_tolerance << "," << std::sqrt(err2)
        << "," << std::log10(1.0 / dt) << "," << std::log10(std::sqrt(err2))
        << "\n";
  }

  gp << " \\\n\tcsv_file every ::1 using 4:5 with linespoints";

  csv.close();
  gp.close();
}

TEST(EulerBeamInextensibleADDMTest, BisshoppAndDrucker) {
  std::string csv_filename = "addm-static-convergence.csv";
  std::string gp_filename = "addm-static-convergence.gp";
  std::string tex_filename = "addm-static-convergence.tex";

  std::ofstream gp(gp_filename);
  std::ofstream csv(csv_filename);

  ASSERT_TRUE(gp.is_open());
  ASSERT_TRUE(csv.is_open());

  const real_t length = 1.;
  const real_t EI = 1.;
  const real_t r_penalty = 1e2;
  const real_t tip_force_y = -1.;
  const double finest_mesh_tol = 5e-7;

  EulerBeam::EulerBeamBCs boundary_conditions = {
      .end = {EulerBeam::left, EulerBeam::right},
      .type = {EulerBeam::clamped_bc, EulerBeam::point_force_bc},
      .vals = {{.position = {0, 0, 0}, .slope = {1, 0, 0}},
               {.force = {0, tip_force_y, 0}}}};

  BisshoppAndDrucker1945Result res =
      BisshoppAndDrucker1945(length, EI, -tip_force_y);

  gp << "if (!exists(\"csv_file\")) csv_file = \"" << csv_filename << "\"\n";
  gp << "if (!exists(\"plot_output\")) plot_output = \"../vector/"
     << tex_filename << "\"\n";

  gp << R"gnuplot(
set datafile separator comma
set terminal cairolatex pdf size 4.8in,3.2in color colortext font ",10"
set output plot_output

set xlabel "$\\log{(1/\\Delta s)}$"
set ylabel "$\\log{\\norm{\\vec {\\epsilon}_{\\rm tip}}}$"
set grid
set key off
plot \
  csv_file every ::1 using 6:7 with linespoints title "ADDM"
)gnuplot";

  csv << std::setprecision(16);
  csv << "nodes,ds,x_error,y_error,tip_error,log_inv_ds,log_error\n";

  const std::vector<size_t> node_counts = {11, 22, 44, 88, 172};
  real_t finest_tip_error = std::numeric_limits<real_t>::infinity();
  real_t previous_tip_error = std::numeric_limits<real_t>::infinity();

  for (const size_t nodes : node_counts) {
    EulerBeamInextensibleADDM beam(length, EI, nodes, boundary_conditions,
                                   r_penalty);
    beam.solve();

    EulerBeamMesh &mesh = beam.get_mesh();
    const auto &centerline = mesh.get_centerline();
    const std::array<real_t, 3> &tip = centerline[nodes - 1];

    const real_t x_error = std::abs(length - tip[0]) - res.A;
    const real_t y_error = std::abs(tip[1]) - res.delta;
    const real_t tip_error = std::sqrt(x_error * x_error + y_error * y_error);
    const real_t ds = mesh.get_ds();

    csv << nodes << "," << ds << "," << x_error << "," << y_error << ","
        << tip_error << "," << std::log10(1.0 / ds) << ","
        << std::log10(tip_error) << "\n";

    EXPECT_LT(tip_error, previous_tip_error);
    previous_tip_error = tip_error;
    finest_tip_error = tip_error;
  }

  EXPECT_LT(finest_tip_error, finest_mesh_tol);
}

TEST(EulerBeamInextensibleADDMTest, Glowinski) {
  GTEST_LOG_(INFO) << "CTEST_FULL_OUTPUT";
  real_t length = 32.6, EI = 700., mu = 7.67, r_penalty = 1e5;
  std::array<real_t, 3> load = {0, -9.81 * mu, 0};

  real_t dt = 1e-2;
  real_t t = 0;
  real_t tf = 10;
  size_t Nt = size_t(ceil(tf / dt));

  real_t dt_save = 0.1;
  size_t Nt_save = size_t(ceil(dt_save / dt));
  (void)t;
  (void)dt_save;
  (void)Nt_save;

  size_t nodes = 61;

  EulerBeam::EulerBeamBCs boundary_conditions = {
      .end = {EulerBeam::left, EulerBeam::right},
      .type = {EulerBeam::simple_bc, EulerBeam::simple_bc},
      .vals = {{
                   .position = {0, 0, 0},
               },
               {
                   .position = {20, 0, 0},
               }}};

  EulerBeamInextensibleMoM static_beam(length, EI, nodes, boundary_conditions,
                                       1e8);
  static_beam.apply_initial_condition();

  ELFF_LOG("Static Solve:");
  static_beam.solve(load);

  boundary_conditions.type[1] = EulerBeam::free_bc;

  EulerBeamInextensibleADDM beam(length, EI, mu, nodes, boundary_conditions,
                                 r_penalty);
  beam.apply_initial_condition(static_beam.get_mesh());
  beam.set_outer_tolerance(1e-10);

  ELFF_LOG("Dynamic Solve:");
  for (size_t ti = 0; ti < Nt; ti++) {
    std::string filename = "glowinski_addm.vtkhdf";

    if (ti == 0) {
      vtkPolyData pd = beam.get_mesh().to_vtk_polydata();
      vtkHDFPolyData hdf_pd(filename, pd);
      hdf_pd.write_new_transient(true, ti * dt);
    } else {
      vtkPolyData pd = beam.get_mesh().to_vtk_polydata();
      vtkHDFPolyData hdf_pd(filename, pd);
      hdf_pd.append_transient(ti * dt);
    }

    beam.solve(dt, load);
  }
}

TEST(EulerBeamInextensibleADDMTest, Huang) {
  GTEST_LOG_(INFO) << "CTEST_FULL_OUTPUT";

  real_t length = 1;
  size_t nodes = 30;
  real_t kappa = 0.1 * M_PI;

  real_t x0 = 0, y0 = 0;

  EulerBeamMesh ic_mesh(nodes, length);

  std::vector<std::array<real_t, 3>> &ic_centerline = ic_mesh.get_centerline();
  std::vector<std::array<real_t, 3>> &ic_slope = ic_mesh.get_slope();
  std::vector<std::array<real_t, 3>> &ic_velocity =
      ic_mesh.get_centerline_velocity();
  std::vector<real_t> &ic_s = ic_mesh.get_curvilinear_axis();

  for (size_t ni = 0; ni < nodes; ++ni) {
    real_t s = ic_s[ni];

    ic_centerline[ni][0] = x0 + (length - s) * std::cos(kappa);
    ic_centerline[ni][1] = y0 + (length - s) * std::sin(kappa);
    ic_centerline[ni][2] = 0.;

    ic_slope[ni][0] = -std::cos(kappa);
    ic_slope[ni][1] = -std::sin(kappa);
    ic_slope[ni][2] = 0.;

    ic_velocity[ni][0] = 0.;
    ic_velocity[ni][1] = 0.;
    ic_velocity[ni][2] = 0.;
  }

  EulerBeam::EulerBeamBCs boundary_conditions = {
      .end = {EulerBeam::left, EulerBeam::right},
      .type = {EulerBeam::free_bc, EulerBeam::simple_bc},
      .vals = {{
                   .position = {0, 0, 0},
               },
               {}}};

  real_t EI = 0.01;
  real_t mu = 1;
  real_t r_penalty = 1e3;

  real_t dt = 0.02;
  real_t t = 0;
  real_t tf = 0.8;
  size_t Nt = size_t(ceil(tf / dt));
  (void)t;

  std::array<real_t, 3> load = {10, 0, 0};

  EulerBeamInextensibleADDM beam(length, EI, mu, nodes, boundary_conditions,
                                 r_penalty);
  beam.apply_initial_condition(ic_mesh);

  for (size_t ti = 0; ti < Nt; ti++) {
    std::string filename = "huang_addm.vtkhdf";

    if (ti == 0) {
      vtkPolyData pd = beam.get_mesh().to_vtk_polydata();
      vtkHDFPolyData hdf_pd(filename, pd);
      hdf_pd.write_new_transient(true, ti * dt);
    } else {
      vtkPolyData pd = beam.get_mesh().to_vtk_polydata();
      vtkHDFPolyData hdf_pd(filename, pd);
      hdf_pd.append_transient(ti * dt);
    }

    beam.solve(dt, load);
  }
}

TEST(EulerBeamInextensibleADDMTest, DynamicBoundaryUpdateRefreshesHistory) {
  const real_t length = 1.0;
  const real_t EI = 1.0;
  const real_t mu = 1.0;
  const real_t r_penalty = 1e3;
  const real_t dt = 1e-2;
  const size_t nodes = 11;

  EulerBeam::EulerBeamBCs boundary_conditions = {
      .end = {EulerBeam::left, EulerBeam::right},
      .type = {EulerBeam::simple_bc, EulerBeam::simple_bc},
      .vals = {{.position = {0.0, 0.0, 0.0}},
               {.position = {length, 0.0, 0.0}}}};

  EulerBeamInextensibleADDM beam(length, EI, mu, nodes, boundary_conditions,
                                 r_penalty);
  beam.apply_initial_condition();

  EulerBeam::EulerBeamTimeDependentBCs time_dependent_bcs = {
      .end = {EulerBeam::left, EulerBeam::right},
      .type = {EulerBeam::simple_bc, EulerBeam::simple_bc},
      .history =
          {
              {.position_history = {{0.0, 0.0, 0.0},
                                    {0.0, 0.0, 0.0},
                                    {0.0, 0.0, 0.0}},
               .velocity_history = {{0.0, 0.0, 0.0},
                                    {0.0, 0.0, 0.0},
                                    {0.0, 0.0, 0.0}},
               .acceleration_history = {{0.0, 0.0, 0.0},
                                        {0.0, 0.0, 0.0},
                                        {0.0, 0.0, 0.0}}},
              {.position_history = {{length, 0.0, 0.0},
                                    {length + 0.1, 0.2, 0.0},
                                    {length + 0.2, 0.3, 0.0}},
               .velocity_history = {{0.0, 0.0, 0.0},
                                    {0.4, 0.5, 0.0},
                                    {0.6, 0.7, 0.0}},
               .acceleration_history = {{0.0, 0.0, 0.0},
                                        {0.8, 0.9, 0.0},
                                        {1.0, 1.1, 0.0}}},
          },
      .time_zero_idx = 0};

  beam.set_time_dependent_boundary_conditions(time_dependent_bcs);

  beam.solve(dt, std::array<real_t, 3>{0.0, 0.0, 0.0});

  auto &mesh_after_step = beam.get_mesh();
  const auto &centerline_after_step = mesh_after_step.get_centerline();
  const auto &velocity_after_step = mesh_after_step.get_centerline_velocity();
  EXPECT_NEAR(centerline_after_step.front()[0], 0.0, 1e-12);
  EXPECT_NEAR(centerline_after_step.front()[1], 0.0, 1e-12);
  EXPECT_NEAR(centerline_after_step.back()[0], length + 0.1, 1e-12);
  EXPECT_NEAR(centerline_after_step.back()[1], 0.2, 1e-12);
  EXPECT_NEAR(velocity_after_step.front()[0], 0.0, 1e-12);
  EXPECT_NEAR(velocity_after_step.front()[1], 0.0, 1e-12);
  EXPECT_NEAR(velocity_after_step.back()[0], 0.4, 1e-12);
  EXPECT_NEAR(velocity_after_step.back()[1], 0.5, 1e-12);

  beam.solve(dt, std::array<real_t, 3>{0.0, 0.0, 0.0});

  const auto &centerline_after_second_step = beam.get_mesh().get_centerline();
  EXPECT_NEAR(centerline_after_second_step.back()[0], length + 0.2, 1e-12);
  EXPECT_NEAR(centerline_after_second_step.back()[1], 0.3, 1e-12);
}

} // namespace ELFF
