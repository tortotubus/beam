

#include <gtest/gtest.h>

#include <format>
#include <iomanip>
#include <string>

#include "io/CXX/vtkHDFPolyData.hpp"
#include "models/beam/EulerBeamDynamicInextensibleMoM.hpp"

namespace beam {

using namespace io::CXX;

TEST(EulerBeamDynamicInextensibleMoMTest, Glowinski)
{

  GTEST_LOG_(INFO) << "CTEST_FULL_OUTPUT";
  real_t length = 32.6, EI = 700., mu = 7.67, r_pentalty = 1e5;
  std::array<real_t, 3> load = { 0, -9.81 * mu, 0 };

  real_t dt = 1e-2;
  real_t t = 0;
  real_t tf = 10;
  size_t Nt = size_t(ceil(tf / dt));

  real_t dt_save = 0.1;
  size_t Nt_save = size_t(ceil(dt_save / dt));

  size_t nodes = 61;

  EulerBeamBCs boundary_conditions = { .end = { left, right },
                                       .type = { simple_bc, simple_bc },
                                       .vals = { {
                                                   .position = { 0, 0, 0 },
                                                 },
                                                 {
                                                   .position = { 20, 0, 0 },
                                                 } } };

  EulerBeamStaticInextensibleMoM static_beam(
    length, EI, nodes, boundary_conditions, r_pentalty);
  static_beam.apply_initial_condition();
  static_beam.solve(load);
  static_beam.plot("Static Inextensible Euler Beam");

  boundary_conditions.type[1] = free_bc;

  EulerBeamDynamicInextensibleMoM dynamic_beam(
    length, EI, mu, nodes, boundary_conditions, r_pentalty);
  dynamic_beam.apply_initial_condition(static_beam.get_mesh());

  for (size_t ti = 0; ti < Nt; ti++) {
    dynamic_beam.solve(dt, load);

    if (ti % Nt_save == 0) {
      vtkPolyData pd = dynamic_beam.get_mesh().to_vtk_polydata();

      std::ostringstream fname_oss;
      fname_oss  << "glowinski_mom_" << ti << ".hdf";
      std::string fname = fname_oss.str();

      vtkHDFPolyData hdf_pd(fname, pd);
      hdf_pd.write_new_static(true);
    }
  }
};

TEST(EulerBeamDynamicInextensibleMoMTest, Huang)
{
  GTEST_LOG_(INFO) << "CTEST_FULL_OUTPUT";

  real_t length = 1;
  size_t nodes = 30;
  real_t kappa = 0.1 * M_PI;

  real_t x0 = 0, y0 = 0;

  EulerBeamMesh ic_mesh(nodes, length);

  std::vector<std::array<real_t, 3>>& ic_centerline = ic_mesh.get_centerline();
  std::vector<std::array<real_t, 3>>& ic_slope = ic_mesh.get_slope();
  std::vector<std::array<real_t, 3>>& ic_velocity =
    ic_mesh.get_centerline_velocity();
  std::vector<real_t>& ic_s = ic_mesh.get_curvilinear_axis();

  for (size_t ni = 0; ni < nodes; ++ni) {
    real_t s = ic_s[ni];

    ic_centerline[ni][0] = x0 + (length - s) * std::cos(kappa);
    ic_centerline[ni][1] = y0 + (length - s) * std::sin(kappa);
    ic_centerline[ni][2] = 0.;

    // Tangent dX/ds:
    ic_slope[ni][0] = -std::cos(kappa);
    ic_slope[ni][1] = -std::sin(kappa);
    ic_slope[ni][2] = 0.;

    ic_velocity[ni][0] = 0.;
    ic_velocity[ni][1] = 0.;
    ic_velocity[ni][2] = 0.;
  }

  ic_mesh.plot_gnuplot("Initial condition");

  EulerBeamBCs boundary_conditions = { .end = { left, right },
                                       .type = { free_bc, simple_bc },
                                       .vals = { {
                                                   .position = { 0, 0, 0 },
                                                 },
                                                 {
                                                   //
                                                 } } };

  real_t EI = 0.01;
  real_t mu = 1;
  real_t r_penalty = 1e5;

  real_t dt = 0.02;
  real_t t = 0;
  real_t tf = 0.8;
  size_t Nt = size_t(ceil(tf / dt));

  std::array<real_t, 3> load = { 10, 0, 0 };

  EulerBeamDynamicInextensibleMoM dynamic_beam(
    length, EI, mu, nodes, boundary_conditions, r_penalty);
  dynamic_beam.apply_initial_condition(ic_mesh);

  for (size_t ti = 0; ti < Nt; ti++) {
    vtkPolyData pd = dynamic_beam.get_mesh().to_vtk_polydata();

    std::ostringstream fname_oss;
    fname_oss  << "huang_mom_" << ti << ".hdf";
    std::string fname = fname_oss.str();

    vtkHDFPolyData hdf_pd(fname, pd);
    hdf_pd.write_new_static(true);

    dynamic_beam.solve(dt, load);
  }
};

}
