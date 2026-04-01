

#include <elff/io/CXX/vtkHDFPolyData.hpp>
#include <elff/models/beam/EulerBeamInextensibleGGL.hpp>
#include <elff/models/beam/EulerBeamInextensibleMoM.hpp>
#include <gtest/gtest.h>

#include <cmath>
#include <string>

namespace ELFF {

using namespace IO::CXX;
using namespace Models;

TEST(EulerBeamInextensibleGGLTest, Glowinski)
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

  EulerBeamInextensibleMoM static_beam(
    length, EI, nodes, boundary_conditions, r_penalty);
  static_beam.apply_initial_condition();

  ELFF_LOG("Static Solve:");
  static_beam.solve(load);

  boundary_conditions.type[1] = EulerBeam::free_bc;

  EulerBeamInextensibleGGL dynamic_beam(
    length, EI, mu, nodes, boundary_conditions, r_penalty);
  dynamic_beam.apply_initial_condition(static_beam.get_mesh());

  ELFF_LOG("Dynamic Solve:");
  for (size_t ti = 0; ti < Nt; ++ti) {
    const std::string filename = "glowinski_ggl.vtkhdf";
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

TEST(EulerBeamInextensibleGGLTest, Huang)
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

  EulerBeamInextensibleGGL dynamic_beam(
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
