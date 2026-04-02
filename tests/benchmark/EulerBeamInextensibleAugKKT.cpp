#include <elff/io/CXX/vtkHDFPolyData.hpp>
#include <elff/models/beam/EulerBeamInextensibleAugKKT.hpp>
#include <gtest/gtest.h>
#include <string>

#include "EulerBeamStaticInextensibleReferences.hpp"

namespace ELFF {

using namespace IO::CXX;
using namespace Models;

TEST(EulerBeamInextensibleAugKKTTest, Glowinski)
{

  GTEST_LOG_(INFO) << "CTEST_FULL_OUTPUT";
  real_t length = 32.6, EI = 700., mu = 7.67, r_penalty = 11e2;
  std::array<real_t, 3> load = { 0, -9.81 * mu, 0 };

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
    .end = { EulerBeam::left, EulerBeam::right },
    .type = { EulerBeam::simple_bc, EulerBeam::simple_bc },
    .vals = { {
                .position = { 0, 0, 0 },
              },
              {
                .position = { 20, 0, 0 },
              } }
  };

  EulerBeamInextensibleAugKKT static_beam(
    length, EI, nodes, boundary_conditions, r_penalty);
  static_beam.apply_initial_condition();

  ELFF_LOG("Static Solve:");
  static_beam.solve(load);

  boundary_conditions.type[1] = EulerBeam::free_bc;

  EulerBeamInextensibleAugKKT dynamic_beam(
    length, EI, mu, nodes, boundary_conditions, r_penalty);

  dynamic_beam.apply_initial_condition(static_beam.get_mesh());

  ELFF_LOG("Dynamic Solve:");
  for (size_t ti = 0; ti < Nt; ti++) {

    std::string filename = "glowinski_augkkt.vtkhdf";

    if (ti == 0) {
      vtkPolyData pd = dynamic_beam.get_mesh().to_vtk_polydata();
      vtkHDFPolyData hdf_pd(filename, pd);
      hdf_pd.write_new_transient(true, ti * dt);
    } else {
      vtkPolyData pd = dynamic_beam.get_mesh().to_vtk_polydata();
      vtkHDFPolyData hdf_pd(filename, pd);
      hdf_pd.append_transient(ti * dt);
    }

    dynamic_beam.solve(dt, load);
  }
};

TEST(EulerBeamInextensibleAugKKTTest, Huang)
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
                .position = { 0, 0, 0 },
              },
              {
              } }
  };

  real_t EI = 0.01;
  real_t mu = 1;
  real_t r_penalty = 1e4;

  real_t dt = 0.02;
  real_t t = 0;
  real_t tf = 0.8;
  size_t Nt = size_t(ceil(tf / dt));
  (void)t;

  std::array<real_t, 3> load = { 10, 0, 0 };

  EulerBeamInextensibleAugKKT dynamic_beam(
    length, EI, mu, nodes, boundary_conditions, r_penalty);

  dynamic_beam.apply_initial_condition(ic_mesh);

  for (size_t ti = 0; ti < Nt; ti++) {

    std::string filename = "huang_augkkt.vtkhdf";

    if (ti == 0) {
      vtkPolyData pd = dynamic_beam.get_mesh().to_vtk_polydata();
      vtkHDFPolyData hdf_pd(filename, pd);
      hdf_pd.write_new_transient(true, ti * dt);
    } else {
      vtkPolyData pd = dynamic_beam.get_mesh().to_vtk_polydata();
      vtkHDFPolyData hdf_pd(filename, pd);
      hdf_pd.append_transient(ti * dt);
    }

    dynamic_beam.solve(dt, load);
  }
};

TEST(EulerBeamInextensibleAugKKTTest, BisshoppAndDrucker)
{
  real_t length = 1., EI = 1., area = 1., r_pentalty = 1e4;
  size_t nodes = 250;

  real_t tip_force_y = -1;

  double comparison_tol = 4e-9;

  EulerBeam::EulerBeamBCs boundary_conditions = {
    .end = { EulerBeam::left, EulerBeam::right },
    .type = { EulerBeam::clamped_bc, EulerBeam::point_force_bc },
    .vals = { { .position = { 0, 0, 0 }, .slope = { 1, 0, 0 } },
              { .force = { 0, tip_force_y, 0 } } }
  };

  EulerBeamInextensibleAugKKT sparse_beam(
    length, EI, nodes, boundary_conditions, r_pentalty);

  sparse_beam.solve();
  sparse_beam.get_mesh().plot_gnuplot("Bisshopp and Drucker AugKKT Sparse");

  EulerBeamMesh& mesh = sparse_beam.get_mesh();
  auto centerline = mesh.get_centerline();
  std::array<real_t, 3> tip = centerline[nodes - 1];

  BisshoppAndDrucker1945Result res =
    BisshoppAndDrucker1945(length, EI, -tip_force_y);

  EXPECT_NEAR(std::abs(length - tip[0]), res.A, comparison_tol);
  EXPECT_NEAR(std::abs(tip[1]), res.delta, comparison_tol);
}


}
