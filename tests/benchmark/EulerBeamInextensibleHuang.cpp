
#include <gtest/gtest.h>

#include <cmath>
#include <string>

#include <elff/io/CXX/vtkHDFPolyData.hpp>
#include <elff/models/beam/EulerBeamInextensibleHuang.hpp>
#include <elff/models/beam/EulerBeamInextensibleMoM.hpp>

namespace ELFF {

using namespace IO::CXX;
using namespace Models;

namespace {

EulerBeamMesh
make_huang_initial_mesh(size_t nodes, real_t length, real_t kappa)
{
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

  return ic_mesh;
}

} // namespace

TEST(EulerBeamInextensibleHuangTest, Huang)
{
  GTEST_LOG_(INFO) << "CTEST_FULL_OUTPUT";

  const real_t length = 1.0;
  const size_t nodes = 30;
  const real_t kappa = 0.1 * M_PI;

  EulerBeamMesh ic_mesh = make_huang_initial_mesh(nodes, length, kappa);

  EulerBeam::EulerBeamBCs boundary_conditions = {
    .end = { EulerBeam::left, EulerBeam::right },
    .type = { EulerBeam::free_bc, EulerBeam::simple_bc },
    .vals = { {
                .position = { 0.0, 0.0, 0.0 },
              },
              {
                .position = { 0.0, 0.0, 0.0 },
              } }
  };

  const real_t EI = 0.01;
  const real_t mu = 1.0;
  const real_t dt = 0.02;
  const real_t tf = 0.8;
  const size_t Nt = static_cast<size_t>(std::ceil(tf / dt));

  const std::array<real_t, 3> load = { 10.0, 0.0, 0.0 };

  EulerBeamInextensibleHuang beam(length, EI, mu, nodes, boundary_conditions);
  beam.set_implicit_bending(true);
  beam.apply_initial_condition(ic_mesh);

  for (size_t ti = 0; ti < Nt; ++ti) {
    const std::string filename = "huang_fd.vtkhdf";

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

TEST(EulerBeamInextensibleHuangTest, Glowinski)
{
  GTEST_LOG_(INFO) << "CTEST_FULL_OUTPUT";

  const real_t length = 32.6;
  const real_t EI = 700.0;
  const real_t mu = 7.67;
  const real_t r_penalty = 1e5;
  const std::array<real_t, 3> load = { 0.0, -9.81 * mu, 0.0 };

  const real_t dt = 1e-2;
  const real_t tf = 10.0;
  const size_t Nt = static_cast<size_t>(std::ceil(tf / dt));
  const size_t nodes = 61;

  EulerBeam::EulerBeamBCs boundary_conditions = {
    .end = { EulerBeam::left, EulerBeam::right },
    .type = { EulerBeam::simple_bc, EulerBeam::simple_bc },
    .vals = { {
                .position = { 0.0, 0.0, 0.0 },
              },
              {
                .position = { 20.0, 0.0, 0.0 },
              } }
  };

  EulerBeamInextensibleMoM static_beam(
    length, EI, nodes, boundary_conditions, r_penalty);
  static_beam.apply_initial_condition();

  ELFF_LOG("Static Solve:");
  static_beam.solve(load);

  boundary_conditions.type[1] = EulerBeam::free_bc;

  EulerBeamInextensibleHuang beam(length, EI, mu, nodes, boundary_conditions);
  beam.set_implicit_bending(true);
  beam.apply_initial_condition(static_beam.get_mesh());

  for (size_t ti = 0; ti < Nt; ++ti) {
    const std::string filename = "glowinski_huang.vtkhdf";

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


} // namespace ELFF
