
#include <gtest/gtest.h>

#include <cmath>
#include <string>

#include <elff/io/CXX/vtkHDFPolyData.hpp>
#include <elff/models/beam/EulerBeamHuang.hpp>

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

void
expect_mesh_finite(EulerBeamMesh& mesh)
{
  const auto& centerline = mesh.get_centerline();
  for (const auto& p : centerline) {
    ASSERT_TRUE(std::isfinite(p[0]));
    ASSERT_TRUE(std::isfinite(p[1]));
    ASSERT_TRUE(std::isfinite(p[2]));
  }
}

} // namespace

TEST(EulerBeamHuangTest, Huang)
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

  EulerBeamHuang beam(length, EI, mu, nodes, boundary_conditions);
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

TEST(EulerBeamHuangTest, HuangClampedSmoke)
{
  GTEST_LOG_(INFO) << "CTEST_FULL_OUTPUT";

  const real_t length = 1.0;
  const size_t nodes = 30;
  const real_t kappa = 0.1 * M_PI;

  EulerBeamMesh ic_mesh = make_huang_initial_mesh(nodes, length, kappa);

  EulerBeam::EulerBeamBCs boundary_conditions = {
    .end = { EulerBeam::left, EulerBeam::right },
    .type = { EulerBeam::free_bc, EulerBeam::clamped_bc },
    .vals = { {
                .position = { 0.0, 0.0, 0.0 },
              },
              {
                .position = { 0.0, 0.0, 0.0 },
                .slope = { -std::cos(kappa), -std::sin(kappa), 0.0 },
              } }
  };

  const real_t EI = 0.01;
  const real_t mu = 1.0;
  const real_t dt = 0.02;
  const std::array<real_t, 3> load = { 10.0, 0.0, 0.0 };

  EulerBeamHuang beam(length, EI, mu, nodes, boundary_conditions);
  beam.apply_initial_condition(ic_mesh);

  for (size_t ti = 0; ti < 5; ++ti) {
    beam.solve(dt, load);
  }

  const auto& centerline = beam.get_mesh().get_centerline();
  ASSERT_NEAR(centerline.back()[0], 0.0, 1e-12);
  ASSERT_NEAR(centerline.back()[1], 0.0, 1e-12);
}

TEST(EulerBeamHuangTest, HuangMirroredSimpleFreeSmoke)
{
  const real_t length = 1.0;
  const size_t nodes = 30;
  const real_t kappa = 0.1 * M_PI;

  EulerBeamMesh ic_mesh = make_huang_initial_mesh(nodes, length, kappa);
  const auto& ic_centerline = ic_mesh.get_centerline();

  EulerBeam::EulerBeamBCs boundary_conditions = {
    .end = { EulerBeam::left, EulerBeam::right },
    .type = { EulerBeam::simple_bc, EulerBeam::free_bc },
    .vals = { {
                .position = { ic_centerline.front()[0], ic_centerline.front()[1], 0.0 },
              },
              {
                .position = { 0.0, 0.0, 0.0 },
              } }
  };

  EulerBeamHuang beam(length, 0.01, 1.0, nodes, boundary_conditions);
  beam.apply_initial_condition(ic_mesh);
  for (size_t ti = 0; ti < 5; ++ti) {
    beam.solve(0.02, std::array<real_t, 3>{ 10.0, 0.0, 0.0 });
  }

  const auto& centerline = beam.get_mesh().get_centerline();
  ASSERT_NEAR(centerline.front()[0], ic_centerline.front()[0], 1e-12);
  ASSERT_NEAR(centerline.front()[1], ic_centerline.front()[1], 1e-12);
  expect_mesh_finite(beam.get_mesh());
}

TEST(EulerBeamHuangTest, HuangFreeFreeSmoke)
{
  const real_t length = 1.0;
  const size_t nodes = 30;
  const real_t kappa = 0.1 * M_PI;

  EulerBeamMesh ic_mesh = make_huang_initial_mesh(nodes, length, kappa);

  EulerBeam::EulerBeamBCs boundary_conditions = {
    .end = { EulerBeam::left, EulerBeam::right },
    .type = { EulerBeam::free_bc, EulerBeam::free_bc },
    .vals = { {
                .position = { 0.0, 0.0, 0.0 },
              },
              {
                .position = { 0.0, 0.0, 0.0 },
              } }
  };

  EulerBeamHuang beam(length, 0.01, 1.0, nodes, boundary_conditions);
  beam.apply_initial_condition(ic_mesh);
  for (size_t ti = 0; ti < 5; ++ti) {
    beam.solve(0.02, std::array<real_t, 3>{ 10.0, 0.0, 0.0 });
  }

  expect_mesh_finite(beam.get_mesh());
}

TEST(EulerBeamHuangTest, HuangSimpleSimpleSmoke)
{
  const real_t length = 1.0;
  const size_t nodes = 30;
  const real_t kappa = 0.1 * M_PI;

  EulerBeamMesh ic_mesh = make_huang_initial_mesh(nodes, length, kappa);
  const auto& ic_centerline = ic_mesh.get_centerline();

  EulerBeam::EulerBeamBCs boundary_conditions = {
    .end = { EulerBeam::left, EulerBeam::right },
    .type = { EulerBeam::simple_bc, EulerBeam::simple_bc },
    .vals = { {
                .position = { ic_centerline.front()[0], ic_centerline.front()[1], 0.0 },
              },
              {
                .position = { ic_centerline.back()[0], ic_centerline.back()[1], 0.0 },
              } }
  };

  EulerBeamHuang beam(length, 0.01, 1.0, nodes, boundary_conditions);
  beam.apply_initial_condition(ic_mesh);
  for (size_t ti = 0; ti < 5; ++ti) {
    beam.solve(0.02, std::array<real_t, 3>{ 10.0, 0.0, 0.0 });
  }

  const auto& centerline = beam.get_mesh().get_centerline();
  ASSERT_NEAR(centerline.front()[0], ic_centerline.front()[0], 1e-12);
  ASSERT_NEAR(centerline.front()[1], ic_centerline.front()[1], 1e-12);
  ASSERT_NEAR(centerline.back()[0], ic_centerline.back()[0], 1e-12);
  ASSERT_NEAR(centerline.back()[1], ic_centerline.back()[1], 1e-12);
  expect_mesh_finite(beam.get_mesh());
}

TEST(EulerBeamHuangTest, Huang3DSmoke)
{
  const real_t length = 1.0;
  const size_t nodes = 30;
  const real_t kappa = 0.1 * M_PI;

  EulerBeamMesh ic_mesh = make_huang_initial_mesh(nodes, length, kappa);
  auto& ic_centerline = ic_mesh.get_centerline();
  auto& ic_slope = ic_mesh.get_slope();
  auto& ic_s = ic_mesh.get_curvilinear_axis();

  for (size_t ni = 0; ni < nodes; ++ni) {
    const real_t s = ic_s[ni];
    ic_centerline[ni][2] = 0.05 * std::sin(M_PI * s / length);
  }

  for (size_t ni = 0; ni < nodes; ++ni) {
    if (ni == 0) {
      ic_slope[ni][2] = (ic_centerline[1][2] - ic_centerline[0][2]) / ic_mesh.get_ds();
    } else if (ni + 1 == nodes) {
      ic_slope[ni][2] =
        (ic_centerline[ni][2] - ic_centerline[ni - 1][2]) / ic_mesh.get_ds();
    } else {
      ic_slope[ni][2] =
        (ic_centerline[ni + 1][2] - ic_centerline[ni - 1][2]) / (2.0 * ic_mesh.get_ds());
    }
  }

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

  EulerBeamHuang beam(length, 0.01, 1.0, nodes, boundary_conditions);
  beam.apply_initial_condition(ic_mesh);
  for (size_t ti = 0; ti < 5; ++ti) {
    beam.solve(0.02, std::array<real_t, 3>{ 10.0, 0.0, 1.0 });
  }

  const auto& centerline = beam.get_mesh().get_centerline();
  ASSERT_TRUE(std::isfinite(centerline.front()[2]));
  ASSERT_TRUE(std::isfinite(centerline.back()[2]));
  ASSERT_NEAR(centerline.back()[2], 0.0, 1e-12);
  expect_mesh_finite(beam.get_mesh());
}

} // namespace ELFF
