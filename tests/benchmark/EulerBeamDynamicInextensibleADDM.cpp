

#include <gtest/gtest.h>

#include <format>
#include <iomanip>
#include <string>

#include "io/CXX/vtkHDFPolyData.hpp"
#include "models/beam/EulerBeamDynamicInextensibleADDM.hpp"

namespace ELFF {

using namespace IO::CXX;
using namespace Models;

TEST(EulerBeamDynamicInextensibleADDMTest, GlowinskiStatic)
{

  //   GTEST_LOG_(INFO) << "CTEST_FULL_OUTPUT";
  real_t length = 32.6, EI = 700., mu = 7.67, r_pentalty = 1e5;
  std::array<real_t, 3> load = { 0, -9.81 * mu, 0 };

  real_t dt = 1e-2;
  real_t t = 0;
  real_t tf = 10;
  size_t Nt = size_t(ceil(tf / dt));

  real_t dt_save = 0.1;
  size_t Nt_save = size_t(ceil(dt_save / dt));

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

  EulerBeamStaticInextensibleADDM static_beam(
    length, EI, nodes, boundary_conditions, r_pentalty);
  static_beam.apply_initial_condition();
  static_beam.solve(load);

  vtkPolyData pds = static_beam.get_mesh().to_vtk_polydata();
  vtkHDFPolyData hdf_pds("glowinski_addm_static.vtkhdf", pds);
  hdf_pds.write_new_transient(true, 0);
};

TEST(EulerBeamDynamicInextensibleADDMTest, GlowinskiDynamic)
{

  //   GTEST_LOG_(INFO) << "CTEST_FULL_OUTPUT";
  real_t length = 32.6, EI = 700., mu = 7.67, r_pentalty = 1e5;
  std::array<real_t, 3> load = { 0, -9.81 * mu, 0 };

  real_t dt = 1e-2;
  real_t t = 0;
  real_t tf = 10;
  size_t Nt = size_t(ceil(tf / dt));

  real_t dt_save = 0.1;
  size_t Nt_save = size_t(ceil(dt_save / dt));

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

  EulerBeamStaticInextensibleADDM static_beam(
    length, EI, nodes, boundary_conditions, r_pentalty);
  static_beam.apply_initial_condition();
  static_beam.solve(load);

  boundary_conditions.type[1] = EulerBeam::free_bc;

  EulerBeamDynamicInextensibleADDM dynamic_beam(
    length, EI, mu, nodes, boundary_conditions, r_pentalty);
  dynamic_beam.apply_initial_condition(static_beam.get_mesh());

  for (size_t ti = 0; ti < 4; ti++) {

    std::string filename = "glowinski_addm.vtkhdf";

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
}
