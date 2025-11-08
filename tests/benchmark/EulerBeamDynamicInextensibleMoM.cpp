

#include <gtest/gtest.h>

#include <format>
#include <string>

#include "models/beam/EulerBeamDynamicInextensibleMoM.hpp"
#include "io/VTKHDFPolyData.hpp"

using beam::EulerBeamDynamicInextensibleMoM;
using beam::EulerBeam;
using beam::EulerBeamMesh;
using beam::EulerBeamBCs;
using beam::real_t;
using beam::left;
using beam::right;
using beam::simple_bc;
using beam::clamped_bc;
using beam::free_bc;
using beam::point_force_bc;
using beam::point_torque_bc;
using beam::VTKPolyData;
using beam::VTKHDFPolyData;

TEST(EulerBeamDynamicInextensibleMoMTest, StaticSolveAndPlot) {

  real_t length = 1., EI = 1., mu = 1., r_pentalty = 1e3;
  std::array<real_t, 3> load = {0,-10,-1};
  real_t dt = 0.05;
  real_t tf = 2.0;
  size_t nodes = 15;

  EulerBeamBCs boundary_conditions = {
    .end = {left, right},
    .type = {clamped_bc, free_bc},
    .vals = {{
      .position = {0,0,0},
      .slope = {1,0,0}
    }, {}}
  };

  EulerBeamDynamicInextensibleMoM static_beam(length, EI, mu, nodes, boundary_conditions, r_pentalty);
  static_beam.apply_initial_condition();

  GTEST_LOG_(INFO) << "CTEST_FULL_OUTPUT";

  

  for (size_t ti = 0; ti < 1000; ti++) {
    static_beam.solve(dt, load);
    
    if (ti % 2 == 0) {
      VTKPolyData pd = static_beam.get_mesh().to_vtk_polydata();

      std::ostringstream fname_oss;
      fname_oss << "eb_dynamic_mom_" << ti << ".hdf"; 
      std::string fname = fname_oss.str();  

      VTKHDFPolyData hdf_pd(fname, pd);
      hdf_pd.write();
    }
  }

  static_beam.plot("Static Inextensible Euler Beam");
};

