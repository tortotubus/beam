

#include <gtest/gtest.h>

#include <format>
#include <string>
#include <iomanip>

#include "models/beam/EulerBeamDynamicInextensibleAugKKT.hpp"
#include "models/beam/EulerBeamStaticInextensibleAugKKT.hpp"
#include "io/CXX/vtkHDFPolyData.hpp"

namespace ELFF {

using namespace io::CXX;
// using ELFF::vtkHDFPolyData;

// TEST(EulerBeamDynamicInextensibleAugKKTTest, StaticSolveAndPlot) {

//   real_t length = 1., EI = 1., mu = 1., r_pentalty = 1e5;
//   std::array<real_t, 3> load = {0,-10,0};

//   real_t dt = 0.1;
//   real_t t = 0;
//   real_t tf = 100;
//   size_t Nt = size_t(ceil(tf/dt));

//   real_t dt_save = 0.5;
//   size_t Nt_save = size_t(ceil(dt_save/dt));

//   size_t nodes = 25;

//   EulerBeamBCs boundary_conditions = {
//     .end = {left, right},
//     .type = {clamped_bc, free_bc},
//     .vals = {{
//       .position = {0,0,0},
//       .slope = {1,0,0}
//     }, {}}
//   };

//   EulerBeamDynamicInextensibleAugKKT static_beam(length, EI, mu, nodes, boundary_conditions, r_pentalty);
//   static_beam.apply_initial_condition();

//   GTEST_LOG_(INFO) << "CTEST_FULL_OUTPUT";

  

//   for (size_t ti = 0; ti < Nt; ti++) {
//     static_beam.solve(dt, load);
    
//     if (ti % Nt_save == 0) {
//       vtkPolyData pd = static_beam.get_mesh().to_vtk_polydata();

//       std::ostringstream fname_oss;
//       fname_oss  << "eb_dynamic_akkt_" << ti << ".hdf"; 
//       std::string fname = fname_oss.str();  

//       vtkHDFPolyData hdf_pd(fname, pd);
//       hdf_pd.write_new_static(true);
//     }
//   }

//   static_beam.plot("Static Inextensible Euler Beam");
// };

// TEST(EulerBeamDynamicInextensibleAugKKTTest, Glowinski) {

//   real_t length = 32.6, EI = 700., mu = 7.67, r_pentalty = 1e5;
//   std::array<real_t, 3> load = {0,-9.81*mu,0};

//   real_t dt = 1e-2;
//   real_t t = 0;
//   real_t tf = 10;
//   size_t Nt = size_t(ceil(tf/dt));

//   real_t dt_save = 0.1;
//   size_t Nt_save = size_t(ceil(dt_save/dt));

//   size_t nodes = 61;

//   EulerBeamBCs boundary_conditions = {
//     .end = {left, right},
//     .type = {simple_bc, simple_bc},
//     .vals = {{
//       .position = {0,0,0},
//     }, {
//       .position = {20,0,0},
//     }}
//   };

//   EulerBeamStaticInextensibleAugKKT static_beam(length, EI, nodes, boundary_conditions, r_pentalty);
//   static_beam.apply_initial_condition();
//   static_beam.solve(load);
//   static_beam.plot("Static Inextensible Euler Beam");

//   boundary_conditions.type[1] = free_bc;

//   EulerBeamDynamicInextensibleAugKKT dynamic_beam(length, EI, mu, nodes, boundary_conditions, r_pentalty);
//   dynamic_beam.apply_initial_condition(static_beam.get_mesh());

//   for (size_t ti = 0; ti < Nt; ti++) {
//     dynamic_beam.solve(dt, load);
    
//     if (ti % Nt_save == 0) {
//       vtkPolyData pd = dynamic_beam.get_mesh().to_vtk_polydata();

//       std::ostringstream fname_oss;
//       fname_oss  << "glowinski_akkt_" << ti << ".hdf"; 
//       std::string fname = fname_oss.str();  

//       vtkHDFPolyData hdf_pd(fname, pd);
//       hdf_pd.write_new_static(true);
//     }
//   }

// };

}