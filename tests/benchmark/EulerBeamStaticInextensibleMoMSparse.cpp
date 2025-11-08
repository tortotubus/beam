

#include <gtest/gtest.h>
#include "models/beam/EulerBeamStaticInextensibleMoMSparse.hpp"

using beam::EulerBeamStaticInextensibleMoMSparse;
using beam::EulerBeam;
using beam::EulerBeamBCs;
using beam::real_t;
using beam::left;
using beam::right;
using beam::clamped_bc;
using beam::free_bc;


// TEST(EulerBeamStaticInextensibleMoMSparseTest, SolveUniformLoadAndPlot) {

//   real_t length = 1., EI = 1., mu = 1., r_pentalty = 1e4;
//   std::array<real_t,3> load = {0.,-10.,0.};
//   size_t nodes = 100;

//   EulerBeamBCs boundary_conditions = {
//     .end = {left, right},
//     .type = {clamped_bc, free_bc},
//     .vals = {{
//       .position = {0,0,0},
//       .slope = {1,0,0}
//     }, {
//       .position = {1,0,0},
//       .slope = {1,0,0}
//     }}
//   };

//   EulerBeamStaticInextensibleMoMSparse static_beam(length, EI, nodes, boundary_conditions, r_pentalty);

//   GTEST_LOG_(INFO) << "CTEST_FULL_OUTPUT";
//   static_beam.solve(load);
//   static_beam.plot("Static Inextensible Euler Beam");
// };