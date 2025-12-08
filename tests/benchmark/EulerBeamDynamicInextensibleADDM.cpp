

#include "models/beam/EulerBeamDynamicInextensibleADDM.hpp"
#include <gtest/gtest.h>

// TEST(EulerBeamDynamicInextensibleADDMTest, StaticSolveAndPlot)
// {

//   real_t length = 1., EI = 1., mu = 1., r_pentalty = 1e3;
//   std::array<real_t, 3> load = { 0., -1., 0. };
//   real_t dt = 0.0005;
//   real_t tf = 1.0;
//   size_t nodes = 15;

//   EulerBeamBCs boundary_conditions = {
//     .end = { left, right },
//     .type = { simple_bc, free_bc },
//     .vals = { { .position = { 0, 0, 0 }, .slope = { 1, 0, 0 } },
//               { .position = { length, 0, 0 }, .slope = { 0, 0, 0 } } }
//   };

//   EulerBeamDynamicInextensibleADDM static_beam(length, EI, mu, nodes, boundary_conditions, r_pentalty);

//   for (size_t ti = 0; ti < 3; ti++)
//     static_beam.solve(dt, load);

//   static_beam.plot("Dynamic Inextensible Euler Beam Final");
// };
