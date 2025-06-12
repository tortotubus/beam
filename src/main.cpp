#include <beam/Euler/InextenstibleEulerBeam.hpp>
#include <beam/LinAlg/Utility.hpp>
#include <iostream>

using namespace beam;

int main() {
  EulerBeamBCs bcs = {
      .end = {left, right},
      .type = {clamped_bc, free_bc},
      .vals = {{
                   .position = {0., 0., 0.},
                   .slope = {1., 0., 0.},
               },
               {
                   .position = {1., 0., 0.},
                   .slope = {1., 0., 0.},
               }},
  };

  double length = 1., EI = 1., load = -10., area = 1., r = 1e4;
  size_t nnodes = 20;

  InextensibleEulerBeam b(length, EI, load, area, nnodes, r, bcs);

  b.apply_initial_condition();
  //b.assemble_residual();
  //b.apply_boundary_conditions();

  //auto res = b.get_residual();
  //std::cout << res << "\n\n";

  //b.assemble_system();
  //auto jac_ad = b.get_jacobian();
  //std::cout << jac_ad << "\n\n";

  //b.assemble_jacobian_fd();
  //auto jac_fd = b.get_jacobian();
  //std::cout << jac_fd << "\n\n";

  b.solve();

  auto centerline = b.get_centerline();
  std::cout << centerline << "\n\n";

  b.gnuplot();

  return 0;
}
