#include <beam/Euler/InextenstibleEulerBeam.hpp>
#include <iostream>

using namespace beam;

int main() {
  EulerBeamBCs bcs = {
      .end = {left, right},
      .type = {clamped_bc, free_bc},
      .vals = {{
                   .position = {0., 0., 0.},
                   .angle = {0., 0., 0.},
               },
               {
                   .position = {1., 0., 0.},
                   .angle = {0., 0., 0.},
               }},
  };

  double length = 1., EI = 1., load = -1., area = 1., r = 1e4;
  size_t nnodes = 2;

  InextensibleEulerBeam b(length, EI, load, area, nnodes, r, bcs);

  Eigen::VectorXd u(10);
  u << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;

  b.set_solution(u);
  b.assemble_residual();
  Eigen::VectorXd res = b.get_residual();
  for (size_t i = 0; i < res.size(); ++i) {
    std::cout << res[i] << " ";
  }
  std::cout << std::endl;

  return 0;
}
