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
  //u << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;
  u << 0,1,1,1,0,0,0,0,0,0;
  b.set_solution(u);

  b.assemble_residual();
  auto res = b.get_residual();
  for (size_t i = 0; i < res.size(); ++i) {
    std::cout << res[i] << " ";
  }
  std::cout << std::endl << std::endl;

  b.assemble_with_jacobian();
  auto jac_ad = b.get_jacobian();
  for (size_t i = 0; i < 10; ++i) {
    for (size_t j = 0; j < 10; ++j) {
      std::cout << jac_ad(i, j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  b.assemble_jacobian_fd();
  auto jac_fd = b.get_jacobian();
  for (size_t i = 0; i < 10; ++i) {
    for (size_t j = 0; j < 10; ++j) {
      std::cout << jac_fd(i, j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  return 0;
}
