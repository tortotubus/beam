#include "library/ibm/IBMeshManager.h"
#include "library/ibm/IBKernels.h"

trace void compute_ibm_eulvel (vector u) {
  foreach_ibnode () {
    node->eulvel.x = 0;
    peskin_cosine_kernel_dimensionless (node) {
      foreach_dimension () {
        node->eulvel.x += weight * u.x[];
      }
    }
  }
}

trace void compute_ibm_forcing () {
  foreach_ibnode () {
    foreach_dimension () {
      node->force.x = (node->lagvel.x - node->eulvel.x) / dt;
    }
  }
}

// trace void spread_ibm_forcing (vector f, scalar rho) {
//   foreach_ibnode () {
//     peskin_cosine_kernel_dimensionless (node) {
//       // double delta_h = weight / sq (Delta);
//       // coord F = node->force;
//       double dV = 16;
//       foreach_dimension() {
//         f.x[] += (dt / (rho[] * dV)) * weight * node->force.x;
//       }
//     }
//   }
// }

trace void check_dimensionless_kernel () {
  foreach_ibnode () {
    double sum = 0.;
    peskin_cosine_kernel_dimensionless (node) {
      sum += weight;
    }
    if (abs (sum - 1.0) >= 1e-2) {
      printf ("Wrong sum: %f\n", sum);
    }
  }
}