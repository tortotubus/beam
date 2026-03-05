#include "library/ibm/IBConfig.h"
#include "library/ibm/IBNode.h"
#include "library/ibm/IBMacros.h"

macro peskin_cosine_kernel_gather_dimensionless (IBNode* node) {
  coord lagpos = node->lagpos;
  coord_periodic_boundary(lagpos);

#if TREE
  foreach_cache (node->stencil) {
#else 
  foreach_neighborhood_coord(lagpos, PESKIN_SUPPORT_RADIUS) {
#endif 
    double weight = 1.0;
    coord kernel_dist = {0};
    coord cell_centre = {
      .x = x, 
      .y = y, 
      .z = z
    };
    foreach_dimension () {
      kernel_dist.x = fabs (lagpos.x-cell_centre.x) / Delta;
      if (kernel_dist.x <= PESKIN_SUPPORT_RADIUS) {
        weight *= .25 * (1 + cos (.5 * pi * kernel_dist.x));
      } else {
        weight = 0;
      }
    }
    // clang-format off
    {...}
    // clang-format on
  }
}

macro peskin_cosine_kernel_spread_dimensionless (IBNode* node) {
  coord lagpos = node->lagpos;
  coord_periodic_boundary(lagpos);
#if TREE
  foreach_cache (node->stencil) {
#else 
  foreach_neighborhood_coord_nonlocal(lagpos, PESKIN_SUPPORT_RADIUS) {
#endif 
    // Check if cell is local, since we should not write into boundary cells 
    // as they do not MPI transfer and will be over-written by the actual 
    // owners
    if (is_local(cell)) {
      double weight = 1.0;
      coord kernel_dist = {0};
      coord cell_centre = {
        .x = x, 
        .y = y, 
        .z = z
      };
      foreach_dimension () {
        kernel_dist.x = fabs (lagpos.x-cell_centre.x) / Delta;
        if (kernel_dist.x <= PESKIN_SUPPORT_RADIUS) {
          weight *= .25 * (1 + cos (.5 * pi * kernel_dist.x));
        } else {
          weight = 0;
        }
      }
      // clang-format off
          {...}
      // clang-format on
    }
  }
}

// macro peskin_cosine_kernel  (IBNode* node) {
//   coord lagpos = node->lagpos;
//   coord_periodic_boundary(lagpos);

// #if TREE
//   foreach_cache (node->stencil) {
// #else 
//   foreach_neighborhood_coord_nonlocal(lagpos, PESKIN_SUPPORT_RADIUS) {
// #endif 
//     double weight = 1.0;
//     coord kernel_dist = {0};
//     coord cell_centre = {
//       .x = x, 
//       .y = y, 
//       .z = z
//     };
//     foreach_dimension () {
//       kernel_dist.x = fabs (lagpos.x-cell_centre.x) / Delta;
//       if (kernel_dist.x <= PESKIN_SUPPORT_RADIUS) {
//         weight *= .25 * (1 + cos (.5 * pi * kernel_dist.x));
//       } else {
//         weight = 0;
//       }
//     }
//     // clang-format off
//         {...}
//     // clang-format on
//   }
// }