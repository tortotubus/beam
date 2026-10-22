#include "library/ibm/IBConfig.h"
#include "library/ibm/IBNode.h"
#include "library/ibm/IBMacros.h"

static inline double peskin_cosine_kernel_1d (double r) {
  double q = fabs (r);
  if (q <= 2)
    return (1. + cos (pi * q / 2)) / (2. * 2);
  return 0.;
}

static inline double roma_three_point_kernel_1d (double r) {
  double q = fabs (r);
  if (q <= 0.5)
    return (1. + sqrt (1. - 3. * q * q)) / 3.;
  if (q <= 1.5) {
    double a = 1. - 3. * sq (1. - q);
    return (5. - 3. * q - sqrt (max (0., a))) / 6.;
  }
  return 0.;
}

// TODO: Test if using function pointer is any slower than other methods of runtime
// dispatching
typedef double (*ib_kernel_1d) (double);
ib_kernel_1d default_kernel = peskin_cosine_kernel_1d;

macro ib_kernel_gather_dimensionless (IBNode* node = node) {
  IBNODE_VARIABLES ();
  foreach_neighbor_coord (2, pos) {
    double weight = 1.0;
    coord kernel_dist = {0};
    coord cell_centre = {.x = x, .y = y, .z = z};
    coord kernel_delta = {0};
    domain_displacement (kernel_delta, d, cell_centre);
    foreach_dimension () {
      kernel_dist.x = fabs (kernel_delta.x) / Delta;
      weight *= default_kernel (kernel_dist.x);
    }
    // clang-format off
    {...}
    // clang-format on
  }
}

macro ib_kernel_spread_dimensionless (IBNode* node = node) {
  IBNODE_VARIABLES ();
#if TREE
  foreach_neighbor_coord_level (2, node->depth, pos) {
    double weight = 1.0;
    coord kernel_dist = {0};
    coord cell_centre = {.x = x, .y = y, .z = z};
    coord kernel_delta = {0};
    domain_displacement (kernel_delta, d, cell_centre);
    foreach_dimension () {
      kernel_dist.x = fabs (kernel_delta.x) / Delta;
      weight *= default_kernel (kernel_dist.x);
    }
    if (is_local (cell)) {
      // clang-format off
      {...}
      // clang-format on
    }
  }
#else
  foreach_neighbor_coord_nonlocal (2, pos) {
    double weight = 1.0;
    coord kernel_dist = {0};
    coord cell_centre = {.x = x, .y = y, .z = z};
    coord kernel_delta = {0};
    domain_displacement (kernel_delta, d, cell_centre);
    foreach_dimension () {
      kernel_dist.x = fabs (kernel_delta.x) / Delta;
      weight *= default_kernel (kernel_dist.x);
    }

    {
      // clang-format off
      {...}
      // clang-format on
    }
  }
#endif
}
