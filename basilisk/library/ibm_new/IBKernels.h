#include "ibm/IBConfig.h"
#include "ibm/IBNode.h"
#include "ibm/IBMacros.h"

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

macro peskin_cosine_kernel_gather_dimensionless (IBNode* node = node) {
  // bool ib_set_dirty = true;
  IBNODE_VARIABLES();
#if TREE
  foreach_neighbor_coord_level (PESKIN_SUPPORT_RADIUS, node->depth, pos) {
    if (!is_local (cell))
      continue;
#else
  foreach_neighbor_coord (PESKIN_SUPPORT_RADIUS, pos) {
#endif
    double weight = 1.0;
    coord kernel_dist = {0};
    coord cell_centre = {.x = x, .y = y, .z = z};
    foreach_dimension () {
      kernel_dist.x = fabs (d.x - cell_centre.x) / Delta;
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

macro roma_three_point_kernel_gather_dimensionless (IBNode* node = node) {
  IBNODE_VARIABLES();
#if TREE
  foreach_neighbor_coord_level (2, node->depth, pos) {
    if (!is_local (cell))
      continue;
#else
  foreach_neighbor_coord (2, pos) {
#endif
    double weight = 1.0;
    coord kernel_dist = {0};
    coord cell_centre = {.x = x, .y = y, .z = z};
    foreach_dimension () {
      kernel_dist.x = fabs (d.x - cell_centre.x) / Delta;
      weight *= roma_three_point_kernel_1d (kernel_dist.x);
    }
    // clang-format off
    {...}
    // clang-format on
  }
}

macro roma_three_point_kernel_spread_dimensionless (IBNode* node = node) {
  IBNODE_VARIABLES();
#if TREE
  foreach_neighbor_coord_level (2, node->depth, pos) {
    double weight = 1.0;
    coord kernel_dist = {0};
    coord cell_centre = {.x = x, .y = y, .z = z};
    foreach_dimension () {
      kernel_dist.x = fabs (d.x - cell_centre.x) / Delta;
      weight *= roma_three_point_kernel_1d (kernel_dist.x);
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
    foreach_dimension () {
      kernel_dist.x = fabs (d.x - cell_centre.x) / Delta;
      weight *= roma_three_point_kernel_1d (kernel_dist.x);
    }

    {
      // clang-format off
          {...}
      // clang-format on
    }
  }
#endif
}

macro roma_three_point_kernel_gather (IBNode* node = node) {
  IBNODE_VARIABLES();
#if TREE
  foreach_neighbor_coord_level (2, node->depth, pos) {
    if (!is_local (cell))
      continue;
#else
  foreach_neighbor_coord (2, pos) {
#endif
    double weight = 1.0;
    coord kernel_dist = {0};
    coord cell_centre = {.x = x, .y = y, .z = z};
    foreach_dimension () {
      kernel_dist.x = fabs (d.x - cell_centre.x);
      weight *= roma_three_point_kernel_1d (kernel_dist.x / Delta);
    }
    // clang-format off
    {...}
    // clang-format on
  }
}

macro roma_three_point_kernel_spread (IBNode* node = node) {
  IBNODE_VARIABLES();
#if TREE
  foreach_neighbor_coord_level (2, node->depth, pos) {
    double weight = 1.0;
    coord kernel_dist = {0};
    coord cell_centre = {.x = x, .y = y, .z = z};
    foreach_dimension () {
      kernel_dist.x = fabs (d.x - cell_centre.x);
      weight *= roma_three_point_kernel_1d (kernel_dist.x / Delta);
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
    foreach_dimension () {
      kernel_dist.x = fabs (d.x - cell_centre.x);
      weight *= roma_three_point_kernel_1d (kernel_dist.x / Delta);
    }

    {
      // clang-format off
          {...}
      // clang-format on
    }
  }
#endif
}

macro peskin_cosine_kernel_spread_dimensionless (IBNode* node = node) {
  // bool ib_set_dirty = false;
  IBNODE_VARIABLES();
#if TREE
  foreach_neighbor_coord_level (PESKIN_SUPPORT_RADIUS, node->depth, pos) {
    double weight = 1.0;
    coord kernel_dist = {0};
    coord cell_centre = {.x = x, .y = y, .z = z};
    foreach_dimension () {
      kernel_dist.x = fabs (d.x - cell_centre.x) / Delta;
      if (kernel_dist.x <= PESKIN_SUPPORT_RADIUS) {
        weight *= .25 * (1 + cos (.5 * pi * kernel_dist.x));
      } else {
        weight = 0;
      }
    }
    if (is_local (cell)) {
      // clang-format off
          {...}
      // clang-format on
    }
  }
#else
  foreach_neighbor_coord_nonlocal (PESKIN_SUPPORT_RADIUS, pos) {
    double weight = 1.0;
    coord kernel_dist = {0};
    coord cell_centre = {.x = x, .y = y, .z = z};
    foreach_dimension () {
      kernel_dist.x = fabs (d.x - cell_centre.x) / Delta;
      if (kernel_dist.x <= PESKIN_SUPPORT_RADIUS) {
        weight *= .25 * (1 + cos (.5 * pi * kernel_dist.x));
      } else {
        weight = 0;
      }
    }

    {
      // clang-format off
          {...}
      // clang-format on
    }
  }
#endif
}

macro peskin_cosine_kernel_gather (IBNode* node = node) {
  // bool ib_set_dirty = true;
  IBNODE_VARIABLES();
#if TREE
  foreach_neighbor_coord_level (PESKIN_SUPPORT_RADIUS, node->depth, pos) {
    if (!is_local (cell))
      continue;
#else
  foreach_neighbor_coord (PESKIN_SUPPORT_RADIUS, pos) {
#endif
    double weight = 1.0;
    coord kernel_dist = {0};
    coord cell_centre = {.x = x, .y = y, .z = z};
    foreach_dimension () {
      kernel_dist.x = fabs (d.x - cell_centre.x);
      if (kernel_dist.x <= Delta * PESKIN_SUPPORT_RADIUS) {
        weight *= (1. + cos (pi * kernel_dist.x / (Delta * PESKIN_SUPPORT_RADIUS))) / (2. * PESKIN_SUPPORT_RADIUS);
      } else {
        weight = 0;
      }
    }
    // clang-format off
    {...}
    // clang-format on
  }
}

macro peskin_cosine_kernel_spread (IBNode* node = node) {
  IBNODE_VARIABLES();
#if TREE
  foreach_neighbor_coord_level (PESKIN_SUPPORT_RADIUS, node->depth, pos) {
    double weight = 1.0;
    coord kernel_dist = {0};
    coord cell_centre = {.x = x, .y = y, .z = z};
    foreach_dimension () {
      kernel_dist.x = fabs (d.x - cell_centre.x);
      if (kernel_dist.x <= Delta *  PESKIN_SUPPORT_RADIUS)
        weight *= (1. + cos (pi * kernel_dist.x / (Delta * PESKIN_SUPPORT_RADIUS))) / (2. * PESKIN_SUPPORT_RADIUS);
      else 
        weight = 0.;
    }

    // if (weight != 0.) {
      if (is_local (cell)) {
        // clang-format off
          {...}
        // clang-format on
      }
    // }
  }
#else
  foreach_neighbor_coord_nonlocal (PESKIN_SUPPORT_RADIUS, pos) {
    double weight = 1.0;
    coord kernel_dist = {0};
    coord cell_centre = {.x = x, .y = y, .z = z};
    foreach_dimension () {
      kernel_dist.x = fabs (d.x - cell_centre.x);
      if (kernel_dist.x <= Delta *  PESKIN_SUPPORT_RADIUS)
        weight *= (1. + cos (pi * kernel_dist.x / (Delta * PESKIN_SUPPORT_RADIUS))) / (2. * PESKIN_SUPPORT_RADIUS);
      else 
        weight = 0.;
    }

    // if (weight != 0.) 
    {
      // clang-format off
          {...}
      // clang-format on
    }
  }
#endif
}
