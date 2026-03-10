#include "library/ibm/IBConfig.h"
#include "library/ibm/IBNode.h"
#include "library/ibm/IBMacros.h"

macro peskin_cosine_kernel_gather_dimensionless (IBNode* node = node) {
  bool ib_set_dirty = true;
  foreach_neighbor_coord (PESKIN_SUPPORT_RADIUS, node->pos) {
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

macro peskin_cosine_kernel_spread_dimensionless (IBNode* node = node) {
  bool ib_set_dirty = false;
  // foreach_neighborhood_coord_level (node->pos, PESKIN_SUPPORT_RADIUS, node->depth) {
  foreach_neighbor_coord_level(PESKIN_SUPPORT_RADIUS, node->depth, node->pos) {
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
    if (is_local(cell)) {
    // clang-format off
          {...}
    // clang-format on
    }
  }
}
