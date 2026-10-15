#pragma once

/**
 * @file NeighborCoord.h
 *
 * Shared neighbor-iteration macros around arbitrary coordinates. These are
 * generic versions of the IBM coordinate-neighbor helpers.
 */

#include "library/common/DomainLocate.h"

/**
 * @brief Iterate leaf cells in a neighbor stencil around @p c using local
 * @c locate().
 */
macro2 foreach_domain_neighbor_coord (int r, coord c) {
  int ig = 0, jg = 0, kg = 0;
  NOT_UNUSED (ig);
  NOT_UNUSED (jg);
  NOT_UNUSED (kg);
  Point point = locate (c.x, c.y, c.z);
  if (point.level >= 0) {
    foreach_neighbor (r) {
      if (allocated (0) && is_leaf (cell)) {
        // clang-format off
        {...}
        // clang-format on
      }
    }
  }
}

/**
 * @brief Iterate leaf cells in a neighbor stencil around @p c, allowing
 * non-local/ghost cells where available.
 */
macro2 foreach_domain_neighbor_coord_nonlocal (int r, coord c) {
  int ig = 0, jg = 0, kg = 0;
  NOT_UNUSED (ig);
  NOT_UNUSED (jg);
  NOT_UNUSED (kg);
  coord d = c;
  domain_wrap_coord (d);
#if TREE
  Point point = domain_locate_nonlocal (d.x, d.y, d.z);
  if (point.level >= 0) {
    foreach_neighbor (r) {
      if (allocated (0) && is_leaf (cell)) {
        // clang-format off
        {...}
        // clang-format on
      }
    }
  }
#else
  Point point = domain_locate_nonlocal (d.x, d.y, d.z);
  if (point.level >= 0) {
    foreach_neighbor (r) {
#if dimension == 1
      if (point.i >= GHOSTS && point.i <= point.n.x + GHOSTS) {
        // clang-format off
        {...}
        // clang-format on
      }
#elif dimension == 2
      if (point.i >= GHOSTS && point.i <= point.n.x + GHOSTS) {
      if (point.j >= GHOSTS && point.j <= point.n.y + GHOSTS) {
        // clang-format off
        {...}
        // clang-format on
      }}
#else
      if (point.i >= GHOSTS && point.i <= point.n.x + GHOSTS) {
      if (point.j >= GHOSTS && point.j <= point.n.y + GHOSTS) {
      if (point.k >= GHOSTS && point.k <= point.n.z + GHOSTS) {
        // clang-format off
        {...}
        // clang-format on
      }}}
#endif
    }
  }
#endif
}

/**
 * @brief Iterate cells in a neighbor stencil around @p c at a requested level.
 */
macro2 foreach_domain_neighbor_coord_level (int r, int l, coord c) {
  int ig = 0, jg = 0, kg = 0;
  NOT_UNUSED (ig);
  NOT_UNUSED (jg);
  NOT_UNUSED (kg);
  coord d = c;
  domain_wrap_coord (d);
#if TREE
  Point point = domain_locate_level (d.x, d.y, d.z, l);
  if (point.level >= 0) {
    foreach_neighbor (r) {
      if (allocated (0) && is_leaf (cell)) {
        // clang-format off
        {...}
        // clang-format on
      }
    }
  }
#else
  Point point = domain_locate_level (d.x, d.y, d.z, l);
  if (point.level >= 0) {
    foreach_neighbor (r) {
      // clang-format off
      {...}
      // clang-format on
    }
  }
#endif
}
