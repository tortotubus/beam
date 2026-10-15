#include "library/common/DomainGeometry.h"
#include "library/ibm/IBLocate.h"

/**
 * @def coord_periodic_boundary
 *
 * @brief Wraps the @c coord struct when periodic boundary conditions are set
 */
macro coord_periodic_boundary (coord c) {
  domain_wrap_coord (c);
}

/**
 * @def periodic_minimum_image_delta
 *
 * @brief Maps a displacement to the nearest periodic image.
 */
macro periodic_minimum_image_delta (double delta, bool is_periodic) {
  if (is_periodic) {
    delta = fmod (delta + 0.5 * L0, L0);
    if (delta < 0.)
      delta += L0;
    delta -= 0.5 * L0;
  }
}

/**
 * @def point_periodic_boundary
 *
 * @brief Wraps the @c Point struct when periodic boundary conditions are set
 */
macro point_periodic_boundary (Point p) {
  domain_wrap_point (p);
}

/**
 * @def foreach_neighbor_coord
 */
macro2 foreach_neighbor_coord(int r, coord c) {
  int ig = 0, jg = 0, kg = 0;
  NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
  coord d = c;
  domain_wrap_coord (d);
  Point point = locate(d.x,d.y,d.z);
  if (point.level >= 0) {
    foreach_neighbor(r) {
      if (allocated(0) && is_leaf(cell)) {
        // clang-format off
        {...}
        // clang-format on
      }
    }
  }
}
/**
 * @def foreach_neighbor_coord_nonlocal
 */
macro2 foreach_neighbor_coord_nonlocal(int r, coord c) {
  int ig = 0, jg = 0, kg = 0;
  NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
  coord d = c;
  domain_wrap_coord (d);
#if TREE
  Point point = domain_locate_nonlocal(d.x,d.y,d.z);
  if (point.level >= 0) {
    foreach_neighbor(r) {
      if (allocated (0) && is_leaf (cell)) {
        // clang-format off
        {...}
        // clang-format on
      }
    }
  }
#else // !TREE
  Point point = domain_locate_nonlocal(d.x,d.y,d.z);
  if (point.level >= 0) {
    foreach_neighbor(r) {
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
#else // dimension == 3
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
#endif // !TREE
}

/**
 * @def foreach_neighbor_coord_level
 */
macro2 foreach_neighbor_coord_level(int r, int l, coord c) {
  int ig = 0, jg = 0, kg = 0;
  NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
  coord d = c;
  domain_wrap_coord (d);
#if TREE
  int effective_l = l;
  if (effective_l > depth())
    effective_l = depth();
  Point point = domain_locate_level(d.x,d.y,d.z,effective_l);
  if (point.level >= 0) {
    foreach_neighbor(r) {
      if (allocated(0) && is_leaf(cell))
      // clang-format off
      {...}
      // clang-format on
    }
  }
#else 
  Point point = domain_locate_level(d.x,d.y,d.z,l);
  if (point.level >= 0) {
    foreach_neighbor(r) {
      // clang-format off
      {...}
      // clang-format on
    }
  }
#endif 

}
