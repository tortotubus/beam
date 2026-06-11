#include "ibm/IBLocate.h"

/**
 * @def coord_periodic_boundary
 *
 * @brief Wraps the @c coord struct when periodic boundary conditions are set
 */
macro coord_periodic_boundary (coord c) {
  if (Period.x) {
    c.x = X0 + fmod (fmod (c.x - X0, L0) + L0, L0);
  }
  if (Period.y) {
    c.y = Y0 + fmod (fmod (c.y - Y0, L0) + L0, L0);
  }
  if (Period.z) {
    c.z = Z0 + fmod (fmod (c.z - Z0, L0) + L0, L0);
  }
}

/**
 * @def point_periodic_boundary
 *
 * @brief Wraps the @c Point struct when periodic boundary conditions are set
 */
macro point_periodic_boundary (Point p) {
  int N = (1 << p.level);
#if dimension >= 1
  if (Period.x) {
    int t = p.i - GHOSTS;
    t = (t % N + N) % N;
    p.i = t + GHOSTS;
  }
#endif
#if dimension >= 2
  if (Period.y) {
    int t = p.j - GHOSTS;
    t = (t % N + N) % N;
    p.j = t + GHOSTS;
  }
#endif
#if dimension >= 3
  if (Period.z) {
    int t = p.k - GHOSTS;
    t = (t % N + N) % N;
    p.k = t + GHOSTS;
  }
#endif
}

/**
 * @def foreach_neighbor_coord
 */
macro2 foreach_neighbor_coord(int r, coord c) {
  int ig = 0, jg = 0, kg = 0;
  NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
  coord d = c;
  coord_periodic_boundary(d);
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
  coord_periodic_boundary(d);
#if TREE
  Point point = locate_nonlocal(d.x,d.y,d.z);
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
  Point point = locate_nonlocal(d.x,d.y,d.z);
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
  coord_periodic_boundary(d);
#if TREE
  int effective_l = l;
  if (effective_l > depth())
    effective_l = depth();
  Point point = locate_level(d.x,d.y,d.z,effective_l);
  if (point.level >= 0) {
    foreach_neighbor(r) {
      if (allocated(0) && is_leaf(cell))
      // clang-format off
      {...}
      // clang-format on
    }
  }
#else 
  Point point = locate_level(d.x,d.y,d.z,l);
  if (point.level >= 0) {
    foreach_neighbor(r) {
      // clang-format off
      {...}
      // clang-format on
    }
  }
#endif 

}
