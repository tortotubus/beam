#pragma once

#include <math.h>

static inline double ib_locate_periodic_coord (double x,
                                               double x0,
                                               bool is_periodic) {
  return is_periodic ? x0 + fmod (fmod (x - x0, L0) + L0, L0) : x;
}

#if TREE

/**
 * @brief 
 */
trace Point locate_nonlocal (double xp = 0., double yp = 0., double zp = 0.) {
  xp = ib_locate_periodic_coord (xp, X0, Period.x);
#if dimension >= 2
  yp = ib_locate_periodic_coord (yp, Y0, Period.y);
#endif
#if dimension >= 3
  zp = ib_locate_periodic_coord (zp, Z0, Period.z);
#endif

  for (int l = depth (); l >= 0; l--) {
    Point point = {0};
    point.level = l;
    int n = 1 << point.level;
    point.i = (xp - X0) / L0 * n + GHOSTS;
#if dimension >= 2
    point.j = (yp - Y0) / L0 * n + GHOSTS;
#endif
#if dimension >= 3
    point.k = (zp - Z0) / L0 * n + GHOSTS;
#endif
    if (point.i >= 0 && point.i < n + 2 * GHOSTS
#if dimension >= 2
        && point.j >= 0 && point.j < n + 2 * GHOSTS
#endif
#if dimension >= 3
        && point.k >= 0 && point.k < n + 2 * GHOSTS
#endif
    ) {
      if (allocated (0) && is_leaf (cell))
        return point;
    } else
      break;
  }
  Point point = {0};
  point.level = -1;
  return point;
}

#else    // !TREE

/**
 * @brief
 */
trace Point locate_nonlocal (double xp = 0., double yp = 0., double zp = 0.) {
  xp = ib_locate_periodic_coord (xp, X0, Period.x);
#if dimension >= 2
  yp = ib_locate_periodic_coord (yp, Y0, Period.y);
#endif
#if dimension >= 3
  zp = ib_locate_periodic_coord (zp, Z0, Period.z);
#endif

  Point point = {0};
  point.level = depth ();
  SET_DIMENSIONS ();
  point.level = -1;
#if _MPI // fixme: can probably simplify with below
  point.i = (int) floor ((xp - X0) / L0 * point.n.x * Dimensions.x + GHOSTS -
                         mpi_coords[0] * point.n.x);
  if (point.i < 0 || point.i >= point.n.x + 2 * GHOSTS)
    return point;
#if dimension >= 2
  point.j = (int) floor ((yp - Y0) / L0 * point.n.x * Dimensions.x + GHOSTS -
                         mpi_coords[1] * point.n.x);
  if (point.j < 0 || point.j >= point.n.y + 2 * GHOSTS)
    return point;
#endif
#if dimension >= 3
  point.k = (int) floor ((zp - Z0) / L0 * point.n.x * Dimensions.x + GHOSTS -
                         mpi_coords[2] * point.n.x);
  if (point.k < 0 || point.k >= point.n.z + 2 * GHOSTS)
    return point;
#endif
#else // !_MPI
  point.i = (int) floor ((xp - X0) / L0 * point.n.x + GHOSTS);
  if (point.i < 0 || point.i >= point.n.x + 2 * GHOSTS)
    return point;
#if dimension >= 2
  point.j = (int) floor ((yp - Y0) / L0 * point.n.x + GHOSTS);
  if (point.j < 0 || point.j >= point.n.y + 2 * GHOSTS)
    return point;
#endif
#if dimension >= 3
  point.k = (int) floor ((zp - Z0) / L0 * point.n.x + GHOSTS);
  if (point.k < 0 || point.k >= point.n.z + 2 * GHOSTS)
    return point;
#endif
#endif // !_MPI
  point.level = depth ();
  return point;
}

#endif // !TREE

#if TREE

/**
 * @brief
 */
trace Point locate_level (double xp = 0.,
                          double yp = 0.,
                          double zp = 0.,
                          int level) {
  xp = ib_locate_periodic_coord (xp, X0, Period.x);
#if dimension >= 2
  yp = ib_locate_periodic_coord (yp, Y0, Period.y);
#endif
#if dimension >= 3
  zp = ib_locate_periodic_coord (zp, Z0, Period.z);
#endif

  {
    Point point = {0};
    point.level = level;
    int n = 1 << point.level;
    point.i = (xp - X0) / L0 * n + GHOSTS;
#if dimension >= 2
    point.j = (yp - Y0) / L0 * n + GHOSTS;
#endif
#if dimension >= 3
    point.k = (zp - Z0) / L0 * n + GHOSTS;
#endif
    if (level >= 0 && level <= depth () && point.i >= 0 &&
        point.i < n + 2 * GHOSTS
#if dimension >= 2
        && point.j >= 0 && point.j < n + 2 * GHOSTS
#endif
#if dimension >= 3
        && point.k >= 0 && point.k < n + 2 * GHOSTS
#endif
    ) {
      return point;
    }
  }
  Point point = {0};
  point.level = -1;
  return point;
}

#else    // !TREE

/**
 * @brief 
 */
trace Point locate_level (double xp = 0.,
                          double yp = 0.,
                          double zp = 0.,
                          int level) {
  xp = ib_locate_periodic_coord (xp, X0, Period.x);
#if dimension >= 2
  yp = ib_locate_periodic_coord (yp, Y0, Period.y);
#endif
#if dimension >= 3
  zp = ib_locate_periodic_coord (zp, Z0, Period.z);
#endif

  Point point = {0};
  point.level = depth ();
  SET_DIMENSIONS ();
  point.level = -1;
#if _MPI // fixme: can probably simplify with below
  point.i = (int) floor ((xp - X0) / L0 * point.n.x * Dimensions.x + GHOSTS -
                         mpi_coords[0] * point.n.x);
  if (point.i < 0 || point.i >= point.n.x + 2 * GHOSTS)
    return point;
#if dimension >= 2
  point.j = (int) floor ((yp - Y0) / L0 * point.n.x * Dimensions.x + GHOSTS -
                         mpi_coords[1] * point.n.x);
  if (point.j < 0 || point.j >= point.n.y + 2 * GHOSTS)
    return point;
#endif
#if dimension >= 3
  point.k = (int) floor ((zp - Z0) / L0 * point.n.x * Dimensions.x + GHOSTS -
                         mpi_coords[2] * point.n.x);
  if (point.k < 0 || point.k >= point.n.z + 2 * GHOSTS)
    return point;
#endif
#else // !_MPI
  point.i = (int) floor ((xp - X0) / L0 * point.n.x + GHOSTS);
  if (point.i < 0 || point.i >= point.n.x + 2 * GHOSTS)
    return point;
#if dimension >= 2
  point.j = (int) floor ((yp - Y0) / L0 * point.n.x + GHOSTS);
  if (point.j < 0 || point.j >= point.n.y + 2 * GHOSTS)
    return point;
#endif
#if dimension >= 3
  point.k = (int) floor ((zp - Z0) / L0 * point.n.x + GHOSTS);
  if (point.k < 0 || point.k >= point.n.z + 2 * GHOSTS)
    return point;
#endif
#endif // !_MPI
  point.level = depth ();
  return point;
}
#endif
