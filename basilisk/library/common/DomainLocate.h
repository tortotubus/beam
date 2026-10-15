#pragma once

/*!
 * @file DomainLocate.h
 *
 * Shared coordinate-to-Point lookup helpers. These mirror the locate helpers
 * currently carried by both particle and IBM libraries, but use common names so
 * existing callsites can migrate gradually.
 */

#include <math.h>

#include "library/common/DomainGeometry.h"

/*!
 * @brief Interpolate from a scalar using only local cells.
 *
 * This intentionally uses @c locate() and never performs the MPI gather used by
 * Basilisk @c interpolate().
 */
trace double domain_interpolate_local (scalar s,
                                       double xp = 0.,
                                       double yp = 0.,
                                       double zp = 0.,
                                       bool linear = true) {
  double val = nodata;
  Point point = locate (xp, yp, zp);
  if (point.level < 0)
    return val;

  int ig = 0, jg = 0, kg = 0;
  NOT_UNUSED (ig);
  NOT_UNUSED (jg);
  NOT_UNUSED (kg);
  POINT_VARIABLES ();

  return linear ? interpolate_linear (point, s, xp, yp, zp) : s[];
}

#if TREE

/*!
 * @brief Locate a leaf containing a coordinate, including non-local tree cells
 * replicated as MPI ghosts.
 */
trace Point domain_locate_nonlocal (double xp = 0.,
                                    double yp = 0.,
                                    double zp = 0.) {
  for (int l = depth (); l >= 0; l--) {
    Point point = {0};
    point.level = l;
    const int n = 1 << point.level;
    point.i = (int) floor ((xp - X0)/domain_extent_x ()*n + GHOSTS);
#if dimension >= 2
    point.j = (int) floor ((yp - Y0)/domain_extent_y ()*n + GHOSTS);
#endif
#if dimension >= 3
    point.k = (int) floor ((zp - Z0)/domain_extent_z ()*n + GHOSTS);
#endif
    if (point.i >= 0 && point.i < n + 2*GHOSTS
#if dimension >= 2
        && point.j >= 0 && point.j < n + 2*GHOSTS
#endif
#if dimension >= 3
        && point.k >= 0 && point.k < n + 2*GHOSTS
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

#else

/*!
 * @brief Locate a coordinate on a multigrid/cartesian grid, allowing local MPI
 * ghost cells.
 */
trace Point domain_locate_nonlocal (double xp = 0.,
                                    double yp = 0.,
                                    double zp = 0.) {
  Point point = {0};
  point.level = depth ();
  SET_DIMENSIONS ();
  point.level = -1;
#if _MPI
  point.i = (int) floor ((xp - X0)/domain_extent_x ()*point.n.x*Dimensions.x +
                         GHOSTS - mpi_coords[0]*point.n.x);
#else
  point.i = (int) floor ((xp - X0)/domain_extent_x ()*point.n.x + GHOSTS);
#endif
  if (point.i < 0 || point.i >= point.n.x + 2*GHOSTS)
    return point;
#if dimension >= 2
#if _MPI
  point.j = (int) floor ((yp - Y0)/domain_extent_y ()*point.n.y*Dimensions.y +
                         GHOSTS - mpi_coords[1]*point.n.y);
#else
  point.j = (int) floor ((yp - Y0)/domain_extent_y ()*point.n.y + GHOSTS);
#endif
  if (point.j < 0 || point.j >= point.n.y + 2*GHOSTS)
    return point;
#endif
#if dimension >= 3
#if _MPI
  point.k = (int) floor ((zp - Z0)/domain_extent_z ()*point.n.z*Dimensions.z +
                         GHOSTS - mpi_coords[2]*point.n.z);
#else
  point.k = (int) floor ((zp - Z0)/domain_extent_z ()*point.n.z + GHOSTS);
#endif
  if (point.k < 0 || point.k >= point.n.z + 2*GHOSTS)
    return point;
#endif
  point.level = depth ();
  return point;
}

#endif

#if TREE

/*!
 * @brief Locate the Point at a requested tree level.
 */
trace Point domain_locate_level (double xp = 0.,
                                 double yp = 0.,
                                 double zp = 0.,
                                 int level) {
  Point point = {0};
  point.level = level;
  const int n = 1 << point.level;
  point.i = (int) floor ((xp - X0)/domain_extent_x ()*n + GHOSTS);
#if dimension >= 2
  point.j = (int) floor ((yp - Y0)/domain_extent_y ()*n + GHOSTS);
#endif
#if dimension >= 3
  point.k = (int) floor ((zp - Z0)/domain_extent_z ()*n + GHOSTS);
#endif
  if (level >= 0 && level <= depth () && point.i >= 0 &&
      point.i < n + 2*GHOSTS
#if dimension >= 2
      && point.j >= 0 && point.j < n + 2*GHOSTS
#endif
#if dimension >= 3
      && point.k >= 0 && point.k < n + 2*GHOSTS
#endif
  )
    return point;

  point.level = -1;
  return point;
}

#else

/*!
 * @brief Locate a coordinate on a fixed grid. For multigrid this is equivalent
 * to @ref domain_locate_nonlocal because the requested level is ignored by the
 * underlying fixed grid.
 */
trace Point domain_locate_level (double xp = 0.,
                                 double yp = 0.,
                                 double zp = 0.,
                                 int level) {
  NOT_UNUSED (level);
  return domain_locate_nonlocal (xp, yp, zp);
}

#endif
