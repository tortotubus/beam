#pragma once

/*!
 * @file DomainGeometry.h
 *
 * Shared coordinate/domain helpers for particle-like and immersed-boundary
 * objects. These utilities are intentionally independent of Particle and IBNode
 * storage so both libraries can migrate to them incrementally.
 */

#include <math.h>

/*!
 * @brief Domain extent in x.
 */
static inline double domain_extent_x () {
  return L0;
}

#if dimension >= 2
/*!
 * @brief Domain extent in y, including Basilisk multigrid aspect ratio.
 */
static inline double domain_extent_y () {
  return L0*(double) Dimensions.y/(double) Dimensions.x;
}
#endif

#if dimension >= 3
/*!
 * @brief Domain extent in z, including Basilisk multigrid aspect ratio.
 */
static inline double domain_extent_z () {
  return L0*(double) Dimensions.z/(double) Dimensions.x;
}
#endif

/*!
 * @brief Wrap a coordinate into periodic directions.
 */
macro domain_wrap_coord (coord c) {
  if (Period.x) {
    const double Lx = domain_extent_x ();
    c.x = X0 + fmod (fmod (c.x - X0, Lx) + Lx, Lx);
  }
#if dimension >= 2
  if (Period.y) {
    const double Ly = domain_extent_y ();
    c.y = Y0 + fmod (fmod (c.y - Y0, Ly) + Ly, Ly);
  }
#endif
#if dimension >= 3
  if (Period.z) {
    const double Lz = domain_extent_z ();
    c.z = Z0 + fmod (fmod (c.z - Z0, Lz) + Lz, Lz);
  }
#endif
}

/*!
 * @brief Wrap a grid point index across periodic directions.
 */
macro domain_wrap_point (Point p) {
  const int N = (1 << p.level);
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

/*!
 * @brief Apply minimum-image convention to a displacement vector.
 */
macro domain_minimum_image_delta (coord dr) {
  if (Period.x) {
    const double Lx = domain_extent_x ();
    dr.x -= Lx*nearbyint (dr.x/Lx);
  }
#if dimension >= 2
  if (Period.y) {
    const double Ly = domain_extent_y ();
    dr.y -= Ly*nearbyint (dr.y/Ly);
  }
#endif
#if dimension >= 3
  if (Period.z) {
    const double Lz = domain_extent_z ();
    dr.z -= Lz*nearbyint (dr.z/Lz);
  }
#endif
}

/*!
 * @brief Compute minimum-image displacement from @p b to @p a.
 */
macro domain_displacement (coord dr, coord a, coord b) {
  foreach_dimension()
    dr.x = a.x - b.x;
  domain_minimum_image_delta (dr);
}

/*!
 * @brief Compute squared minimum-image distance between two coordinates.
 */
macro domain_distance2 (double r2, coord a, coord b) {
  coord domain_distance_dr = {0., 0., 0.};
  domain_displacement (domain_distance_dr, a, b);
  r2 = 0.;
  foreach_dimension()
    r2 += sq (domain_distance_dr.x);
}
