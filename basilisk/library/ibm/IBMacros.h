#include "library/ibm/IBLocate.h"

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
 * @def foreach_neighborhood_coord
 *
 * @brief Loops through the L1 neighborhood of cells point centered on
 * coordinate @c c
 */
macro2 foreach_neighborhood_coord (coord c, size_t radius) {

  // Wrap the coordinate according to any periodic boundary conditions
  coord centre_coord = c;
  coord_periodic_boundary (centre_coord);
  Point centre_point = locate (centre_coord.x, centre_coord.y, centre_coord.z);

  // Set point level to the located Point
  Point point = {0};
  Point iterpoint = {0};

  point.level = centre_point.level;
  iterpoint.level = centre_point.level;

  if (iterpoint.level != -1) {
#if dimension == 1
    for (iterpoint.i = -radius + centre_point.i;
         iterpoint.i <= radius + centre_point.i;
         iterpoint.i++) {
      point = iterpoint;
      if (!allocated (0, 0, 0)) {
        continue;
      } else {
        // point_periodic_boundary (point);
        int ig = point.i;
        // clang-format off
      {...}
        // clang-format on
      }
    }
#elif dimension == 2
    for (iterpoint.i = -radius + centre_point.i;
         iterpoint.i <= radius + centre_point.i;
         iterpoint.i++) {
      for (iterpoint.j = -radius + centre_point.j;
           iterpoint.j <= radius + centre_point.j;
           iterpoint.j++) {
        point.i = iterpoint.i;
        point.j = iterpoint.j;
        if (!allocated (0, 0, 0)) {
          continue;
        } else {
          // point_periodic_boundary (point);
          int ig = point.i;
          int jg = point.j;
          // clang-format off
          {...}
          // clang-format on
        }
      }
    }
#else // dimension == 3
    for (iterpoint.i = -radius + centre_point.i;
         iterpoint.i <= radius + centre_point.i;
         iterpoint.i++) {
      for (iterpoint.j = -radius + centre_point.j;
           iterpoint.j <= radius + centre_point.j;
           point.j++) {
        for (iterpoint.k = -radius + centre_point.k;
             iterpoint.k <= radius + centre_point.k;
             iterpoint.k++) {
          point = iterpoint;
          if (!allocated (0, 0, 0)) {
            continue;
          } else {
            // point_periodic_boundary (point);
            int ig = point.i;
            int jg = point.j;
            int kg = point.k;
            // clang-format off
            {...}
            // clang-format on
          }
        }
      }
    }
#endif
  }
}

macro2 foreach_neighborhood_coord_nonlocal (coord c, size_t radius) {

  coord centre_coord = c;
  coord_periodic_boundary (centre_coord);
  Point centre_point = locate_nonlocal (centre_coord.x, centre_coord.y, centre_coord.z);

  // Set point level to the located Point
  Point point = {0};
  Point iterpoint = {0};

  point.level = centre_point.level;
  iterpoint.level = centre_point.level;

  if (iterpoint.level > 0) {
#if dimension == 1
    for (iterpoint.i = -radius + centre_point.i;
         iterpoint.i <= radius + centre_point.i;
         iterpoint.i++) {
      point = iterpoint;
      if (!allocated (0, 0, 0) || !is_local (cell)) {
        continue;
      } else {
        // point_periodic_boundary (point);
        int ig = point.i;
        // clang-format off
      {...}
        // clang-format on
      }
    }
#elif dimension == 2
    for (iterpoint.i = -radius + centre_point.i;
         iterpoint.i <= radius + centre_point.i;
         iterpoint.i++) {
      for (iterpoint.j = -radius + centre_point.j;
           iterpoint.j <= radius + centre_point.j;
           iterpoint.j++) {
        point.i = iterpoint.i;
        point.j = iterpoint.j;
        if (!allocated (0, 0, 0) || !is_local (cell)) {
          continue;
        } else {
          // point_periodic_boundary (point);
          int ig = point.i;
          int jg = point.j;
          // clang-format off
          {...}
          // clang-format on
        }
      }
    }
#else // dimension == 3
    for (iterpoint.i = -radius + centre_point.i;
         iterpoint.i <= radius + centre_point.i;
         iterpoint.i++) {
      for (iterpoint.j = -radius + centre_point.j;
           iterpoint.j <= radius + centre_point.j;
           point.j++) {
        for (iterpoint.k = -radius + centre_point.k;
             iterpoint.k <= radius + centre_point.k;
             iterpoint.k++) {
          point = iterpoint;
          if (!allocated (0, 0, 0)) {
            continue;
          } else {
            // point_periodic_boundary (point);
            int ig = point.i;
            int jg = point.j;
            int kg = point.k;
            // clang-format off
            {...}
            // clang-format on
          }
        }
      }
    }
#endif
  }
}

/**
 * @def foreach_neighborhood_coord
 *
 * @brief Loops through the L1 neighborhood of cells point centered on
 * coordinate @c c at level @c level where it is up to the user to check
 * that a point exists
 */
macro2 foreach_neighborhood_coord_level (coord c, size_t radius, int lvl) {

  coord centre_coord = c;
  coord_periodic_boundary (centre_coord);
  Point centre_point =
    locate_level (centre_coord.x, centre_coord.y, centre_coord.z, lvl);

  Point point = {0};
  Point iterpoint = {0};

  point.level = centre_point.level;
  iterpoint.level = centre_point.level;

  if (iterpoint.level != -1) {
#if dimension == 1
    for (iterpoint.i = -radius + centre_point.i;
         iterpoint.i <= radius + centre_point.i;
         iterpoint.i++) {
      point.i = iterpoint.i;
      point.level = centre_point.level;
      if (allocated (0) && is_leaf (cell)) {
        // point_periodic_boundary (point);
        int ig = point.i;
        // clang-format off
        {...}
        // clang-format on
      }
    }
#elif dimension == 2
    for (iterpoint.i = -radius + centre_point.i;
         iterpoint.i <= radius + centre_point.i;
         iterpoint.i++) {
      for (iterpoint.j = -radius + centre_point.j;
           iterpoint.j <= radius + centre_point.j;
           iterpoint.j++) {
        point.i = iterpoint.i;
        point.j = iterpoint.j;
        point.level = centre_point.level;
        if (allocated (0) && is_leaf (cell)) {
          // point_periodic_boundary (point);
          int ig = point.i;
          int jg = point.j;
          // clang-format off
          {...}
          // clang-format on
        }
      }
    }
#else // dimension == 3
    for (iterpoint.i = -radius + centre_point.i;
         iterpoint.i <= radius + centre_point.i;
         iterpoint.i++) {
      for (iterpoint.j = -radius + centre_point.j;
           iterpoint.j <= radius + centre_point.j;
           point.j++) {
        for (iterpoint.k = -radius + centre_point.k;
             iterpoint.k <= radius + centre_point.k;
             iterpoint.k++) {
          point.i = iterpoint.i;
          point.j = iterpoint.j;
          point.k = iterpoint.k;
          point.level = centre_point.level;
          if (allocated (0) && is_leaf (cell)) {
            // point_periodic_boundary (point);
            int ig = point.i;
            int jg = point.j;
            int kg = point.k;
            // clang-format off
            {...}
            // clang-format on
          }
        }
      }
    }
#endif
  }
}
