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
  // printf ("Periodic coord: %f %f %f\n",
  //         centre_coord.x,
  //         centre_coord.y,
  //         centre_coord.z);
  Point centre_point = locate (centre_coord.x, centre_coord.y, centre_coord.z);
  // printf ("Periodic point: %d %d\n", centre_point.i, centre_point.j);

  // Set point level to the located Point
  Point point = {0};
  point.level = centre_point.level;

  if (point.level != -1) {
#if dimension == 1
    for (point.i = -radius + centre.i; point.i <= radius + centre.i;
         point.i++) {
      if (!allocated (0, 0, 0)) {
        fprintf (stderr, "warning: Neighborhood not fully resolved\n");
      } else {
        // point_periodic_boundary (point);
        int ig = point.i;
        // clang-format off
      {...}
        // clang-format on
      }
    }
#elif dimension == 2
    for (point.i = -radius + centre_point.i; point.i <= radius + centre_point.i;
         point.i++) {
      for (point.j = -radius + centre_point.j;
           point.j <= radius + centre_point.j;
           point.j++) {
        if (!allocated (0, 0, 0)) {
          fprintf (stderr,
                   "warning: Neighborhood not fully resolved: %d %d\n",
                   point.i,
                   point.j);
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
    for (point.i = -radius + centre_point.i; point.i <= radius + centre_point.i;
         point.i++) {
      for (point.j = -radius + centre_point.j;
           point.j <= radius + centre_point.j;
           point.j++) {
        for (point.k = -radius + centre_point.k;
             point.k <= radius + centre_point.k;
             point.k++) {
          if (!allocated (0, 0, 0)) {
            fprintf (stderr, "warning: Neighborhood not fully resolved\n");
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

  // Wrap the coordinate according to any periodic boundary conditions
  coord centre_coord = c;
  coord_periodic_boundary (centre_coord);
  // printf ("Periodic coord: %f %f %f\n",
  //         centre_coord.x,
  //         centre_coord.y,
  //         centre_coord.z);

  int n = 1 << lvl;
  Point centre_point = {
    .level = lvl,
#if dimension == 1
    .i = (centre_coord.x - X0) / L0 * n + GHOSTS,
#elif dimension == 2
    .i = (centre_coord.x - X0) / L0 * n + GHOSTS,
    .j = (centre_coord.y - Y0) / L0 * n + GHOSTS
#else
    .i = (centre_coord.x - X0) / L0 * n + GHOSTS,
    .j = (centre_coord.y - Y0) / L0 * n + GHOSTS,
    .k = (centre_coord.z - Z0) / L0 * n + GHOSTS
#endif
  };

  // printf ("Periodic point: %d %d\n", centre_point.i, centre_point.j);

  // Set point level to the located Point
  Point point = {0};
  point.level = centre_point.level;

#if dimension == 1
  for (point.i = -radius + centre.i; point.i <= radius + centre.i; point.i++) {
    int ig = point.i;
    // clang-format off
    {...}
    // clang-format on
  }
#elif dimension == 2
  for (point.i = -radius + centre_point.i; point.i <= radius + centre_point.i;
       point.i++) {
    for (point.j = -radius + centre_point.j; point.j <= radius + centre_point.j;
         point.j++) {
      // point_periodic_boundary (point);
      int ig = point.i;
      int jg = point.j;
      // clang-format off
      {...}
      // clang-format on
    }
  }
#else // dimension == 3
  for (point.i = -radius + centre_point.i; point.i <= radius + centre_point.i;
       point.i++) {
    for (point.j = -radius + centre_point.j; point.j <= radius + centre_point.j;
         point.j++) {
      for (point.k = -radius + centre_point.k;
           point.k <= radius + centre_point.k;
           point.k++) {
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
#endif
}