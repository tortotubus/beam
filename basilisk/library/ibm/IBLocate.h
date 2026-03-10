// #if TREE
// trace Point locater (double xp = 0., double yp = 0., double zp = 0.) {
//   return locate (xp, yp, zp);
// }
// #else // multigrid
trace Point locater (double xp = 0., double yp = 0., double zp = 0.) {
  Point point = {0};
  point.level = depth ();
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
    return point;
  } else {
    point.level = -1;
    return point;
  }
}
// #endif

trace Point locate_nonlocal (double xp = 0., double yp = 0., double zp = 0.) {
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

trace Point locate_level (double xp = 0.,
                          double yp = 0.,
                          double zp = 0.,
                          int level) {
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
