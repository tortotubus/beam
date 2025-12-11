#define BGHOSTS 2

trace void synchronize (scalar* list) {
  for (scalar s in list)
    s.dirty = true;
  boundary (list);
}

// clang-format off
#define ACROSS_PERIODIC(a, b) (fabs(a - b) > L0 / 2.)
#define PERIODIC_1DIST(a, b) (fabs(a - L0 - b) > L0 / 2. ? a + L0 - b : a - L0 - b)
#define GENERAL_1DIST(a, b) (ACROSS_PERIODIC(a, b) ? PERIODIC_1DIST(a, b) : a - b)
#define PERIODIC_1DAVG(a, b) (fabs(a - L0 - b) > L0 / 2. ? a + L0 + b : a - L0 + b)
#define GENERAL_1DAVG(a, b) (ACROSS_PERIODIC(a, b) ? PERIODIC_1DAVG(a, b) : a + b)
#define GENERAL_SQNORM(a, b) (sq(GENERAL_1DIST(a.x, b.x)) + sq(GENERAL_1DIST(a.y, b.y)) + sq(GENERAL_1DIST(a.z, b.z)))

#define POS_PBC_X(X) ((!Period.x) ? (X) : (((X - (X0 + L0 / 2)) > L0 / 2.) ? (X) - L0 : (X)))
#define POS_PBC_Y(Y) ((!Period.y) ? (Y) : (((Y - (Y0 + L0 / 2)) > L0 / 2.) ? (Y) - L0 : (Y)))
#define POS_PBC_Z(Z) ((!Period.z) ? (Z) : (((Z - (Z0 + L0 / 2)) > L0 / 2.) ? (Z) - L0 : (Z)))

// clang-format on


#define PESKIN_SUPPORT_RADIUS 2

// macro peskin_kernel (coord delta_centre, Point point = point) {
//   coord dist; 
// #if dimension == 1
// #elif dimension == 2 
//   dist.x = x - delta_centre.x;
//   dist.y = y - delta_centre.y;
//   bool cutoff = fabs(dist.x) <= Delta * 2 && fabs(dist.y) <= Delta * 2;
//   double peskin_weight = 0.;
//   if (cutoff) {
//     peskin_weight = (1. + cos(0.5 * pi * dist.x / Delta)) * (1. + cos(0.5 * pi * dist.x / Delta)) / 16.;
//   } else {
//     peskin_weight = 0.;
//   }
// #else 
// #endif 
// }

macro peskin_kernel  (coord delta_centre, Point point = point) {
  if (point.level >= 0) {
    coord dist;
#if dimension == 1
    dist.x = GENERAL_1DIST (x, delta_centre.x);
    if (fabs (dist.x) <= Delta * PESKIN_SUPPORT_RADIUS) {
      double peskin_weight =
        (1. + cos (pi * dist.x / (Delta * PESKIN_SUPPORT_RADIUS))) /
        (2. * Delta * PESKIN_SUPPORT_RADIUS);
      NOT_UNUSED (peskin_weight);
      // clang-format off
      {...}
      // clang-format on
    }
#elif dimension == 2
    dist.x = GENERAL_1DIST (x, delta_centre.x);
    dist.y = GENERAL_1DIST (y, delta_centre.y);
    if (fabs (dist.x) <= Delta * PESKIN_SUPPORT_RADIUS &&
        fabs (dist.y) <= Delta * PESKIN_SUPPORT_RADIUS) {
      double peskin_weight =
        (1. + cos (pi * dist.x / (Delta * PESKIN_SUPPORT_RADIUS))) *
        (1. + cos (pi * dist.y / (Delta * PESKIN_SUPPORT_RADIUS))) /
        sq (2. * PESKIN_SUPPORT_RADIUS);
      NOT_UNUSED (peskin_weight);
      // clang-format off
      {...}
      // clang-format on
    }
#else
    dist.x = GENERAL_1DIST (x, delta_centre.x);
    dist.y = GENERAL_1DIST (y, delta_centre.y);
    dist.z = GENERAL_1DIST (z, delta_centre.z);
    if (fabs (dist.x) <= Delta * PESKIN_SUPPORT_RADIUS &&
        fabs (dist.y) <= Delta * PESKIN_SUPPORT_RADIUS &&
        fabs (dist.z) <= Delta * PESKIN_SUPPORT_RADIUS) {
      double peskin_weight =
        (1. + cos (pi * dist.x / (Delta * PESKIN_SUPPORT_RADIUS))) *
        (1. + cos (pi * dist.y / (Delta * PESKIN_SUPPORT_RADIUS))) *
        (1. + cos (pi * dist.z / (Delta * PESKIN_SUPPORT_RADIUS))) /
        cube (2. * Delta * PESKIN_SUPPORT_RADIUS);
      NOT_UNUSED (peskin_weight);
      // clang-format off
      {...}
      // clang-format on
    }
#endif
  }
}

// Old locate-based peskin stencil

macro foreach_neighbor_of_coord3 (coord c, unsigned int size) {
  {
    // int ig = 0;
    // NOT_UNUSED(ig);
    // int jg = 0;
    // NOT_UNUSED(jg);
    // int kg = 0;
    // NOT_UNUSED(kg);

    double delta_smallest = (L0 / (1 << grid->maxdepth));

#if dimension == 1
    for (int _k = -size; _k <= size; _k++) {
      Point point = locate (POS_PBC_X (c.x + _k * delta_smallest));

      if (point.level == -1) {
        // fprintf (stderr,
        //          "Coord (%f) is located outside of grid across "
        //          "non-periodic boundary.\n",
        //          c.x + _k * delta_smallest);
      } else {
        // clang-format off
      {...}
        // clang-format on
      }
    }
#elif dimension == 2
    for (int _k = -size; _k <= size; _k++) {
      for (int _l = -size; _l <= size; _l++) {
        Point point = locate (POS_PBC_X (c.x + _k * delta_smallest),
                              POS_PBC_Y (c.y + _l * delta_smallest));

        if (point.level == -1) {
          // fprintf (stderr,
          //          "Coord (%f, %f) is located outside of grid across "
          //          "non-periodic boundary.\n",
          //          c.x + _k * delta_smallest,
          //          c.y + _l * delta_smallest);
        } else {
          // clang-format off
        {...}
          // clang-format on
        }
      }
    }
#else // dimension == 3
    for (int _l = -size; _l <= size; _l++) {
      for (int _m = -size; _m <= size; _m++) {
        for (int _n = -size; _n <= size; _n++) {
          Point point = locate (POS_PBC_X (c.x + _l * delta_smallest),
                                POS_PBC_Y (c.y + _m * delta_smallest),
                                POS_PBC_Z (c.z + _n * delta_smallest));
          if (point.level == -1) {
            // fprintf (stderr,
            //          "Coord (%f, %f, %f) is located outside of grid across "
            //          "non-periodic boundary.\n",
            //          c.x + _l * delta_smallest,
            //          c.y + _m * delta_smallest,
            //          c.z + _n * delta_smallest);
          } else {
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


macro2 foreach_neighbor_of_coord2 (coord c, unsigned int size) {
  {
    int ig = 0;
    NOT_UNUSED (ig);
    int jg = 0;
    NOT_UNUSED (jg);
    int kg = 0;
    NOT_UNUSED (kg);

#if dimension == 1
    Point point = locate (POS_PBC_X (c.x));

    if (point.level == -1) {
      fprintf (stderr,
               "Coord (%f) is located outside of grid across "
               "non-periodic boundary.\n",
               c.x);
      // abort ();
    } else {
      // foreach_neighbor (size) {
      //   // clang-format off
      // {...}
      //   // clang-format on
      // }
    }
#elif dimension == 2
    Point point = locate (POS_PBC_X (c.x), POS_PBC_Y (c.y));

    if (point.level == -1) {
      fprintf (stderr,
               "Coord (%f, %f) is located outside of grid across "
               "non-periodic boundary.\n",
               c.x,
               c.y);
      // abort ();
    } else {
      foreach_neighbor (size) {
        // clang-format off
    {...}
        // clang-format on
      }
    }
#else // dimension == 3
    Point point = locate (POS_PBC_X (c.x), POS_PBC_Y (c.y), POS_PBC_Z (c.z));
    if (point.level == -1) {
      fprintf (stderr,
               "Coord (%f, %f, %f) is located outside of grid across "
               "non-periodic boundary.\n",
               c.x,
               c.y,
               c.z);
      // abort ();
    } else {
    }
#endif
  }
}

macro2 foreach_neighbor_of_coord  (coord c, unsigned int size) {
  {
    int ig = 0;
    NOT_UNUSED (ig);
    int jg = 0;
    NOT_UNUSED (jg);
    int kg = 0;
    NOT_UNUSED (kg);

    Point point = {0};

    const int _lc = depth();
    const int _nc = 1 << _lc;

    const int _ic = (c.x - X0)/L0*_nc + GHOSTS;
    const int _jc = (c.y - Y0)/L0*_nc + GHOSTS;
    const int _kc = (c.z - Z0)/L0*_nc + GHOSTS;

    for (int _in = -size; _in <= size; _in++) {
      for (int _jn = -size; _jn <= size; _jn++) {

        point.i = _ic + _in;
        point.j = _jc + _jn;
        point.level = _lc;

        if (allocated (0, 0, 0) && is_local(cell)) {
          // clang-format off
          {...}
          // clang-format on
        }
      }
    }
  }
}