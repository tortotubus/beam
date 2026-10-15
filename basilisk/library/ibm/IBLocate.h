#pragma once

/* Compatibility wrappers for the shared domain-location helpers. IBM node
 * coordinates are material coordinates and may be outside the primary box;
 * only the copy used for Eulerian lookup is wrapped. */

#include "library/common/DomainLocate.h"

trace Point locate_nonlocal (double xp = 0.,
                             double yp = 0.,
                             double zp = 0.) {
  coord c = {.x = xp, .y = yp, .z = zp};
  domain_wrap_coord (c);
  return domain_locate_nonlocal (c.x, c.y, c.z);
}

trace Point locate_level (double xp = 0.,
                          double yp = 0.,
                          double zp = 0.,
                          int level) {
  coord c = {.x = xp, .y = yp, .z = zp};
  domain_wrap_coord (c);
  return domain_locate_level (c.x, c.y, c.z, level);
}
