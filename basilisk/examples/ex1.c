#define LEVEL 6
#define LENGTH 2.
#define T_END 1.

#include "grid/quadtree.h"
#include "navier-stokes/centered.h"

// #include "../interface/ibm.h"

int
main()
{
  L0 = 1;
  origin(-.5 * L0, -.5 * L0);
  init_grid(1 << LEVEL);
  periodic(left);
  DT = 5.e-3;
  run();
}

u[]

event
impose_u(i++)
{
  foreach () {
    u.x[] = 1.;
    u.y[] = 0.;
  }
}

event
progress_output(i++)
{
  if (pid() == 0) {
  }
}

event
movies(i += 5; t <= T_END)
{
  scalar omega[];
  vorticity(u, omega);
  output_ppm(omega, file = "vort.mp4", linear = true);
}

event
output(t = T_END)
{
  if (pid() == 0) {
  }
}

event
end(t = T_END)
{
  return 0;
}