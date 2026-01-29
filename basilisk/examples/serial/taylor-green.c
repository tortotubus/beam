/**
# Taylor--Green vortices

[Taylor--Green
vortices](http://en.wikipedia.org/wiki/Taylor%E2%80%93Green_vortex)
are one of the few exact non-trivial solutions of the incompressible
Euler equations. In this test case, we use this solution as initial
condition and check whether the numerical scheme can respect the
balance between non-linear advection terms and pressure
gradients. Numerical diffusion will in particular introduce
dissipation. This dissipation can be quantified and is a useful
measure of the accuracy of the numerical scheme.

We solve the incompressible Euler equations on a Cartesian (multi)grid
either with the centered Navier-Stokes solver or with the "all Mach"
solver (in incompressible mode). */

#include "grid/multigrid.h"
  #if ALL_MACH
  # include "all-mach.h"
  # include "bcg.h"
  # define u q

  event tracer_advection (i++)
    advection ((scalar *){q}, uf, dt, (scalar *){g});
#else
  #include "library/ibm/navier-stokes/centered-splitting.h" 
#endif


int
main()
{
  /**
  Space and time are dimensionless. The domain is unity, centered on
  the origin and periodic in all directions. */

  size(1. [0]);
  DT = HUGE[0];
  origin(-0.5, -0.5);
  foreach_dimension()
    periodic(right);

  /**
  We check convergence with spatial resolution from 32^2^ to
  256^2^. */
  for (N = 32; N <= 256; N *= 2)
    run();
}

event
init(i = 0)
{
  /**
  This is the initial Taylor--Green solution for velocity and
  pressure. */
  
  foreach () {
    u.x[] = -cos(2. * pi * x) * sin(2. * pi * y);
    u.y[] = sin(2. * pi * x) * cos(2. * pi * y);
    p[] = -(cos(4. * pi * x) + cos(4. * pi * y)) / 4.;
  } 

  /**
  We also need to define the initial centered pressure gradient (this
  improves the accuracy of the initial conditions). */

  // foreach()
  //   foreach_dimension()
  //     g.x[] = - (p[1] - p[-1])/(2.*Delta);
}

event logfile (i++) {

  /**
  We log the evolution of the maximum divergence and of the total
  kinetic energy. */
  
  scalar div[], ke[];
  foreach() {
    div[] = (u.x[1,0] - u.x[-1,0] + u.y[0,1] - u.y[0,-1])/(2.*Delta);
    ke[] = sq(u.x[]) + sq(u.y[]);
  }
  printf ("%d %d %g %g %g\n", N, i, t, normf(div).max, statsf(ke).sum);
}

event
error(t = 2)
{ 
  /**
  At $t=2$ we compute the error on the norm of the velocity. */
  
  scalar e[];
  foreach () {
    double u0 = -cos(2. * pi * x) * sin(2. * pi * y);
    double v0 = sin(2. * pi * x) * cos(2. * pi * y);
    e[] = norm(u) - sqrt(sq(u0) + sq(v0));
  }
  norm n = normf(e);
  fprintf(stderr, "%d %g %g %g\n", N, n.avg, n.rms, n.max);
}
