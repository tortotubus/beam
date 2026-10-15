#include "run.h"
#include "common/timestep.h"
#include "bcg.h"
#if EMBED
#include "viscosity-embed.h"
#else
#include "viscosity.h"
#endif

#include "ibm/IBMeshManager.h"
#if TREE
#include "ibm/IBAdapt.h"
#endif

scalar p[];
vector u[], g[];
scalar pf[];
face vector uf[];

vector ibmf[];

#ifndef ADVECT_LAG_RK2
#define ADVECT_LAG_RK2 1
#endif

(const) face vector mu = zerof, a = zerof, alpha = unityf;
(const) scalar rho = unity;
bool stokes = false;

face vector mu_alpha[], mu_beta[];
mgstats mgp = {0}, mgpf = {0}, mgu = {0};

#if EMBED
#define neumann_pressure(i)                                                              \
  (alpha.n[i] ? a.n[i] * fm.n[i] / alpha.n[i] : a.n[i] * rho[] / (cm[] + SEPS))
#else
#define neumann_pressure(i) (alpha.n[i] ? a.n[i] * fm.n[i] / alpha.n[i] : 0.)
#endif

p[right] = neumann (neumann_pressure (ghost));
p[left] = neumann (-neumann_pressure (0));

#if AXI
uf.n[bottom] = 0.;
uf.t[bottom] = dirichlet (0); // since uf is multiplied by the metric which
                              // is zero on the axis of symmetry
p[top] = neumann (neumann_pressure (ghost));
#else // !AXI
#if dimension > 1
p[top] = neumann (neumann_pressure (ghost));
p[bottom] = neumann (-neumann_pressure (0));
#endif
#if dimension > 2
p[front] = neumann (neumann_pressure (ghost));
p[back] = neumann (-neumann_pressure (0));
#endif
#endif // !AXI

#if TREE && EMBED
void pressure_embed_gradient (Point point, scalar p, coord* g) {
  foreach_dimension () g->x = rho[] / (cm[] + SEPS) * (a.x[] + a.x[1]) / 2.;
}
#endif // TREE && EMBED

event defaults (i = 0) {
  init_ibsolver ();
  ibmeshmanager_init (0);

  if (is_constant (a.x)) {
    a = new face vector;
    foreach_face () a.x[] = 0.;
  }

  mgp = (mgstats) {0};
  mgpf = (mgstats) {0};
  mgu = (mgstats) {0};

  CFL = 0.8;

  p.nodump = pf.nodump = true;

  if (alpha.x.i == unityf.x.i) {
    alpha = fm;
    rho = cm;
  } else if (!is_constant (alpha.x)) {
    face vector alphav = alpha;
    foreach_face () alphav.x[] = fm.x[];
  }

#if TREE
  uf.x.refine = refine_face_solenoidal;

#if EMBED
  uf.x.refine = refine_face;
  foreach_dimension () uf.x.prolongation = refine_embed_face_x;
  for (scalar s in{p, pf, u, g}) {
    s.restriction = restriction_embed_linear;
    s.refine = s.prolongation = refine_embed_linear;
    s.depends = list_add (s.depends, cs);
  }
  for (scalar s in{p, pf})
    s.embed_gradient = pressure_embed_gradient;
#endif // EMBED
#endif // TREE

  foreach ()
    foreach_dimension () dimensional (u.x[] == Delta / t);

  foreach ()
    foreach_dimension () ibmf.x[] = 0.;
}

event default_display (i = 0) display ("squares (color = 'u.x', spread = -1);");

double dtmax;


event init (i = 0) {

  trash ({uf});
  foreach_face () uf.x[] = fm.x[] * face_value (u.x, 0);

#if _MPI
  ibmeshmanager_update_pid ();
  ibmeshmanager_boundary ();
#endif

  event ("properties");

  dtmax = DT;
  event ("stability");

}

event set_dtmax (i++, last) dtmax = DT;

event stability (i++, last) {
  dt = dtnext (stokes ? dtmax : timestep (uf, dtmax));
}

event vof (i++, last);
event tracer_advection (i++, last);
event tracer_diffusion (i++, last);

event properties (i++, last) {
}

void prediction () {
  vector du;
  foreach_dimension () {
    scalar s = new scalar;
    du.x = s;
  }

  if (u.x.gradient)
    foreach ()
      foreach_dimension () {
#if EMBED
        if (!fs.x[] || !fs.x[1])
          du.x[] = 0.;
        else
#endif
          du.x[] = u.x.gradient (u.x[-1], u.x[], u.x[1]) / Delta;
      }
  else
    foreach ()
      foreach_dimension () {
#if EMBED
        if (!fs.x[] || !fs.x[1])
          du.x[] = 0.;
        else
#endif
          du.x[] = (u.x[1] - u.x[-1]) / (2. * Delta);
      }

  trash ({uf});
  foreach_face () {
    double un = dt * (u.x[] + u.x[-1]) / (2. * Delta), s = sign (un);
    int i = -(s + 1.) / 2.;
    uf.x[] =
      u.x[i] + (g.x[] + g.x[-1]) * dt / 4. + s * (1. - s * un) * du.x[i] * Delta / 2.;
#if dimension > 1
    if (fm.y[i, 0] && fm.y[i, 1]) {
      double fyy = u.y[i] < 0. ? u.x[i, 1] - u.x[i] : u.x[i] - u.x[i, -1];
      uf.x[] -= dt * u.y[i] * fyy / (2. * Delta);
    }
#endif
#if dimension > 2
    if (fm.z[i, 0, 0] && fm.z[i, 0, 1]) {
      double fzz = u.z[i] < 0. ? u.x[i, 0, 1] - u.x[i] : u.x[i] - u.x[i, 0, -1];
      uf.x[] -= dt * u.z[i] * fzz / (2. * Delta);
    }
#endif
    uf.x[] *= fm.x[];
  }

  delete ((scalar*) {du});
}

event advection_term (i++, last) {
  if (!stokes) {
    prediction ();
    mgpf = project (uf, pf, alpha, dt / 2., mgpf.nrelax);
    advection ((scalar*) {u}, uf, dt, (scalar*) {g});
  }
}

static void correction (double dt) {
  foreach ()
    foreach_dimension () u.x[] += dt * g.x[];
}

event viscous_term (i++, last) {
  if (constant (mu.x) != 0.) {
    correction (dt);
    mgu = viscosity (u, mu, rho, dt, mgu.nrelax);
    correction (-dt);
  }

  /**
  Match the centered solver ordering: reset the acceleration field before
  force events add their contributions. The projection step must still see the
  acceleration that was used to build uf. */

  if (!is_constant (a.x)) {
    face vector af = a;
    trash ({af});
    foreach_face () af.x[] = 0.;
  }
}

#include "ibm/IBKernels.h"

event interface_force_velocity_coupled (i++, last) {
  trash ({ibmf});
  foreach () {
    foreach_dimension () {
      ibmf.x[] = 0.;
    }
  }

  foreach_dimension () u.x.dirty = true;
  boundary ((scalar*) {u});

  foreach_ibnode_per_ibmesh () {
    if (mesh->model.type != IB_MODEL_VELOCITY_COUPLED)
      continue;

    foreach_dimension () {
      ibval (nvel.x) = 0.;
    }

    peskin_cosine_kernel_gather_dimensionless (node) {
      foreach_dimension () {
        ibval (nvel.x) += weight * u.x[];
      }
    }
  }

#if _MPI
  {
    IBscalar* list = NULL;
    foreach_dimension () {
      list = iblist_add (list, nvel.x);
    }
    ibmeshmanager_boundary (list);
  }
#endif

  ibmeshmanager_evaluate_velocity_coupled_midpoints (dt);

  foreach_ibnode_per_ibmesh () {
    if (mesh->model.type != IB_MODEL_VELOCITY_COUPLED)
      continue;

    peskin_cosine_kernel_spread_dimensionless (node) {
      foreach_dimension () {
        ibmf.x[] += -weight * ibval (nforce.x) * ibval (nweight) / dv ();
      }
    }
  }

  foreach_dimension () ibmf.x.dirty = true;
  boundary ((scalar*) {ibmf});

  face vector af = a;
  foreach_face ()
    if (fm.x[] > 1.e-20)
      af.x[] += 0.5 * alpha.x[] * (ibmf.x[] + ibmf.x[-1]);
}

event acceleration (i++, last) {

  trash ({uf});
  foreach_face () uf.x[] = fm.x[] * (face_value (u.x, 0) + dt * a.x[]);
}

void centered_gradient (scalar p, vector g) {

  face vector gf[];
  foreach_face () gf.x[] = fm.x[] * a.x[] - alpha.x[] * (p[] - p[-1]) / Delta;

  trash ({g});
  foreach ()
    foreach_dimension () g.x[] = (gf.x[] + gf.x[1]) / (fm.x[] + fm.x[1] + SEPS);
}

event projection (i++, last) {
  mgp = project (uf, p, alpha, dt, mgp.nrelax);
  centered_gradient (p, g);

  correction (dt);
}

event advance_velocity_coupled_interface (i++, last) {
  foreach_dimension () u.x.dirty = true;
  boundary ((scalar*) {u});

  foreach_ibnode_per_ibmesh () {
    if (mesh->model.type != IB_MODEL_VELOCITY_COUPLED)
      continue;

    foreach_dimension () {
      ibval (nvel.x) = 0.;
    }

    peskin_cosine_kernel_gather_dimensionless (node) {
      foreach_dimension () {
        ibval (nvel.x) += weight * u.x[];
      }
    }
  }

#if _MPI
  {
    IBscalar* list = NULL;
    foreach_dimension () {
      list = iblist_add (list, nvel.x);
    }
    ibmeshmanager_boundary (list);
  }
#endif

#if ADVECT_LAG_RK2
  ibmeshmanager_evaluate_velocity_coupled_midpoints (dt);

  foreach_dimension () u.x.dirty = true;
  boundary ((scalar*) {u});

  foreach_ibnode_per_ibmesh () {
    if (mesh->model.type != IB_MODEL_VELOCITY_COUPLED)
      continue;

    foreach_dimension () {
      ibval (nvel.x) = 0.;
    }

    peskin_cosine_kernel_gather_dimensionless (node) {
      foreach_dimension () {
        ibval (nvel.x) += weight * u.x[];
      }
    }
  }

#if _MPI
  {
    IBscalar* list = NULL;
    foreach_dimension () {
      list = iblist_add (list, nvel.x);
    }
    ibmeshmanager_boundary (list);
  }
#endif
#endif

  ibmeshmanager_advance_velocity_coupled_positions (dt);
}

event end_timestep (i++, last);

#if TREE
event adapt (i++, last) {
#if _MPI
  ibmeshmanager_update_pid ();
#endif
#if EMBED
  fractions_cleanup (cs, fs);
  foreach_face () if (uf.x[] && !fs.x[]) uf.x[] = 0.;
#endif
  event ("properties");
}
#endif
