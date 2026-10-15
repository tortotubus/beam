#include "run.h"
#include "library/timestep.h"
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

IBvector gravity;
IBvector eulvel;
IBscalar eulrho;
IBscalar sumw2;
IBvector dforce;

(const) face vector mu = zerof, a = zerof, alpha = unityf;
(const) scalar rho = unity;
bool stokes = false;

face vector mu_alpha[], mu_beta[];
mgstats mgp = {0}, mgpf = {0}, mgu_a = {0}, mgu_b = {0};

double alpha_split = 0.5;
double beta_split = 0.5;
double ib_force_relaxation = 1.0;
int ib_richardson_iters = 2;

#if EMBED
#define neumann_pressure(i)                                                              \
  (alpha.n[i] ? a.n[i] * fm.n[i] / alpha.n[i] : a.n[i] * rho[] / (cm[] + SEPS))
#else
#define neumann_pressure(i) (a.n[i] * fm.n[i] / alpha.n[i])
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
  new_ibvector (eulvel);
  new_ibvector (gravity);
  new_ibvector (dforce);

  foreach_dimension () {
    ibnodump (eulvel.x) = true;
    ibnodump (gravity.x) = true;
    ibnodump (dforce.x) = true;
  }

  new_ibscalar (eulrho);
  new_ibscalar (sumw2);

  ibnodump (eulrho) = true;
  ibnodump (sumw2) = true;

  ibmeshmanager_init (0);

  mgp = (mgstats) {0};
  mgpf = (mgstats) {0};
  mgu_a = (mgstats) {0};
  mgu_b = (mgstats) {0};

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

// event init_ib (i = 0) {

// }

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

// #if TREE 
//   adapt_wavelet_ibm(NULL,NULL,0,1,all,true);
// #endif
}

event set_dtmax (i++, last) dtmax = DT;

event stability (i++, last) {
  dt = dtnext (stokes ? dtmax : timestep (uf, dtmax));
}

event vof (i++, last);
event tracer_advection (i++, last);
event tracer_diffusion (i++, last);

event properties (i++, last) {
  if (!is_constant (mu.x)) {
    foreach_face () {
      mu_alpha.x[] = mu.x[] * alpha_split;
      mu_beta.x[] = mu.x[] * beta_split;
    }
  }
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

event alpha_viscous_term (i++, last) {
  if (constant (mu.x) != 0. && alpha_split != 0.) {
    correction (dt);
    mgu_a = viscosity (u, mu_alpha, rho, dt, mgu_a.nrelax);
    correction (-dt);
  }
}

#include "ibm/IBKernels.h"

event advance_interface (i++) {
  foreach_ibnode () {
    foreach_dimension () ibval (nforce.x) += ibval (gravity.x);
  }
  ibmeshmanager_advance_positions (dt);
}

event interface_force_richardson (i++, last) {
  vector dibmf[];
  vector uib[];

  trash ({ibmf});

  // Working Eulerian velocity and current Eulerian IBM increment
  foreach () {
    foreach_dimension () {
      uib.x[] = u.x[];
      ibmf.x[] = 0.;  // accumulated Eulerian IBM force density
      dibmf.x[] = 0.; // current Richardson Eulerian force-density increment
    }
  }

  // node->f is the accumulated Lagrangian force density
  foreach_ibnode () {
    foreach_dimension () {
      ibval (nforce.x) = 0.;
    }
  }

  for (int k = 0; k < ib_richardson_iters; k++) {

    // Reset current Eulerian increment
    foreach () {
      foreach_dimension () {
        dibmf.x[] = 0.;
      }
    }

    // Refresh working field for gathers
    foreach_dimension () uib.x.dirty = true;
    boundary ((scalar*) {uib});

    // Gather current interpolated velocity and diagonal self-mobility proxy
    foreach_ibnode () {
      foreach_dimension () {
        ibval (eulvel.x) = 0.;
        ibval (dforce.x) = 0.; // current force-density increment
      }
      double lsumw2 = 0.;
      peskin_cosine_kernel_gather_dimensionless (node) {
        lsumw2 += sq (weight);
        foreach_dimension () {
          ibval (eulvel.x) += weight * uib.x[];
        }
      }
      ibval (sumw2) = lsumw2;
      ibval (eulrho) = 1.;
    }

#if _MPI
    {
      IBscalar* list = NULL;
      foreach_dimension () {
        list = iblist_add (list, eulvel.x);
      }
      list = iblist_add (list, sumw2);
      list = iblist_add (list, nweight);
      list = iblist_add (list, eulrho);
      ibmeshmanager_boundary (list);
    }
#endif

    // Compute Richardson increment in an approximate Lagrangian force density
    //
    // delta f_density = -omega * rho_L * slip / (dt * q * sumw2 / dV)
    //
    // where q = nodal quadrature measure and rho_L is a gathered fluid density
    // proxy near the marker. The extra dV factor keeps node->f in the same
    // Lagrangian density units as the DLMFD implementation.
    foreach_ibnode () {

#if dimension == 1
      double dV = (L0 / (1 << node->depth));
#elif dimension == 2
      double dV = sq (L0 / (1 << node->depth));
#else
      double dV = cube (L0 / (1 << node->depth));
#endif
      foreach_dimension () {
        double slip = ibval (nvel.x) - ibval (eulvel.x);
        // Current Richardson increment: force density on the object
        ibval (dforce.x) = -ib_force_relaxation * (ibval (eulrho) * dV * slip) /
                           (dt * ibval (nweight) * ibval (sumw2) + SEPS);
        // Accumulate total force density
        ibval (nforce.x) += ibval (dforce.x);
      }
    }

#if _MPI
    {
      IBscalar* list = NULL;
      foreach_dimension () list = iblist_add (list, dforce.x);
      ibmeshmanager_boundary (list);
    }
#endif

    // Spread force density -> Eulerian correction using nodal measure q
    //
    // nodal force = force_density * q, then divide by cell measure to recover
    // the Eulerian force density used in the fluid update.
    foreach_ibnode () {
      peskin_cosine_kernel_spread_dimensionless (node) {
        foreach_dimension () {
          dibmf.x[] += -weight * ibval (dforce.x) * ibval (nweight) / dv ();
        }
      }
    }

    foreach_dimension () dibmf.x.dirty = true;
    boundary ((scalar*) {dibmf});

    // Apply current Eulerian increment and accumulate total Eulerian forcing
    foreach () {
      foreach_dimension () {
        uib.x[] += dt * dibmf.x[] / (rho[] + SEPS);
        ibmf.x[] += dibmf.x[];
      }
    }
  }

  // Commit corrected velocity back to main field
  foreach () {
    foreach_dimension () {
      u.x[] = uib.x[];
    }
  }

  foreach_dimension () u.x.dirty = true;
  boundary ((scalar*) {u});
}

event beta_viscous_term (i++, last) {
  if (constant (mu.x) != 0. && beta_split != 0.) {
    correction (dt);
    mgu_b = viscosity (u, mu_beta, rho, dt, mgu_b.nrelax);
    correction (-dt);
  }
}

event acceleration (i++, last) {

  trash ({uf});
  foreach_face () uf.x[] = fm.x[] * (face_value (u.x, 0) + dt * a.x[]);

  /**
  We reset the acceleration field (if it is not a constant). */

  if (!is_constant (a.x)) {
    face vector af = a;
    trash ({af});
    foreach_face () af.x[] = 0.;
  }
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

  foreach_ibnode () {
    foreach_dimension () {
      ibval (eulvel.x) = 0;
    }
    double lsumw2 = 0.;
    peskin_cosine_kernel_gather_dimensionless (node) {
      lsumw2 += weight * weight;
      foreach_dimension () {
        ibval (eulvel.x) += weight * u.x[];
      }
    }
    ibval (sumw2) = lsumw2;
  }
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
