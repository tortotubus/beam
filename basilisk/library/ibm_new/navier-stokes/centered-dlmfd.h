
#include "run.h"
#include "timestep.h"
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
vector tmp_force[];
vector tmp_vel[];

IBvector gravity;
IBvector eulvel; 
IBvector rhs;
IBvector res;
IBvector w;
IBvector Ay;

(const) face vector mu = zerof, a = zerof, alpha = unityf;
(const) scalar rho = unity;
bool stokes = false;

face vector mu_alpha[], mu_beta[];
mgstats mgp = {0}, mgpf = {0}, mgu_a = {0}, mgu_b = {0};

double alpha_split = 0.5;
double beta_split = 0.5;

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
  new_ibvector (rhs);
  new_ibvector (res);
  new_ibvector (w);
  new_ibvector (Ay);

  ibmeshmanager_init (0);

  if (is_constant (a.x)) {
    a = new face vector;
    foreach_face () a.x[] = 0.;
  }

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
}

event default_display (i = 0) display ("squares (color = 'u.x', spread = -1);");

double dtmax;

event init_ib (i = 0) {
#if _MPI
  ibmeshmanager_update_pid ();
  ibmeshmanager_boundary ();
#endif

#if TREE
  adapt_wavelet_ibm (NULL, NULL, 10);
#endif
}

event init (i = 0) {
  trash ({uf});
  foreach_face () uf.x[] = fm.x[] * face_value (u.x, 0);

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

event advance_lagrangian_mesh (i++, last) {
  foreach_ibnode () foreach_dimension () node->f.x += ibval (gravity.x);

  ibmeshmanager_advance_positions (dt);
}

event interpolate_eulerian_velocities (i++, last) {
  foreach_ibnode () {
    foreach_dimension () {
      ibval (eulvel.x) = 0;
    }
    peskin_cosine_kernel_gather_dimensionless (node) {
      foreach_dimension () {
        ibval (eulvel.x) += weight * u.x[];
      }
    }
  }
}

event compute_constraint_rhs (i++, last) {
  foreach_ibnode () {
    foreach_dimension () {
      ibval (rhs.x) = ibval (eulvel.x) - node->vel.x;
    }
  }
}

void ib_matvec_Aw (double dt) {

  foreach () {
    foreach_dimension () {
      tmp_force.x[] = 0.;
      tmp_vel.x[] = 0.;
    }
  }

  // 1. g = J^T w on grid
  foreach_ibnode () {
    peskin_cosine_kernel_spread_dimensionless (node) {
      foreach_dimension () {
        tmp_force.x[] += (weight / dv ()) * ibval (w.x) * ibval (nweight);
      }
    }
  }

  // 2. c = M_f^{-1} g ~ g / rho_f
  foreach () {
    foreach_dimension () {
      tmp_vel.x[] = tmp_force.x[] / rho[]; // should be rho_f here
    }
  }

  // 3. d = J c, then scale by dt
  foreach_ibnode () {
    foreach_dimension () ibval (Ay.x) = 0.;

    peskin_cosine_kernel_gather_dimensionless (node) foreach_dimension () ibval (Ay.x) +=
      weight * tmp_vel.x[];

    foreach_dimension () ibval (Ay.x) *= dt;
  }
}

double cgtol = 1e-5;
int cgiter = 0;
int cgiter_max = 50;

void ib_solve_lambda_CG (double dt) {
  foreach_ibnode () {
    foreach_dimension () {
      node->f.x = 0.;
      ibval (res.x) = ibval (rhs.x);
      ibval (w.x) = ibval (res.x);
    }
  }

  double rr_old = 0.;

  foreach_ibnode () {
    foreach_dimension () {
      rr_old += sq (ibval (res.x)) * ibval (nweight);
    }
  }

  if (sqrt (rr_old) < cgtol)
    return;

  for (cgiter = 0; cgiter < cgiter_max; cgiter++) {
    ib_matvec_Aw (dt);
    double wy = 0.;

    foreach_ibnode () {
      foreach_dimension () {
        wy += ibval (w.x) * ibval (Ay.x) * ibval (nweight);
      }
    }

    if (fabs (wy) < 1e-30) {
      break;
    }

    double alpha = rr_old / wy;
    double rr_new = 0.;

    // \lambda_{k+1} = \lambda_{k} + \alpha w_{k}\f), and \f(r_{k+1} = r_{k}
    // - \alpha y_k
    foreach_ibnode () {
      foreach_dimension () {
        node->f.x += alpha * ibval (w.x);      // lambda update
        ibval (res.x) -= alpha * ibval (Ay.x); // residual update
        rr_new += sq (ibval (res.x)) * ibval (nweight);
      }
    }

    if (sqrt (rr_new) < cgtol) {
      break;
    }

    double beta = rr_new / rr_old;
    rr_old = rr_new;

    // w_{k+1} = r_{k+1} + \beta w_{k}
    foreach_ibnode () {
      foreach_dimension () {
        ibval (w.x) = ibval (res.x) + beta * ibval (w.x);
      }
    }
  }
}

event solve_lambda_CG (i++, last) {
  ib_solve_lambda_CG (dt);
}

event spread_eulerian_forcing (i++, last) {
  face vector ae = a;

  foreach ()
    foreach_dimension () ibmf.x[] = 0.;

  foreach_ibnode () peskin_cosine_kernel_spread_dimensionless (node) foreach_dimension ()
    ibmf.x[] -= weight / dv () * node->f.x * ibval (nweight);

  // foreach_face ()
  //   if (fm.x[] > 1e-20)
  //     ae.x[] += .5 * alpha.x[] * (ibmf.x[] + ibmf.x[-1]);

  foreach () {
    foreach_dimension () {
      u.x[] += dt * ibmf.x[];
    }
  }
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

event end_timestep (i++, last) {}

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
