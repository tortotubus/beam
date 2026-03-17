/**
# Incompressible Navier--Stokes solver (centered formulation)

We wish to approximate numerically the incompressible,
variable-density Navier--Stokes equations
$$
\partial_t\mathbf{u}+\nabla\cdot(\mathbf{u}\otimes\mathbf{u}) =
\frac{1}{\rho}\left[-\nabla p + \nabla\cdot(2\mu\mathbf{D})\right] +
\mathbf{a}
$$
$$
\nabla\cdot\mathbf{u} = 0
$$
with the deformation tensor
$\mathbf{D}=[\nabla\mathbf{u} + (\nabla\mathbf{u})^T]/2$.

The scheme implemented here is close to that used in Gerris ([Popinet,
2003](/src/references.bib#popinet2003), [Popinet,
2009](/src/references.bib#popinet2009), [Lagrée et al,
2011](/src/references.bib#lagree2011)).

We will use the generic time loop, a CFL-limited timestep, the
Bell-Collela-Glaz advection scheme and the implicit viscosity
solver. If embedded boundaries are used, a different scheme is used
for viscosity. */

#include "run.h"
#include "timestep.h"
#include "bcg.h"
#if EMBED
#include "viscosity-embed.h"
#else
#include "viscosity.h"
#endif

#include "library/ibm/IBMeshManager.h"
#if TREE
#include "library/ibm/IBAdapt.h"
#endif

/**
The primary variables are the centered pressure field $p$ and the
centered velocity field $\mathbf{u}$. The centered vector field
$\mathbf{g}$ will contain pressure gradients and acceleration terms.

We will also need an auxilliary face velocity field $\mathbf{u}_f$ and
the associated centered pressure field $p_f$. */

scalar p[];
vector u[], g[];
scalar pf[];
face vector uf[];

vector ibmf[];
vector tmp_force[];
vector tmp_vel[];

IBvector gravity;

IBscalar nweight;
IBvector eulvel;
// IBvector lagvel;
// IBvector force;
IBvector rhs;
IBvector res;
IBvector w;
IBvector Ay;

/**
In the case of variable density, the user will need to define both the
face and centered specific volume fields ($\alpha$ and $\alpha_c$
respectively) i.e. $1/\rho$. If not specified by the user, these
fields are set to one i.e. the density is unity.

Viscosity is set by defining the face dynamic viscosity $\mu$; default
is zero.

The face field $\mathbf{a}$ defines the acceleration term; default is
zero.

The statistics for the (multigrid) solution of the pressure Poisson
problems and implicit viscosity are stored in *mgp*, *mgpf*, *mgu*
respectively.

If *stokes* is set to *true*, the velocity advection term
$\nabla\cdot(\mathbf{u}\otimes\mathbf{u})$ is omitted. This is a
reference to [Stokes flows](http://en.wikipedia.org/wiki/Stokes_flow)
for which inertia is negligible compared to viscosity. */

(const) face vector mu = zerof, a = zerof, alpha = unityf;
(const) scalar rho = unity;
bool stokes = false;

face vector mu_alpha[], mu_beta[];
mgstats mgp = {0}, mgpf = {0}, mgu_a = {0}, mgu_b = {0};

double alpha_split = 0.5;
double beta_split = 0.5;

/**
## Boundary conditions

For the default symmetric boundary conditions, we need to ensure that
the normal component of the velocity is zero after projection. This
means that, at the boundary, the acceleration $\mathbf{a}$ must be
balanced by the pressure gradient. Taking care of boundary orientation
and staggering of $\mathbf{a}$, this can be written */

#if EMBED
#define neumann_pressure(i)                                                    \
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

/**
For [embedded boundaries on trees](/src/embed-tree.h), we need to
define the pressure gradient for prolongation of pressure close to
embedded boundaries. */

#if TREE && EMBED
void pressure_embed_gradient (Point point, scalar p, coord* g) {
  foreach_dimension () g->x = rho[] / (cm[] + SEPS) * (a.x[] + a.x[1]) / 2.;
}
#endif // TREE && EMBED

/**
## Initial conditions */

event defaults (i = 0) {
  new_ibscalar (nweight);
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

  /**
  We reset the multigrid parameters to their default values. */

  mgp = (mgstats) {0};
  mgpf = (mgstats) {0};
  mgu_a = (mgstats) {0};
  mgu_b = (mgstats) {0};

  CFL = 0.8;

  /**
  The pressures are never dumped. */

  p.nodump = pf.nodump = true;

  /**
  The default density field is set to unity (times the metric and the
  solid factors). */

  if (alpha.x.i == unityf.x.i) {
    alpha = fm;
    rho = cm;
  } else if (!is_constant (alpha.x)) {
    face vector alphav = alpha;
    foreach_face () alphav.x[] = fm.x[];
  }

  /**
  On trees, refinement of the face-centered velocity field needs to
  preserve the divergence-free condition. */

#if TREE
  uf.x.refine = refine_face_solenoidal;

  /**
  When using [embedded boundaries](/src/embed.h), the restriction and
  prolongation operators need to take the boundary into account. */

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

  /**
  We set the dimensions of the velocity field. */

  foreach ()
    foreach_dimension () dimensional (u.x[] == Delta / t);
}

/**
We had some objects to display by default. */

event default_display (i = 0) display ("squares (color = 'u.x', spread = -1);");

/**
After user initialisation, we initialise the face velocity and fluid
properties. */

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

  /**
  We update fluid properties. */

  event ("properties");

  /**
  We set the initial timestep (this is useful only when restoring from
  a previous run). */

  dtmax = DT;
  event ("stability");
}

/**
## Time integration

The timestep for this iteration is controlled by the CFL condition,
applied to the face centered velocity field $\mathbf{u}_f$; and the
timing of upcoming events. */

event set_dtmax (i++, last) dtmax = DT;

event stability (i++, last) {
  dt = dtnext (stokes ? dtmax : timestep (uf, dtmax));
}

/**
If we are using VOF or diffuse tracers, we need to advance them (to
time $t+\Delta t/2$) here. Note that this assumes that tracer fields
are defined at time $t-\Delta t/2$ i.e. are lagging the
velocity/pressure fields by half a timestep. */

event vof (i++, last);
event tracer_advection (i++, last);
event tracer_diffusion (i++, last);

/**
The fluid properties such as specific volume (fields $\alpha$ and
$\alpha_c$) or dynamic viscosity (face field $\mu_f$) -- at time
$t+\Delta t/2$ -- can be defined by overloading this event. */

event properties (i++, last) {
  if (!is_constant (mu.x)) {
    foreach_face () {
      mu_alpha.x[] = mu.x[] * alpha_split;
      mu_beta.x[] = mu.x[] * beta_split;
    }
  }
}

/**
### Predicted face velocity field

For second-order in time integration of the velocity advection term
$\nabla\cdot(\mathbf{u}\otimes\mathbf{u})$, we need to define the face
velocity field $\mathbf{u}_f$ at time $t+\Delta t/2$. We use a version
of the Bell-Collela-Glaz [advection scheme](/src/bcg.h) and the
pressure gradient and acceleration terms at time $t$ (stored in vector
$\mathbf{g}$). */

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
    uf.x[] = u.x[i] + (g.x[] + g.x[-1]) * dt / 4. +
             s * (1. - s * un) * du.x[i] * Delta / 2.;
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

/**
### Advection term

We predict the face velocity field $\mathbf{u}_f$ at time $t+\Delta
t/2$ then project it to make it divergence-free. We can then use it to
compute the velocity advection term, using the standard
Bell-Collela-Glaz advection scheme for each component of the velocity
field. */

event advection_term (i++, last) {
  if (!stokes) {
    prediction ();
    mgpf = project (uf, pf, alpha, dt / 2., mgpf.nrelax);
    advection ((scalar*) {u}, uf, dt, (scalar*) {g});
  }
}

/**
### Viscous term

We first define a function which adds the pressure gradient and
acceleration terms. */

static void correction (double dt) {
  foreach ()
    foreach_dimension () u.x[] += dt * g.x[];
}

/**
The viscous term is computed implicitly. We first add the pressure
gradient and acceleration terms, as computed at time $t$, then call
the implicit viscosity solver. We then remove the acceleration and
pressure gradient terms as they will be replaced by their values at
time $t+\Delta t$. */

event alpha_viscous_term (i++, last) {
  if (constant (mu.x) != 0. && alpha_split != 0.) {
    correction (dt);
    mgu_a = viscosity (u, mu_alpha, rho, dt, mgu_a.nrelax);
    correction (-dt);
  }
}

/**
### Acceleration term

The acceleration term $\mathbf{a}$ needs careful treatment as many
equilibrium solutions depend on exact balance between the acceleration
term and the pressure gradient: for example Laplace's balance for
surface tension or hydrostatic pressure in the presence of gravity.

To ensure a consistent discretisation, the acceleration term is
defined on faces as are pressure gradients and the centered combined
acceleration and pressure gradient term $\mathbf{g}$ is obtained by
averaging.

The (provisionary) face velocity field at time $t+\Delta t$ is
obtained by interpolation from the centered velocity field. The
acceleration term is added. */

#include "library/ibm/IBKernels.h"

event advance_lagrangian_mesh (i++, last) {
  foreach_ibnode() 
    foreach_dimension()
      node->f.x += ibval(gravity.x);

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
      ibval (rhs.x) =  ibval (eulvel.x) - node->vel.x;
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
        tmp_force.x[] += (weight / dv()) * ibval (w.x) * ibval(nweight);
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
    foreach_dimension () 
      ibval (Ay.x) = 0.;
    
    peskin_cosine_kernel_gather_dimensionless (node) 
      foreach_dimension()
        ibval (Ay.x) += weight * tmp_vel.x[];
    
    foreach_dimension () 
      ibval (Ay.x) *= dt;
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
      rr_old += sq (ibval (res.x)) * ibval(nweight);
    }
  }

  if (sqrt (rr_old) < cgtol)
    return;

  for (cgiter = 0; cgiter < cgiter_max; cgiter++) {
    ib_matvec_Aw (dt);
    double wy = 0.;

    foreach_ibnode () {
      foreach_dimension () {
        wy += ibval (w.x) * ibval (Ay.x) * ibval(nweight);
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
        node->f.x += alpha * ibval (w.x); // lambda update
        ibval (res.x) -= alpha * ibval (Ay.x);  // residual update
        rr_new += sq (ibval (res.x)) * ibval(nweight);
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

  foreach()
    foreach_dimension()
      ibmf.x[] = 0.;
    
  foreach_ibnode () 
    peskin_cosine_kernel_spread_dimensionless (node) 
      foreach_dimension()
        ibmf.x[] -= weight / dv() * node->f.x * ibval(nweight);
    
  // foreach_face () 
  //   if (fm.x[] > 1e-20) 
  //     ae.x[] += .5 * alpha.x[] * (ibmf.x[] + ibmf.x[-1]);  

  foreach() {
    foreach_dimension() {
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

/**
## Approximate projection

This function constructs the centered pressure gradient and
acceleration field *g* using the face-centered acceleration field *a*
and the cell-centered pressure field *p*. */

void centered_gradient (scalar p, vector g) {

  /**
  We first compute a face field $\mathbf{g}_f$ combining both
  acceleration and pressure gradient. */

  face vector gf[];
  foreach_face () gf.x[] = fm.x[] * a.x[] - alpha.x[] * (p[] - p[-1]) / Delta;

  /**
  We average these face values to obtain the centered, combined
  acceleration and pressure gradient field. */

  trash ({g});
  foreach ()
    foreach_dimension () g.x[] = (gf.x[] + gf.x[1]) / (fm.x[] + fm.x[1] + SEPS);
}

/**
To get the pressure field at time $t + \Delta t$ we project the face
velocity field (which will also be used for tracer advection at the
next timestep). Then compute the centered gradient field *g*. */

event projection (i++, last) {
  mgp = project (uf, p, alpha, dt, mgp.nrelax);
  centered_gradient (p, g);

  /**
  We add the gradient field *g* to the centered velocity field. */

  correction (dt);
}

/**
Some derived solvers need to hook themselves at the end of the
timestep. */

event end_timestep (i++, last) {}

/**
## Adaptivity

After mesh adaptation fluid properties need to be updated. When using
[embedded boundaries](/src/embed.h) the fluid fractions and face
fluxes need to be checked for inconsistencies. */

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

/**
## See also

* [Double projection](double-projection.h)
* [Performance monitoring](perfs.h)
*/
