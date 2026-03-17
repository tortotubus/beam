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

scalar pi[];
vector ui[], gi[];
scalar pfi[];
face vector ufi[];


IBvector gravity;

IBscalar nweight;
IBvector eulvel;
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
mgstats mgp = {0}, mgpf = {0}, mgu = {0};

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
  mgu = (mgstats) {0};

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

event properties (i++, last);

void prediction (vector uu, vector gg, double dtt) {
  vector du;
  foreach_dimension () {
    scalar s = new scalar;
    du.x = s;
  }

  if (uu.x.gradient)
    foreach ()
      foreach_dimension () {
#if EMBED
        if (!fs.x[] || !fs.x[1])
          du.x[] = 0.;
        else
#endif
          du.x[] = u.x.gradient (uu.x[-1], uu.x[], uu.x[1]) / Delta;
      }
  else
    foreach ()
      foreach_dimension () {
#if EMBED
        if (!fs.x[] || !fs.x[1])
          du.x[] = 0.;
        else
#endif
          du.x[] = (uux[1] - uu.x[-1]) / (2. * Delta);
      }

  trash ({uf});
  foreach_face () {
    double un = dtt * (uu.x[] + uu.x[-1]) / (2. * Delta), s = sign (un);
    int i = -(s + 1.) / 2.;
    uf.x[] = uu.x[i] + (gg.x[] + gg.x[-1]) * dtt / 4. +
             s * (1. - s * un) * duu.x[i] * Delta / 2.;
#if dimension > 1
    if (fm.y[i, 0] && fm.y[i, 1]) {
      double fyy = uu.y[i] < 0. ? uu.x[i, 1] - uu.x[i] : uu.x[i] - uu.x[i, -1];
      uf.x[] -= dtt * uu.y[i] * fyy / (2. * Delta);
    }
#endif
#if dimension > 2
    if (fm.z[i, 0, 0] && fm.z[i, 0, 1]) {
      double fzz = uu.z[i] < 0. ? uu.x[i, 0, 1] - uu.x[i] : uu.x[i] - uu.x[i, 0, -1];
      uf.x[] -= dtt * uu.z[i] * fzz / (2. * Delta);
    }
#endif
    uf.x[] *= fm.x[];
  }

  delete ((scalar*) {du});
}

event copy(i++, last) {
  foreach() {
    pi[] = p[];
    pfi[] = pf[];
    foreach_dimension() {
      ui.x[] = u.x[];
      gi.x[] = g.x[];
    }
  }
  foreach_face() {
    ufi.x[] = uf.x[];
  }
}

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

event advance_lagrangian_mesh (i++, last) {
  foreach_ibnode() 
    foreach_dimension()
      node->f.x += ibval(gravity.x);

  ibmeshmanager_advance_positions (dt);

  foreach_ibnode() 
    foreach_dimension()
      node->f.x -= ibval(gravity.x);
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
