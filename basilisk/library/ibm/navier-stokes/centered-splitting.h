/**
 * Switch from custom Helmholtz (poisson) diffusion to Basilisk's viscosity.h
 *
 * Key idea:
 *   viscosity(u, mu, rho, dt, ...) advances u by an *implicit* diffusion step:
 *       rho (u^{n+1} - u^n)/dt = ∇·(2 mu D(u^{n+1}))   (conceptually)
 *   i.e. it solves a variable-coefficient vector Helmholtz-like system using
 *   Basilisk's operators/metrics and (if you include viscosity-embed.h)
 * embedded boundary support.
 *
 * To match your split scheme:
 *   - predictor diffusion: apply viscosity() to ustar with dt_eff =
 * alpha_split*dt
 *   - corrector diffusion: apply viscosity() to ustarstar with dt_eff =
 * beta_split*dt
 *
 * IMPORTANT:
 *   viscosity() only does diffusion. Your forcing f still enters explicitly in
 *   the corrector RHS (as a "kick") before the diffusion solve, unless you
 *   redesign to treat forcing semi-implicitly.
 */

#include "run.h"
#include "timestep.h"
#include "bcg.h" 
#include "viscosity.h" 

/* Primary unknowns */
vector u[];
scalar p[];

/* Working fields */
vector ustar[];
vector ustarstar[];
vector f[]; // cell-centered force density (force/volume)
vector g[];
face vector uf[];
scalar pf[];

/* Basilisk-style material properties */
(const) face vector mu = zerof; 
(const) face vector alpha = unityf; 
(const) face vector a = zerof; 
(const) scalar rho = unity;

/* Split controls */
double alpha_split = 0.5;
double beta_split = 0.5;

/* Multigrid stats */
mgstats mgp = {0};
mgstats mgpf = {0};
mgstats mgu = {0}; // viscosity solver stats (like centered.h)

/**
## Boundary conditions

For the default symmetric boundary conditions, we need to ensure that
the normal component of the velocity is zero after projection. This
means that, at the boundary, the acceleration $\mathbf{a}$ must be
balanced by the pressure gradient. Taking care of boundary orientation
and staggering of $\mathbf{a}$, this can be written */

#define neumann_pressure(i) (a.n[i] * fm.n[i] / alpha.n[i])
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

static void build_uf_from_u (vector uc, face vector uf) {
  trash ({uf});
  foreach_face () { uf.x[] = fm.x[] * face_value (uc.x, 0); }
  boundary({uf});
}


void centered_gradient (scalar p, vector g)
{

  /**
  We first compute a face field $\mathbf{g}_f$ combining both
  acceleration and pressure gradient. */

  face vector gf[];
  foreach_face()
    gf.x[] = fm.x[]*a.x[] - alpha.x[]*(p[] - p[-1])/Delta;

  /**
  We average these face values to obtain the centered, combined
  acceleration and pressure gradient field. */

  trash ({g});
  foreach()
    foreach_dimension()
      g.x[] = (gf.x[] + gf.x[1])/(fm.x[] + fm.x[1] + SEPS);
}

event defaults (i = 0)
{

  /**
  We reset the multigrid parameters to their default values. */
  
  mgp = (mgstats){0};
  mgpf = (mgstats){0};
  mgu = (mgstats){0};  
  
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
  }
  else if (!is_constant(alpha.x)) {
    face vector alphav = alpha;
    foreach_face()
      alphav.x[] = fm.x[];
  }

  /**
  On trees, refinement of the face-centered velocity field needs to
  preserve the divergence-free condition. */

#if TREE
  uf.x.refine = refine_face_solenoidal;
#endif // TREE

  /**
  We set the dimensions of the velocity field. */

  foreach()
    foreach_dimension()
      dimensional (u.x[] == Delta/t);
}


double dtmax;

event init (i = 0) {

  trash ({uf});
  foreach_face()
    uf.x[] = fm.x[]*face_value (u.x, 0);
  
  /* Let user fill alpha/rho/mu */
  event ("properties"); 
  centered_gradient (p, g);

  dtmax = DT;
  event ("stability");
}

event set_dtmax (i++, last) { dtmax = DT; }

event face_velocity (i++, last) { build_uf_from_u (u, uf); }

event stability (i++, last) {
  dt = dtnext (timestep (uf, dtmax));
}

/* User hook: define alpha, rho, mu (and possibly reset them after adapt) */
event properties (i++, last);

/**
 * (1) Predictor: u*
 *   - copy u^n
 *   - build and (optionally) half-step-project advecting uf
 *   - explicit advection update of ustar
 *   - optional implicit diffusion using viscosity() with dt_eff =
 * alpha_split*dt
 */
event prediction (i++,last) {
  foreach ()
    foreach_dimension () ustar.x[] = u.x[];

  // build_uf_from_u (u, uf);

  /* Half-step projection for advecting velocity (centered.h style) */
  mgpf = project (uf, pf, alpha, dt / 2., mgpf.nrelax);

  /* Explicit advection: advances ustar using uf */
  advection ((scalar*) {ustar}, uf, dt, (scalar*) {g});

  /* Alpha-split diffusion: implicit diffusion over dt_eff */
  if (alpha_split > 0.) {
    /* viscosity() modifies its first argument in-place.
       It assumes RHS is the current field; so this is:
         ustar <- ustar + dt_eff/rho * ∇·(2 mu D(ustar_{new}))  (implicit)
       which corresponds to solving a Helmholtz-like diffusion step. */
    mgu = viscosity (ustar, mu, rho, alpha_split * dt, mgu.nrelax);
  }

  boundary ({ustar});
}

/* IBM hooks */
event compute_interface_force (i++);
event advance_interface_kinematics (i++);

/* Default: no forcing */
event spread_interface_force (i++) {
  foreach ()
    foreach_dimension () f.x[] = 0.;
  boundary ({f});
}

/**
 * (5) Corrector: u**
 *
 * Whiteboard says:
 *    (rho/dt - beta*mu Δ) u** = (rho/dt) u* + f
 *
 * With viscosity(), easiest consistent equivalent is:
 *   A) explicit "kick" by forcing: ustarstar = ustar + dt * f/rho
 *   B) implicit diffusion over dt_eff = beta_split*dt: viscosity(ustarstar,...)
 *
 * This matches the idea of splitting forcing and diffusion, while using
 * Basilisk's diffusion operator.
 */
event correction (i++,last) {
  /* Start from u* */
  foreach ()
    foreach_dimension () ustarstar.x[] = ustar.x[];

  /* Explicit forcing kick: u <- u + dt*(f/rho) */
  foreach ()
    foreach_dimension () ustarstar.x[] += dt * (f.x[] / (rho[] + SEPS));

  /* Beta-split diffusion */
  if (beta_split > 0.) {
    mgu = viscosity (ustarstar, mu, rho, beta_split * dt, mgu.nrelax);
  }

  boundary ({ustarstar});
}

/**
 * (6) Projection to enforce incompressibility and update u^{n+1}
 */
event projection (i++,last) {
  build_uf_from_u (ustarstar, uf);

  mgp = project (uf, p, alpha, dt, mgp.nrelax);

  /* Centered pressure correction */
  foreach ()
    foreach_dimension () {
      double dp = (p[1] - p[-1]) / (2. * Delta);
      /* Use centered 1/rho; alpha is face 1/rho used in project(). */
      u.x[] = ustarstar.x[] - dt * (1. / (rho[] + SEPS)) * dp;
    }

  centered_gradient (p, g);   // refresh g for next timestep
  
  boundary ({u});
}

event end_timestep (i++, last);

#if TREE
event adapt (i++, last) {
  event ("properties");
}
#endif
