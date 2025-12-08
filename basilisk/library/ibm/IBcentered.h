#include "navier-stokes/centered.h"

#include "library/ibm/IBMesh.h"
#include "library/ibm/IBMeshManager.h"

/* ====================================================================================================================
 * Centered Navier-Stokes Solver Events
 * ====================================================================================================================
 */

/**
 * @brief
 *
 * See
 * [centered](https://basilisk.fr/src/navier-stokes/centered.h#time-integration).
 */
event defaults (i = 0) {
  if (is_constant (a.x)) {
    a = new face vector;
    foreach_face () a.x[] = 0.;
  }
}

event acceleration (i++) {
  face vector ae = a;

  // 0. Clear the forcing term on the Eulerian grid
  foreach () {
    if (cm[] > 1e-20) {
      foreach_dimension () {
        forcing.x[] = 0.;
      }
    }
  }

  // 1. Advance lagrangian mesh
  ib_mesh_manager_advance_lagrangian_mesh ();

  // 2. Generate the stencil cache
  ib_mesh_manager_generate_stencil_cache ();

  // 3. Interpolate \f(u*\f) with \f(J u*\f)
  ib_mesh_manager_interpolate_eulerian_velocities ();

  // 4. Build the right-hand side vector \f(b = v_{\rm lag} - J u^* \f)
  ib_mesh_manager_compute_constraint_rhs ();

  // 5. Solve A \lambda = b via CG/Uzawa on nodes
  ib_mesh_manager_solve_lambda_CG (dt);

  // 6. Spread lambda to grid forcing \f(J^T \lambda\f)
  ib_mesh_manager_spread_eulerian_forcing (forcing);

  // 7. Add forcing to the acceleration field
  foreach_face () {
    if (fm.x[] > 1e-20) {
      ae.x[] += .5 * alpha.x[] * (forcing.x[] + forcing.x[-1]);
    }
  }
}

event acceleration(i++) {
  ib_mesh_manager_sort_triplets ();
}

/** At the end of the simulation, we free the allocated memory.*/
event cleanup (t = end) {
  ib_mesh_manager_free ();
}
 