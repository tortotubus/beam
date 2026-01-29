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

// vector forcing[];
// scalar stencils[];

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
  // ib_mesh_manager_advance_lagrangian_mesh ();

  // 1a. Update node ownership
  ibnode_list_set_pid ();

  // 2. Generate the stencil cache
  ibnode_list_generate_stencil_cache ();

  // 3. Interpolate \f(u*\f) with \f(J u*\f)
  ibnode_list_interpolate_eulerian_velocities ();

  // 4. Build the right-hand side vector \f(b = v_{\rm lag} - J u^* \f)
  ibnode_list_compute_constraint_rhs ();

  // 5. Solve A \lambda = b via CG/Uzawa on nodes
  ibnode_list_solve_lambda_CG (dt);

  // 6. Spread lambda to grid forcing \f(J^T \lambda\f)
  ibnode_list_spread_eulerian_forcing ();

  // 7. Add forcing to the acceleration field
  foreach_face () {
    if (fm.x[] > 1e-20) {
      ae.x[] += .5 * alpha.x[] * (forcing.x[] + forcing.x[-1]);
    }
  }

  ibnode_list_tag_stencil (stencils);
}

/** At the end of the simulation, we free the allocated memory.*/
event cleanup (t = end) {
  ib_mesh_manager_free ();
  ibnode_list_free ();
}
 