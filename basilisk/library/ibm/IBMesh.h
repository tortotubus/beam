#include "navier-stokes/centered.h"

#include "interface/ibm/immersedboundary.h"
#include "library/ibm/IBNode.h"

#include <stdio.h>

/* ====================================================================================================================
 * IBMesh
 * ====================================================================================================================
 */

#define RHO_F 1

typedef struct {
  ib_force_structure_model_t model;
  IBNode* nodes;
  int nn;
  int pid;
  bool isactive;
} IBMesh;

/*
 * Macro defintions
 */

macro foreach_ibnode (IBMesh* mesh = mesh) {
  for (int node_index = 0; node_index < mesh->nn; node_index++) {
    IBNode* node = &mesh->nodes[node_index];
    NOT_UNUSED (node);

    // clang-format off
    {...}
    // clang-format on
  }
}

/*
 * Functions
 */

/*!
 * @brief Default constructor for an IBMesh
 *
 * @param mesh Mesh to initialize
 *
 * @memberof IBMesh
 */
void ib_mesh_init (IBMesh* mesh) {
  mesh->model = NULL;
  mesh->nodes = NULL;
  mesh->nn = 0;
  mesh->isactive = false;
}

/*!
 * @brief Destructor for the IBMesh struct
 * @memberof IBMesh
 */
void ib_mesh_free (IBMesh* mesh) {
  if (mesh->nodes) {
    foreach_ibnode (mesh) {
      ibnode_free (node);
    }
    free (mesh->nodes);
    mesh->nodes = NULL;
  }
  mesh->nn = 0;
  mesh->isactive = false;
}

/*!
 * @brief Construct the IBMesh with a number of IBNodes
 *
 * @param mesh The IBMesh to initialize
 * @param nn The number of nodes to initialize
 *
 * @memberof IBMesh
 */
void ib_mesh_init_nn (IBMesh* mesh, int nn) {
  // mesh->model_type = NONE;
  mesh->nodes = (IBNode*) calloc (nn, sizeof (IBNode));
  mesh->nn = nn;
  mesh->isactive = true;

  foreach_ibnode (mesh) {
    ibnode_init (node);
  }
}

/*!
 * @brief Set the model for an IBMesh
 *
 * @param mesh The IBMesh
 * @param model The pointer to the force-coupled structural model
 *
 * @memberof IBMesh
 */
void ib_mesh_model (IBMesh* mesh, ib_force_structure_model_t model) {
#if DEBUG
  printf ("ib_mesh_model\n");
#endif
  // Free the previous mesh, if any
  ib_mesh_free (mesh);

  // Set the model pointer
  mesh->model = model;

  // Get the number of nodes from the C++ model
  int nn = ib_force_structure_model_get_number_of_nodes (model);

  // Initialize the basilisk C mesh
  ib_mesh_init_nn (mesh, nn);

  // Get the current-time mesh from the C++ model
  ib_structure_mesh_t cpp_mesh = ib_force_structure_model_get_current (model);

  foreach_ibnode (mesh) {
    foreach_dimension () {
      node->lagpos.x = cpp_mesh.position[node_index].x;
      node->lagvel.x = cpp_mesh.velocity[node_index].x;
    }
  }

  // Free the C++ mesh
  ib_structure_mesh_free (&cpp_mesh);
}

/*!
 * @brief Helper function to compute the centroid of a given IBMesh
 *
 * @param mesh IBMesh to compute the centroid of
 *
 * @memberof IBMesh
 */
vertex_t ib_mesh_get_centroid (IBMesh* mesh) {
  vertex_t centroid = {.x = 0, .y = 0, .z = 0};

  foreach_ibnode (mesh) {
    foreach_dimension () {
      centroid.x += node->lagpos.x;
    }
  }

  foreach_dimension () {
    if (centroid.x != 0)
      centroid.x /= (double) mesh->nn;
  }

  return centroid;
}

/*!
 * @brief Generate the stencils for each node of a given mesh.
 *
 * @param mesh
 *
 * @memberof IBMeshManager
 */
trace void ib_mesh_generate_stencil_cache (IBMesh* mesh) {
#if DEBUG
  printf ("ib_mesh_generate_stencil_cache\n");
#endif
  foreach_ibnode (mesh) {
    node->stencil.n = 0;
    coord lagpos = {
      .x = node->lagpos.x, .y = node->lagpos.y, .z = node->lagpos.z};
    foreach_neighbor_of_coord (lagpos, PESKIN_SUPPORT_RADIUS) {
      if (point.level >= 0 && point.level != grid->maxdepth) {
        fprintf (
          stderr,
          "warning: lagrangian stencil not fully resolved at level %d of "
          "%d.\n",
          point.level,
          grid->maxdepth);
      }

      cache_append (&node->stencil, point, 0);

      // #if _MPI
      // #if dimension == 1
      //       if (point.level >= 0 && _k == 0) {
      //         node->pid = cell.pid;
      //       } else {
      //         node->pid = -1;
      //       }
      // #elif dimension == 2
      //       if (point.level >= 0 && _k == 0 && _l == 0) {
      //         node->pid = cell.pid;
      //       } else {
      //         node->pid = -1;
      //       }
      // #else // dimension == 3
      //       if (point.level >= 0 && _l == 0 && _m == 0 && _n == 0) {
      //         node->pid = cell.pid;
      //       } else {
      //         node->pid = -1;
      //       }
      // #endif
      // #endif
    }
  }
}

vector forcing[];
vector tmp_force[];
vector tmp_vel[];

/*!
 * @brief Spread the IBM force at the Lagrangian nodes back
 * onto the Eulerian vector field
 *
 * @param mesh The mesh for to gather forcing from
 * @param forcing Eulerian vector field to hold the result in
 *
 * @memberof IBMesh
 */
trace void ib_mesh_spread_eulerian_forcing (IBMesh* mesh, vector forcing) {
#if DEBUG
  printf ("ib_mesh_spread_eulerian_forcing\n");
#endif
#if _MPI
  synchronize ((scalar*) {forcing});
#endif
  foreach_ibnode (mesh) {
    foreach_cache (node->stencil) {
      coord lagpos = {
        .x = node->lagpos.x, .y = node->lagpos.y, .z = node->lagpos.z};
      peskin_kernel (lagpos) {
        foreach_dimension () {
          forcing.x[] += peskin_weight * node->force.x * node->weight;
#if dimension == 1
          forcing.x[] += (peskin_weight /      (Delta)) * node->force.x * node->weight;
#elif dimension == 2
          forcing.x[] += 1 * (peskin_weight /   sq (Delta)) * node->force.x * node->weight;
#else
          forcing.x[] += (peskin_weight / cube (Delta)) * node->force.x * node->weight;
#endif
        }
      }
    }
  }

#if _MPI
  synchronize ((scalar*) {forcing});
#endif
}

/*!
 * @brief Interpolate Eulerian velocities from the Eulerian mesh to the
 * Lagrangian nodes
 *
 * @param mesh The mesh for to gather forcing from
 *
 * @memberof IBMesh
 */
trace void ib_mesh_interpolate_eulerian_velocities (IBMesh* mesh) {
#if DEBUG
  printf ("ib_mesh_interpolate_eulerian_velocities\n");
#endif

#if _MPI
  synchronize ((scalar*) {u});
#endif
  foreach_ibnode (mesh) {
    foreach_dimension () {
      node->eulvel.x = 0.;
    }
    foreach_cache (node->stencil) {
      coord lagpos = {
        .x = node->lagpos.x, .y = node->lagpos.y, .z = node->lagpos.z};
      peskin_kernel (lagpos) {
        foreach_dimension () {
          node->eulvel.x += peskin_weight * u.x[];
        }
      }
    }
  }

#if _MPI
  synchronize ((scalar*) {forcing});
#endif
}

scalar stencils[];

/*!
 * @brief Injects the stencils field with noise so that we may call adapt_wavlet
 * to ensure sufficient cells on the same level near our IBNodes
 *
 * @param mesh The immersed boundary mesh
 *
 * @memberof IBMesh
 */
void ib_mesh_tag_stencil (IBMesh* mesh) {
#if DEBUG
  printf ("tag_stencil\n");
#endif

  foreach_ibnode (mesh) {
    foreach_cache (node->stencil) {
      if (point.level >= 0) {
        coord lagpos = {
          .x = node->lagpos.x, .y = node->lagpos.y, .z = node->lagpos.z};
        peskin_kernel (lagpos) {
          double length = 0;
          foreach_dimension () {
            length += dist.x;
          }
          length = sq (length);
#if dimension == 1
          stencils[] = length / (2. * Delta) * (2. + noise ());
#elif dimension == 2
          stencils[] = length / sq (2. * Delta) * (2. + noise ());
#else
          stencils[] = length / cube (2. * Delta) * (2. + noise ());
#endif
        }
      }
    }
  }
}

/*!
 * @brief Update the position and velocity of the lagrangian points according to
 * their kinematic law(s)
 *
 * @param mesh
 *
 * @memberof IBMesh
 */
trace void ib_mesh_advance_lagrangian_mesh (IBMesh* mesh) {
#if DEBUG
  printf ("ib_mesh_advance_lagrangian_mesh\n");
#endif

  vertex_t* force = calloc (mesh->nn, sizeof (vertex_t));

  foreach_ibnode (mesh) {
    foreach_dimension () {
      force[node_index].x = -node->force.x * node->weight + node->gravity.x;
    }
  }

  ib_structure_mesh_t cpp_mesh =
    ib_force_structure_model_get_next (mesh->model, force, mesh->nn, dt);

  free (force);

  foreach_ibnode (mesh) {
    foreach_dimension () {
      node->lagpos.x = cpp_mesh.position[node_index].x;
      node->lagvel.x = cpp_mesh.velocity[node_index].x;
    }
  }

  ib_structure_mesh_free (&cpp_mesh);
}

/*!
 * @brief Compute the kinematic no-slip residual
 *
 * @memberof IBMesh
 */
trace void ib_mesh_compute_constraint_rhs (IBMesh* mesh) {
  foreach_ibnode (mesh) {
    foreach_dimension () {
      node->rhs.x = node->lagvel.x - node->eulvel.x; // b = v_lag - J u*
    }
  }
}

/*!
 * @memberof IBMesh
 */
trace void ib_mesh_matvec_Aw (IBMesh* mesh, double dt) {

  foreach () {
    foreach_dimension () {
      tmp_force.x[] = 0;
      tmp_vel.x[] = 0;
    }
  }

  // 1. g = J^T w on grid
  foreach_ibnode (mesh) {
    foreach_cache (node->stencil) {
      coord lagpos = {
        .x = node->lagpos.x, .y = node->lagpos.y, .z = node->lagpos.z};
      peskin_kernel (lagpos) {
        foreach_dimension () {
#if dimension == 1
          tmp_force.x[] += (peskin_weight / (Delta)) * node->w.x * node->weight;
#elif dimension == 2
          tmp_force.x[] +=
            (peskin_weight / sq (Delta)) * node->w.x * node->weight;
#else
          tmp_force.x[] +=
            (peskin_weight / cube (Delta)) * node->w.x * node->weight;
#endif
        }
      }
    }
  }

#if _MPI
  synchronize ((scalar*) {tmp_force});
#endif

  // 2. c = M_f^{-1} g ~ g / (rho_f)
  foreach () {
    foreach_dimension () {
#if dimension == 1
      tmp_vel.x[] = tmp_force.x[] / (1 * RHO_F);
#elif dimension == 2
      tmp_vel.x[] = tmp_force.x[] / (1 * RHO_F);
#else
      tmp_vel.x[] = tmp_force.x[] / (1 * RHO_F);
#endif
    }
  }

  // 3. d = J c, then scale by dt
  foreach_ibnode (mesh) {
    foreach_dimension () {
      node->Ay.x = 0.;
    }
    foreach_cache (node->stencil) {
      coord lagpos = {
        .x = node->lagpos.x, .y = node->lagpos.y, .z = node->lagpos.z};
      peskin_kernel (lagpos) {
        foreach_dimension () {
          node->Ay.x += peskin_weight * tmp_vel.x[];
        }
      }
    }
    foreach_dimension () {
      node->Ay.x *= dt;
    }
  }
}

/*!
 * @memberof IBMesh
 */
trace void ib_mesh_solve_lambda_CG (IBMesh* mesh, double dt) {
  const int maxiter = 50000;
  const double tol = 1e-9;

  // Initial guess
  foreach_ibnode (mesh) {
    foreach_dimension () {
      node->force.x = 0.;        // lambda^0 = 0
      node->res.x = node->rhs.x; // res^0 = b
      node->w.x = node->res.x;   // w^0 = res^0
    }
  }

  double rr_old = 0.;
  foreach_ibnode (mesh) {
    foreach_dimension () {
      rr_old += sq (node->res.x) * node->weight;
    }
  }

  // printf("=======================================================\n");
  if (sqrt (rr_old) < tol) {
    return;
  }

  for (int k = 0; k < maxiter; k++) {
    // printf("k: %d : ||R|| = %f\n", k, rr_old);
    ib_mesh_matvec_Aw (mesh, dt);

    double wy = 0.;
    foreach_ibnode (mesh) {
      foreach_dimension () {
        wy += node->w.x * node->Ay.x * node->weight;
      }
    }

    if (fabs (wy) < 1e-30) {
      break;
    }
    double alpha = rr_old / wy;
    double rr_new = 0.;

    // \f(\lambda_{k+1} = \lambda_{k} + \alpha w_{k}\f), and \f(r_{k+1} = r_{k}
    // - \alpha y_k\f).
    foreach_ibnode (mesh) {
      foreach_dimension () {
        node->force.x += alpha * node->w.x; // lambda update
        node->res.x -= alpha * node->Ay.x;  // residual update
        rr_new += sq (node->res.x) * node->weight;
      }
    }

    // Check tolerance
    if (sqrt (rr_new) < tol) {
      break;
    }

    double beta = rr_new / rr_old;
    rr_old = rr_new;

    // w_{k+1} = r_{k+1} + \beta w_{k}
    foreach_ibnode (mesh) {
      foreach_dimension () {
        node->w.x = node->res.x + beta * node->w.x;
      }
    }
  }
}
