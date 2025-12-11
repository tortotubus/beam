#include "library/ibm/IBMesh.h"

/*
 * Forward declaration of this header's functions
 */
void ib_mesh_manager_init (int nm);
void ib_mesh_manager_free ();
void ib_mesh_manager_generate_stencil_caches ();
void ib_mesh_manager_tag_stencils ();

/*
 * Type definitions
 */

/*!
 * @brief Small struct to help us perform MPI reductions
 */
typedef struct {
  int pid;
  int nn;
  int mesh_idx;
} IBMeshTriplet;

// typedef struct {
//   int pid;
//   int node_idx;
//   int mesh_idx;
// } IBNodeTriplet;


/*!
 * @brief Management struct to hold all IBMeshes in the simulation
 */
typedef struct {
  IBMesh* meshes; /*!< The actual immersed boundary meshes */
  int nm;         /*!< Total nuber of meshes*/

  IBMeshTriplet* triplets; /*!< Meta-information used for MPI Allgatherv */
  int* nodes_per_pid;      /*!< Per-proc count of the nodes across all meshes */
} IBMeshManager;

/*
 * Globals
 */

IBMeshManager ib_mesh_manager;

/*
 * Macro definitions
 */

#define ibmesh(mesh_index) (&ib_mesh_manager.meshes[mesh_index])

macro foreach_ibmesh () {
  for (int mesh_index = 0; mesh_index < ib_mesh_manager.nm; mesh_index++) {
    IBMesh* mesh = ibmesh (mesh_index);
    NOT_UNUSED (mesh);
    // clang-format off
    {...}
    // clang-format on
  }
}

/*
 * Function definitions
 */

/*!
 * @brief Compare the mesh triplets
 *
 * @param a The first mesh triplet
 * @param b The second mesh triplet
 *
 * @returns 1 if a > b, -1 if a < b
 *
 * @memberof IBMeshManager
 */
int ib_mesh_manager_compare_triplets (const void* a, const void* b) {

  const IBMeshTriplet* aa = (IBMeshTriplet*) a;
  const IBMeshTriplet* bb = (IBMeshTriplet*) b;

  if (aa->pid < bb->pid)
    return -1;
  if (aa->pid > bb->pid)
    return 1;
  if (aa->mesh_idx < bb->mesh_idx)
    return -1;
  if (aa->mesh_idx > bb->mesh_idx)
    return 1;

  return 0;
}

/*!
 * @brief Set the owners of each IBMesh
 * @memberof IBMeshManager
 */
// void ib_mesh_manager_set_owners () {
//   foreach_ibmesh () {
//     mesh->pid = 0;
//   }
// }

/*!
 * @brief Allocate the mesh triplets
 *
 * @memberof IBMeshManager
 */
void ib_mesh_manager_init_triplets () {
  if (!ib_mesh_manager.triplets) {
    // Malloc the triplets
    ib_mesh_manager.triplets =
      (IBMeshTriplet*) malloc (ib_mesh_manager.nm * sizeof (IBMeshTriplet));
  }

  if (!ib_mesh_manager.nodes_per_pid) {
    // Malloc the nodes-per-pid member
    ib_mesh_manager.nodes_per_pid = (int*) calloc (npe (), sizeof (int));
  }

  // Reset counts every time
  for (int p = 0; p < npe (); p++)
    ib_mesh_manager.nodes_per_pid[p] = 0;

  // Copy information from the mesh array into the triplets
  for (int mi = 0; mi < ib_mesh_manager.nm; mi++) {
    int pid = ib_mesh_manager.meshes[mi].pid;
    int nn = ib_mesh_manager.meshes[mi].nn;

    ib_mesh_manager.triplets[mi].pid = pid;
    ib_mesh_manager.triplets[mi].nn = nn;
    ib_mesh_manager.triplets[mi].mesh_idx = mi;

    ib_mesh_manager.nodes_per_pid[pid] += nn;
  }
}

/*!
 * @brief Sort the mesh triplets
 *
 * @memberof IBMeshManager
 */
void ib_mesh_manager_sort_triplets () {
  // Allocate the triplets if not already present
  ib_mesh_manager_init_triplets ();

  // Sort the triplets
  qsort (ib_mesh_manager.triplets,
         ib_mesh_manager.nm,
         sizeof (IBMeshTriplet),
         ib_mesh_manager_compare_triplets);
}

/*!
 * @brief Compute the counts, displacements for allgatherv
 *
 * @memberof IBMeshManager
 */
void ib_mesh_manager_allgatherv_layout (int** recvcounts,
                                        int** displs,
                                        int* total,
                                        int comp_per_node) {
  int num_proc = npe ();
  *recvcounts = (int*) malloc (num_proc * sizeof (int));
  *displs = (int*) malloc (num_proc * sizeof (int));

  int offset = 0;

  for (int p = 0; p < num_proc; p++) {
    int n_nodes_p = ib_mesh_manager.nodes_per_pid[p];
    (*recvcounts)[p] = n_nodes_p * comp_per_node;
    (*displs)[p] = offset;
    offset += *recvcounts[p];
  }

  *total = offset;
}

/*!
 * @brief Initialize the immersed boundary mesh manager
 *
 * Although the meshes are default constructed and can be changed after calling
 * this, the number of meshes must remain fixed throughout the simulation. Since
 * meshes can be set to active or inactive, it is reccomended to choose the
 * number of meshes as the maximum expected to be needed throughout the life of
 * the simulation.
 *
 * @param nm The number of meshes you plan to have
 *
 * @memberof IBMeshManager
 */
void ib_mesh_manager_init (int nm) {
#if DEBUG
  printf ("ib_mesh_manager_init\n");
#endif
  if (ib_mesh_manager.meshes != NULL) {
    ib_mesh_manager_free ();
  }

  ib_mesh_manager.meshes = NULL;
  ib_mesh_manager.nm = nm;
  ib_mesh_manager.triplets = NULL;
  ib_mesh_manager.nodes_per_pid = NULL;

  ib_mesh_manager.meshes = (IBMesh*) calloc (nm, sizeof (IBMesh));

  // Default constructor for each of our meshes
  foreach_ibmesh () {
    ib_mesh_init (mesh);
  }
}

/*!
 * @brief Free all members in the immersed boundary mesh manager, including the
 * array of meshes, their nodes, and stencil nodes.
 *
 * @memberof IBMeshManager
 */
void ib_mesh_manager_free () {
#if DEBUG
  printf ("ib_mesh_manager_free\n");
#endif

  // Call the IBMesh destructor on each mesh
  foreach_ibmesh () {
    ib_mesh_free (mesh);
  }

  // Free the owned arrays
  free (ib_mesh_manager.meshes);
  free (ib_mesh_manager.triplets);
  free (ib_mesh_manager.nodes_per_pid);

  ib_mesh_manager.meshes = NULL;
  ib_mesh_manager.nm = 0;
  ib_mesh_manager.triplets = NULL;
  ib_mesh_manager.nodes_per_pid = NULL;
}

/*!
 * @brief For force-coupled meshes, we advance the mesh in time according to
 * it's own kinematics using the previous force at time \f(t^{n}\f)
 *
 * @memberof IBMeshManager
 */
trace void ib_mesh_manager_advance_lagrangian_mesh () {
#if DEBUG
  printf ("ib_mesh_manager_advance_lagrangian_mesh");
#endif
  foreach_ibmesh () {
    if (pid () == mesh->pid)
      ib_mesh_advance_lagrangian_mesh (mesh);
  }

#if _MPI
  ib_mesh_manager_sort_triplets ();

  int* recvcounts = NULL;
  int* displs = NULL;
  int total = 0;
  int components = 2 * dimension;

  ib_mesh_manager_allgatherv_layout (&recvcounts, &displs, &total, components);

  double* buf = calloc (total, sizeof (double));
  int local_pid = pid ();
  int idx = displs[local_pid];

  for (int ti = 0; ti < ib_mesh_manager.nm; ti++) {
    int tpid = ib_mesh_manager.triplets[ti].pid;
    int tidx = ib_mesh_manager.triplets[ti].mesh_idx;
    if (tpid == pid ()) {
      foreach_ibnode (ibmesh (tidx)) {
        foreach_dimension () {
          buf[idx++] = node->lagpos.x;
          buf[idx++] = node->lagvel.x;
        }
      }
    }
  }

  MPI_Allgatherv (MPI_IN_PLACE,      // sendbuf
                  0,                 // sendcount
                  MPI_DATATYPE_NULL, // datatype_send
                  buf,               // recvbuf
                  recvcounts,        // counts_recv
                  displs,            // displacements
                  MPI_DOUBLE,        // datatype_recv
                  MPI_COMM_WORLD     // communicator
  );

  for (int ti = 0; ti < ib_mesh_manager.nm; ti++) {
    int tpid = ib_mesh_manager.triplets[ti].pid;
    int tidx = ib_mesh_manager.triplets[ti].mesh_idx;
    if (tpid != pid ()) {
      idx = displs[tpid];
      foreach_ibnode (ibmesh (tidx)) {
        foreach_dimension () {
          node->lagpos.x = buf[idx++];
          node->lagvel.x = buf[idx++];
        }
      }
    }
  }
#endif
}

/*!
 * @brief Generate stencils for every node of every mesh
 *
 * @memberof IBMeshManager
 */
trace void ib_mesh_manager_generate_stencil_cache () {
#if DEBUG
  printf ("ib_mesh_manager_generate_stencil_cache\n");
#endif
  foreach_ibmesh () {
    if (mesh->isactive) {
      ib_mesh_generate_stencil_cache (mesh);
    }
  }
}

/*!
 * @brief
 *
 * @memberof IBMeshManager
 */
trace void ib_mesh_manager_spread_eulerian_forcing (vector forcing) {
  foreach_ibmesh () {
    if (mesh->isactive) {
      ib_mesh_spread_eulerian_forcing (mesh, forcing);
    }
  }
}

/*!
 * @brief
 *
 * @memberof IBMeshManager
 */
trace void ib_mesh_manager_interpolate_eulerian_velocities () {
  foreach_ibmesh () {
    if (mesh->isactive) {
      ib_mesh_interpolate_eulerian_velocities (mesh);
    }
  }

#if 0

  int total_nodes = 0;
  int nprocs = npe (); 
  int components = dimension;
  
  for (int p = 0; p < nprocs; p++)
    total_nodes += ib_mesh_manager.nodes_per_pid[p];

  double* buf = calloc (total_nodes * components, sizeof (double));

  int idx = 0;
  
  foreach_ibmesh () {
    foreach_ibnode(mesh) {
      foreach_dimension() {
        buf[idx++] = node->eulvel.x;
      }
    }
  }

  MPI_Allreduce (MPI_IN_PLACE,
                 buf,
                 total_nodes * components,
                 MPI_DOUBLE,
                 MPI_SUM,
                 MPI_COMM_WORLD);

  idx = 0;
  foreach_ibmesh () {
    foreach_ibnode(mesh) {
      foreach_dimension() {
        node->eulvel.x = buf[idx++];
      }
    }
  }

  free (buf);

#endif
}

/*!
 * @brief
 *
 * @memberof IBMeshManager
 */
trace void ib_mesh_manager_tag_stencil () {
#if DEBUG
  printf ("tag_stencils\n");
#endif

  foreach () {
    stencils[] = 0.;
  }

  foreach_ibmesh () {
    if (mesh->isactive) {
      ib_mesh_tag_stencil (mesh);
    }
  }
}

/*!
 * @brief
 *
 * @memberof IBMeshManager
 */
trace void ib_mesh_manager_compute_constraint_rhs () {
#if DEBUG
  printf ("tag_stencils\n");
#endif
  foreach_ibmesh () {
    if (mesh->isactive) {
      ib_mesh_compute_constraint_rhs (mesh);
    }
  }
}

/*!
 * @brief
 *
 * @memberof IBMeshManager
 */

// 0.475972
trace void ib_mesh_manager_solve_lambda_CG (double dt) {
#if DEBUG
  printf ("tag_stencils\n");
#endif

  const int maxiter = 50000;
  const double tol = 1e-9;

  // Initial guess
  foreach_ibmesh () {
    foreach_ibnode (mesh) {
      foreach_dimension () {
        node->force.x = 0;
        node->res.x = node->rhs.x;
        node->w.x = node->res.x;
      }
    }
  }

  double rr_old = 0.;
  foreach_ibmesh () {
    foreach_ibnode (mesh) {
      foreach_dimension () {
        rr_old += sq (node->res.x)* node->weight;
      }
    }
  }

#if _MPI
  double rr_old_local = rr_old;
  MPI_Allreduce (
    &rr_old_local, &rr_old, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  if (sqrt (rr_old) < tol) {
    return;
  }

  for (int k = 0; k < maxiter; k++) {
    foreach_ibmesh () {
      ib_mesh_matvec_Aw (mesh, dt);
    }

    double wy = 0.;
    foreach_ibmesh () {
      foreach_ibnode (mesh) {
        foreach_dimension () {
          wy += node->w.x * node->Ay.x * node->weight;
        }
      }
    }

#if _MPI
    double wy_local = wy;
    MPI_Allreduce (&wy_local, &wy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    if (fabs (wy) < 1e-30)
      break;

    double alpha = rr_old / wy;
    double rr_new = 0.;

    // \f(\lambda_{k+1} = \lambda_{k} + \alpha w_{k}\f), and \f(r_{k+1} = r_{k}
    // - \alpha y_k\f).
    foreach_ibmesh () {
      foreach_ibnode (mesh) {
        foreach_dimension () {
          node->force.x += alpha * node->w.x; // lambda udate
          node->res.x -= alpha * node->Ay.x;  // residual update
          rr_new += sq (node->res.x)* node->weight;
        }
      }
    }

#if _MPI
    double rr_new_local = rr_new;
    MPI_Allreduce (
      &rr_new_local, &rr_new, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    // Check tolerance
    if (sqrt (rr_new) < tol)
      break;

    double beta = rr_new / rr_old;
    rr_old = rr_new;

    // w_{k+1} = r_{k+1} + \beta w_{k}
    foreach_ibmesh () {
      foreach_ibnode (mesh) {
        foreach_dimension () {
          node->w.x = node->res.x + beta * node->w.x  ;
        }
      }
    }
  }
}