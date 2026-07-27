#pragma once

#include <limits.h>

#include "library/ibm/IBNode.h"
#include "library/ibm/IBMempool.h"
#include "library/ibm/IBMesh.h"
#include "library/ibm/IBMeshManagerMPI.h"

#ifndef MPI_AUTO_BOUNDARY
#define MPI_AUTO_BOUNDARY 0
#endif

#ifndef MPI_DEBUG
#define MPI_DEBUG 0
#endif

/**
 * Type definitions
 */

/**
 * @struct IBMeshManager
 */
typedef struct {
  IBMesh* meshes;
  IBMempool pool;
  int nm;
  bool dirty;
#if _MPI
  IBNodeList local;
  IBExchangeList* snd_migrate;
  IBExchangeList* rcv_migrate;
  IBExchangeList* snd_boundary;
  IBExchangeList* rcv_boundary;
#if !TREE
  MPI_Comm cartcomm;
#endif
#endif
} IBMeshManager;

// ============================================================================
// Globals
// ============================================================================

IBMeshManager ibmm = {0};
static bool ibmeshmanager_use_velocity_midpoint_op = false;

// ============================================================================
// Function declarations
// ============================================================================

void ibmeshmanager_init (int mesh_count);
void ibmeshmanager_free ();
int ibmeshmanager_add_mesh ();
void ibmeshmanager_delete_mesh (int mesh_id);
void ibmeshmanager_add_nodes (int mesh_id, int count);
void ibmeshmanager_delete_all_nodes (int mesh_id);
void ibmeshmanager_set_model (int mesh_id, IBMeshModel model);
void ibmeshmanager_advance_positions (double dt);
void ibmeshmanager_advance_force_coupled_positions (double dt);
void ibmeshmanager_evaluate_velocity_coupled_midpoints (double dt);
void ibmeshmanager_advance_velocity_coupled_positions (double dt);

#if _MPI
int _ibmeshmanager_get_pid (Point p);
void ibmeshmanager_update_pid ();
void ibmeshmanager_boundary (IBscalar* list = iball);
static void ibmeshmanager_check_exchange_counts (const char* kind,
                                                 IBExchangeList* snd,
                                                 IBExchangeList* rcv);
#endif

// ============================================================================
// Macros
// ============================================================================

#ifndef IBMESHMANAGER_POOL_SIZE_BYTES
#define IBMESHMANAGER_POOL_SIZE_BYTES (1 << 19)
#endif

/**
 * @def foreach_ibmesh
 * @brief Loops through all meshes
 * @relates IBMeshManager
 */
macro foreach_ibmesh () {
  for (int mesh_id = 0; mesh_id < ibmm.nm; mesh_id++) {
    IBMesh* mesh = &ibmm.meshes[mesh_id];
    NOT_UNUSED (mesh);
    // clang-format off
    {...}
    // clang-format on
  }
}

/**
 * @def foreach_ibnode
 * @brief Loops through all nodes in the node pool
 * @relates IBMeshManager
 */
macro foreach_ibnode (bool local_only = false) {
#if _MPI
  {
    IBNodeList nlist = local_only ? ibmm.local : ibmm.pool.active;
    for (size_t node_id = 0; node_id < nlist.size; node_id++) {
      IBNode* node = nlist.ptrs[node_id];
      NOT_UNUSED (node);

      coord pos = {0};
      foreach_dimension () pos.x = ibval (npos.x);
      coord_periodic_boundary (pos);
      NOT_UNUSED (pos);

      // bool ib_set_dirty = false;
      // NOT_UNUSED (ib_set_dirty);

      // clang-format off
      {...}
      // clang-format on
    } 
  }
#else
  for (size_t node_id = 0; node_id < ibmm.pool.active.size; node_id++) {
    IBNode* node = ibmm.pool.active.ptrs[node_id];
    NOT_UNUSED (node);

    coord pos = {0};
    foreach_dimension () pos.x = ibval (npos.x);
    coord_periodic_boundary (pos);
    NOT_UNUSED (pos);

    // bool ib_set_dirty = false;
    // NOT_UNUSED (ib_set_dirty);

    // clang-format off
    {...}
    // clang-format on
  }
#endif
}

/**
 * @def foreach_ibnode_per_ibmesh
 * @brief Loops through all nodes of all meshes
 * @relates IBMeshManager
 */
macro foreach_ibnode_per_ibmesh () {
  for (size_t mesh_id = 0; mesh_id < ibmm.nm; mesh_id++) {
    IBMesh* mesh = &ibmm.meshes[mesh_id];
    NOT_UNUSED (mesh);
    for (size_t node_id = 0; node_id < mesh->nodes.size; node_id++) {
      IBNode* node = mesh->nodes.ptrs[node_id];

      NOT_UNUSED (node);
      coord pos = {0};
      foreach_dimension()
        pos.x = ibval(npos.x);
      coord_periodic_boundary(pos);      
      NOT_UNUSED (pos);

      // bool ib_set_dirty = false;
      // NOT_UNUSED (ib_set_dirty);
      // clang-format off
      {...}
      // clang-format on
    }
  }
}

/* Function definitions */

/**
 * @brief Initialize the immersed boundary mesh manager
 * @param mesh_count The number of meshes you plan to have
 * @relates IBMeshManager
 */
void ibmeshmanager_init (int mesh_count) {

  init_ibsolver ();

  if (!ibmm.meshes) {
    foreach_ibmesh () ibmesh_free (mesh);
    free (ibmm.meshes);
  }

  ibmm.meshes = NULL;
  ibmm.nm = mesh_count;
  ibmm.meshes = (IBMesh*) calloc (mesh_count, sizeof (IBMesh));
  ibmm.pool =
    ibmempool_init (IBMESHMANAGER_POOL_SIZE_BYTES, nibvar * sizeof (real));

#if _MPI
  ibmm.dirty = true;
  ibnodelist_init (&ibmm.local, 0);
  ibmm.snd_boundary =
    (IBExchangeList*) calloc (npe (), sizeof (IBExchangeList));
  ibmm.rcv_boundary =
    (IBExchangeList*) calloc (npe (), sizeof (IBExchangeList));
  ibmm.snd_migrate = (IBExchangeList*) calloc (npe (), sizeof (IBExchangeList));
  ibmm.rcv_migrate = (IBExchangeList*) calloc (npe (), sizeof (IBExchangeList));
  for (int i = 0; i < npe (); i++) {
    ibexchangelist_init (&ibmm.snd_boundary[i], 0);
    ibexchangelist_set_pid (&ibmm.snd_boundary[i], i);
    ibexchangelist_init (&ibmm.rcv_boundary[i], 0);
    ibexchangelist_set_pid (&ibmm.rcv_boundary[i], i);
    ibexchangelist_init (&ibmm.snd_migrate[i], 0);
    ibexchangelist_set_pid (&ibmm.snd_migrate[i], i);
    ibexchangelist_init (&ibmm.rcv_migrate[i], 0);
    ibexchangelist_set_pid (&ibmm.rcv_migrate[i], i);
  }

#if !TREE
  int dims[dimension];
  int periods[dimension];
  dims[0] = Dimensions.x;
  periods[0] = Period.x;
#if dimension >= 2
  dims[1] = Dimensions.y;
  periods[1] = Period.y;
#endif
#if dimension >= 3
  dims[2] = Dimensions.z;
  periods[2] = Period.z;
#endif
  MPI_Cart_create (MPI_COMM_WORLD, dimension, dims, periods, 0, &ibmm.cartcomm);
#endif // !TREE

#endif // _MPI

  foreach_ibmesh () {
    ibmesh_init (mesh);
  }
}

/**
 * @brief Free all members in the immersed boundary mesh manager.
 * @relates IBMeshManager
 */
void ibmeshmanager_free () {
  IBMempool* pool = &ibmm.pool;
  foreach_ibmesh () {
    ibmesh_free (mesh);
  }
  free (ibmm.meshes);
  ibmempool_free (&ibmm.pool);
  ibmm.meshes = NULL;
  ibmm.nm = 0;

#if _MPI
  ibnodelist_free (&ibmm.local);
  for (int i = 0; i < npe (); i++) {
    ibexchangelist_free (&ibmm.snd_boundary[i]);
    ibexchangelist_free (&ibmm.rcv_boundary[i]);
    ibexchangelist_free (&ibmm.snd_migrate[i]);
    ibexchangelist_free (&ibmm.rcv_migrate[i]);
  }
  free (ibmm.snd_boundary);
  ibmm.snd_boundary = NULL;
  free (ibmm.rcv_boundary);
  ibmm.rcv_boundary = NULL;
  free (ibmm.snd_migrate);
  ibmm.snd_migrate = NULL;
  free (ibmm.rcv_migrate);
  ibmm.rcv_migrate = NULL;
#endif
}

/**
 * @brief Creates a new mesh returning the index
 * @relates IBMeshManager
 */
int ibmeshmanager_add_mesh () {
  int mesh_id = ibmm.nm;
  ibmm.nm++;
  ibmm.meshes = (IBMesh*) realloc (ibmm.meshes, ibmm.nm * sizeof (IBMesh));
  assert (ibmm.meshes);
  ibmesh_init (&ibmm.meshes[mesh_id]);
  return mesh_id;
}

/**
 * @brief Deletes the mesh from the manager, freeing the mesh object, and
 * marking as free its nodes in the pool.
 * @relates IBMeshManager
 */
void ibmeshmanager_delete_mesh (int mesh_id) {
  assert (mesh_id >= 0 && mesh_id < ibmm.nm);
  IBMempool* pool = &ibmm.pool;
  IBMesh* mesh = &ibmm.meshes[mesh_id];
  ibmesh_delete_all_nodes (mesh, pool);
  ibmesh_free (mesh);

  int last = ibmm.nm - 1;
  if (mesh_id != last)
    ibmm.meshes[mesh_id] = ibmm.meshes[last];
  ibmm.nm--;
}

/**
 * @brief Bulk adds nodes to a given mesh
 * @relates IBMeshManager
 */
void ibmeshmanager_add_nodes (int mesh_id, int count) {
  assert (ibmm.nm > mesh_id);
  IBMesh* mesh = &ibmm.meshes[mesh_id];
  IBMempool* pool = &ibmm.pool;
  ibmesh_add_nodes (mesh, pool, count);
}

/**
 * @brief Deletes all nodes of a given mesh
 * @relates IBMeshManager
 */
void ibmeshmanager_delete_all_nodes (int mesh_id) {
  IBMesh* mesh = &ibmm.meshes[mesh_id];
  IBMempool* pool = &ibmm.pool;
  ibmesh_delete_all_nodes (mesh, pool);
}

/**
 * @brief Set the mesh model
 * @relates IBMeshManager
 */
void ibmeshmanager_set_model (int mesh_id, IBMeshModel model) {
  IBMesh* mesh = &ibmm.meshes[mesh_id];
  IBMempool* pool = &ibmm.pool;
  ibmesh_set_model (mesh, pool, model);
}

/**
 * @brief
 * @relates IBMeshManager
 */
trace void ibmeshmanager_advance_positions_filtered (double dt, int model_type) {
#if _MPI
  ibmeshmanager_update_pid ();
#endif

#if _MPI
  // Synchronize inputs
  {
    size_t stride = dimension;
    size_t n = ibmm.pool.active.size * stride;
    double* f_global = calloc (n, sizeof (double));
    double* vel_global = calloc (n, sizeof (double));

    foreach_ibnode_per_ibmesh () {
      if (model_type != IB_MODEL_INVALID && mesh->model.type != model_type)
        continue;
      if (node->pid == pid ()) {
        int di = 0;
        foreach_dimension () {
          f_global[node_id * stride + di] = ibval (nforce.x);
          vel_global[node_id * stride + di] = ibval (nvel.x);
          di++;
        }
      }
    }

    MPI_Allreduce (
      MPI_IN_PLACE, f_global, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (
      MPI_IN_PLACE, vel_global, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // if (pid () == 0) {
    //   printf ("[proc %d]: ", pid ());
    //   for (int i = 0; i < n; i++) {
    //     printf ("%f ", f_global[i]);
    //   }
    //   printf ("\n");
    // }

    foreach_ibnode_per_ibmesh () {
      if (model_type != IB_MODEL_INVALID && mesh->model.type != model_type)
        continue;
      const int owner = mesh->pid >= 0 ? mesh->pid : 0;
      if (owner == pid ()) {
        int di = 0;
        foreach_dimension () {
          ibval (nforce.x) = f_global[node_id * stride + di];
          ibval (nvel.x) = vel_global[node_id * stride + di];
          di++;
        }
      }
    }

    free (f_global);
    free (vel_global);
  }
#endif

  // Call the models
  foreach_ibmesh () {
    if (model_type != IB_MODEL_INVALID && mesh->model.type != model_type)
      continue;
    switch (mesh->model.type) {
    case IB_MODEL_VELOCITY_COUPLED: {
      const int owner = mesh->pid >= 0 ? mesh->pid : 0;
      if (owner == pid ()) {
        IBVelocityCoupledModelOps* ops = mesh->model.velocity_ops;
        void* ctx = mesh->model.ctx;
        if (ibmeshmanager_use_velocity_midpoint_op && ops->midpoint)
          ops->midpoint (ctx, mesh, dt);
        else
          ops->advance (ctx, mesh, dt);
      }
      break;
    }
    case IB_MODEL_FORCE_COUPLED: {
      const int owner = mesh->pid >= 0 ? mesh->pid : 0;
      if (owner == pid ()) {
        IBForceCoupledModelOps* ops = mesh->model.force_ops;
        void* ctx = mesh->model.ctx;
        ops->advance (ctx, mesh, dt);
      }
      break;
    }
    default:
      break;
    }
  }

#if _MPI
  // Syncronize outputs
  foreach_ibmesh () {
    if (model_type != IB_MODEL_INVALID && mesh->model.type != model_type)
      continue;
    switch (mesh->model.type) {
    case IB_MODEL_VELOCITY_COUPLED: {
      const int owner = mesh->pid >= 0 ? mesh->pid : 0;
      const int nn = (int) mesh->nodes.size;
      const int stride = 3 * dimension + 1;
      double* buff = malloc ((size_t) nn * stride * sizeof (double));

      if (owner == pid ()) {
        for (int ni = 0; ni < nn; ni++) {
          IBNode* node = mesh->nodes.ptrs[ni];
          int k = ni * stride;

          foreach_dimension () {
            buff[k++] = ibval (npos.x);
          }
          foreach_dimension () {
            buff[k++] = ibval (nvel.x);
          }
          foreach_dimension () {
            buff[k++] = ibval (nforce.x);
          }
          buff[k++] = ibval (nweight);
        }
      }

      MPI_Bcast (buff, nn * stride, MPI_DOUBLE, owner, MPI_COMM_WORLD);

      if (owner != pid ()) {
        for (int ni = 0; ni < nn; ni++) {
          IBNode* node = mesh->nodes.ptrs[ni];
          int k = ni * stride;
          foreach_dimension () {
            ibval (npos.x) = buff[k++];
          }
          foreach_dimension () {
            ibval (nvel.x) = buff[k++];
          }
          foreach_dimension () {
            ibval (nforce.x) = buff[k++];
          }
          ibval (nweight) = buff[k++];
        }
      }

      free (buff);
      break;
    }
    case IB_MODEL_FORCE_COUPLED: {
      const int owner = mesh->pid >= 0 ? mesh->pid : 0;
      const int nn = (int) mesh->nodes.size;
      const int stride = 2 * dimension;
      double* buff = malloc ((size_t) nn * stride * sizeof (double));

      if (owner == pid ()) {
        // Pack
        for (int ni = 0; ni < nn; ni++) {
          IBNode* node = mesh->nodes.ptrs[ni];
          int k = ni * stride;

          foreach_dimension () {
            buff[k++] = ibval (npos.x);
          }
          foreach_dimension () {
            buff[k++] = ibval (nvel.x);
          }
        }
      }

      MPI_Bcast (buff, nn * stride, MPI_DOUBLE, owner, MPI_COMM_WORLD);

      if (owner != pid ()) {
        // Free whatever model object we have, we will not call it later
        // ibmeshmodel_destroy (&mesh->model);

        // Unpack
        for (int ni = 0; ni < nn; ni++) {
          IBNode* node = mesh->nodes.ptrs[ni];
          int k = ni * stride;
          foreach_dimension () {
            ibval (npos.x) = buff[k++];
          }
          foreach_dimension () {
            ibval (nvel.x) = buff[k++];
          }
        }
      }

      free (buff);
      break;
    }

    default: {
      break;
    }
    }
  }

  ibmm.dirty = true;
  ibmeshmanager_update_pid ();
#endif
}

trace void ibmeshmanager_advance_positions (double dt) {
  ibmeshmanager_advance_positions_filtered (dt, IB_MODEL_INVALID);
}

trace void ibmeshmanager_advance_force_coupled_positions (double dt) {
  ibmeshmanager_advance_positions_filtered (dt, IB_MODEL_FORCE_COUPLED);
}

trace void ibmeshmanager_evaluate_velocity_coupled_midpoints (double dt) {
  ibmeshmanager_use_velocity_midpoint_op = true;
  ibmeshmanager_advance_positions_filtered (dt, IB_MODEL_VELOCITY_COUPLED);
  ibmeshmanager_use_velocity_midpoint_op = false;
}

trace void ibmeshmanager_advance_velocity_coupled_positions (double dt) {
  ibmeshmanager_use_velocity_midpoint_op = false;
  ibmeshmanager_advance_positions_filtered (dt, IB_MODEL_VELOCITY_COUPLED);
}

#if _MPI

trace inline int _ibmeshmanager_get_pid (Point p) {
#if TREE
  {
    Point point = p;
    int ig = 0, jg = 0, kg = 0;
    NOT_UNUSED (ig);
    NOT_UNUSED (jg);
    NOT_UNUSED (kg);
    POINT_VARIABLES ();

    if (allocated (0))
      return cell.pid;
    else
      return -1;
  }
#else
  /**
   * Here since we:
   *  - cell.pid is not defined
   *  - is_local returns true always
   *  - Domain info is handled by MPI_Cart
   *  - MPI_Comm cartcomm is not easily accessible
   * we instead recreate our own MPI_Cart usign the same parameters and keep it
   * static so that we may query about MPI/Phyiscal domain neighbors.
   */
  {
    Point point = p;

    // if (!is_boundary (point))
    //   return pid ();

    int coords[dimension];
    coords[0] = mpi_coords[0];
    if (point.i < GHOSTS)
      coords[0]--;
    else if (point.i >= point.n.x + GHOSTS)
      coords[0]++;

#if dimension > 1
    coords[1] = mpi_coords[1];
    if (point.j < GHOSTS)
      coords[1]--;
    else if (point.j >= point.n.y + GHOSTS)
      coords[1]++;
#endif

#if dimension > 2
    coords[2] = mpi_coords[2];
    if (point.k < GHOSTS)
      coords[2]--;
    else if (point.k >= point.n.z + GHOSTS)
      coords[2]++;
#endif

    if (!Period.x && (coords[0] < 0 || coords[0] >= Dimensions.x))
      return -1;
    if (Period.x) {
      if (coords[0] < 0)
        coords[0] += Dimensions.x;
      else if (coords[0] >= Dimensions.x)
        coords[0] -= Dimensions.x;
    }
#if dimension > 1
    if (!Period.y && (coords[1] < 0 || coords[1] >= Dimensions.y))
      return -1;
    if (Period.y) {
      if (coords[1] < 0)
        coords[1] += Dimensions.y;
      else if (coords[1] >= Dimensions.y)
        coords[1] -= Dimensions.y;
    }
#endif
#if dimension > 2
    if (!Period.z && (coords[2] < 0 || coords[2] >= Dimensions.z))
      return -1;
    if (Period.z) {
      if (coords[2] < 0)
        coords[2] += Dimensions.z;
      else if (coords[2] >= Dimensions.z)
        coords[2] -= Dimensions.z;
    }
#endif

    int owner = -1;
    MPI_Cart_rank (ibmm.cartcomm, coords, &owner);
    return owner;
  }
#endif
}

static void ibmeshmanager_check_exchange_counts (const char* kind,
                                                 IBExchangeList* snd,
                                                 IBExchangeList* rcv) {
  const int np = npe ();
  unsigned long long* send_counts =
    (unsigned long long*) calloc ((size_t) np, sizeof (unsigned long long));
  unsigned long long* recv_counts =
    (unsigned long long*) calloc ((size_t) np, sizeof (unsigned long long));
  unsigned long long* all_send_counts =
    (unsigned long long*) calloc ((size_t) np * np, sizeof (unsigned long long));
  unsigned long long* all_recv_counts =
    (unsigned long long*) calloc ((size_t) np * np, sizeof (unsigned long long));
  assert (send_counts && recv_counts && all_send_counts && all_recv_counts);

  for (int peer = 0; peer < np; peer++) {
    send_counts[peer] = (unsigned long long) snd[peer].nodes.size;
    recv_counts[peer] = (unsigned long long) rcv[peer].nodes.size;
  }

  MPI_Allgather (send_counts, np, MPI_UNSIGNED_LONG_LONG, all_send_counts, np,
                 MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
  MPI_Allgather (recv_counts, np, MPI_UNSIGNED_LONG_LONG, all_recv_counts, np,
                 MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);

  int mismatch_count = 0;
  for (int src = 0; src < np; src++) {
    for (int dst = 0; dst < np; dst++) {
      const unsigned long long sent = all_send_counts[src * np + dst];
      const unsigned long long posted = all_recv_counts[dst * np + src];
      if (sent != posted)
        mismatch_count++;
    }
  }

  if (mismatch_count) {
    if (pid () == 0) {
      int printed = 0;
      fprintf (stderr,
               "ERROR: inconsistent IB %s exchange counts: %d rank-pair "
               "mismatches\n",
               kind, mismatch_count);
      for (int src = 0; src < np && printed < 32; src++) {
        for (int dst = 0; dst < np && printed < 32; dst++) {
          const unsigned long long sent = all_send_counts[src * np + dst];
          const unsigned long long posted = all_recv_counts[dst * np + src];
          if (sent != posted) {
            fprintf (stderr,
                     "  %s src=%d dst=%d send_nodes=%llu "
                     "dst_recv_nodes_from_src=%llu\n",
                     kind, src, dst, sent, posted);
            printed++;
          }
        }
      }
    }
    free (send_counts);
    free (recv_counts);
    free (all_send_counts);
    free (all_recv_counts);
    MPI_Abort (MPI_COMM_WORLD, 4);
  }

  free (send_counts);
  free (recv_counts);
  free (all_send_counts);
  free (all_recv_counts);
}

/**
 * @brief Update all IBNode pids
 *
 * Locally queries the eulerian grid cell enclosed by a box
 * and sets the pid of the node to the pid of that cell
 *
 * @relates IBMeshManager
 */
trace void ibmeshmanager_update_pid () {
  if (!ibmm.dirty)
    return;

  ibnodelist_clear (&ibmm.local);

  for (int peer = 0; peer < npe (); peer++) {
    ibexchangelist_clear (&ibmm.rcv_boundary[peer]);
    ibexchangelist_clear (&ibmm.snd_boundary[peer]);
    ibexchangelist_clear (&ibmm.rcv_migrate[peer]);
    ibexchangelist_clear (&ibmm.snd_migrate[peer]);
  }

  const int nnode = (int) ibmm.pool.active.size;
  int* node_pids = (int*) malloc ((size_t) nnode * sizeof (int));
  const int support_mask_words = (npe () + 63) / 64;
  const size_t support_mask_count = (size_t) nnode * support_mask_words;
  unsigned long long* node_support_masks =
    (unsigned long long*) calloc (support_mask_count,
                                  sizeof (unsigned long long));
  assert (node_pids && node_support_masks);
  for (int i = 0; i < nnode; i++)
    node_pids[i] = -1;

  /* Pass 1: local owner and support candidates only. */
  foreach_ibnode () {
    IBNODE_VARIABLES();
    // coord d = {0};
    // foreach_dimension () d.x = ibval (npos.x);

    // coord_periodic_boundary (d);
#if TREE
    Point point = locate_nonlocal (pos.x, pos.y, pos.z);
    int ig = 0, jg = 0, kg = 0;
    NOT_UNUSED (ig);
    NOT_UNUSED (jg);
    NOT_UNUSED (kg);
    POINT_VARIABLES ();

    if (point.level >= 0) {
      if (allocated (0)) {
        if (is_local (cell)) {
          node_pids[node_id] = pid ();
        }
      }
    }
    point = locate_level (pos.x, pos.y, pos.z, node->depth);
    if (point.level >= 0) {
      bool has_local_support = false;
      foreach_neighbor (PESKIN_SUPPORT_RADIUS) {
        if (!has_local_support && allocated (0) && is_local (cell))
          has_local_support = true;
      }
      if (has_local_support) {
        const int word = pid () / 64;
        const int bit = pid () % 64;
        node_support_masks[node_id * support_mask_words + word] |=
          1ull << bit;
      }
    }
#else
    Point point = locate (d.x, d.y, d.z);
    if (point.level >= 0) {
      if (!is_boundary (point)) {
        node_pids[node_id] = pid ();
      }
    }
    point = locate_nonlocal (pos.x, pos.y, pos.z);
    if (point.level >= 0) {
      const int word = pid () / 64;
      const int bit = pid () % 64;
      node_support_masks[node_id * support_mask_words + word] |= 1ull << bit;
    }
#endif
  }

  // return;

  /* Synchronize globally: each node gets one owner and one support mask. */
  mpi_all_reduce_array (node_pids, MPI_INT, MPI_MAX, nnode);
  if (support_mask_count > 0) {
    assert (support_mask_count <= (size_t) INT_MAX);
    MPI_Allreduce (MPI_IN_PLACE, node_support_masks, (int) support_mask_count,
                   MPI_UNSIGNED_LONG_LONG, MPI_BOR, MPI_COMM_WORLD);
  }

  // return;

  /* Pass 2: build local/snd_boundary/rcv_boundary from resolved owners. */
  foreach_ibnode () {

    // Check if the owner has actually changed
    int old_pid = node->pid;
    int new_pid = node_pids[node_id];
#if _MPI
    if (new_pid < 0 || new_pid >= npe ()) {
      IBNODE_VARIABLES ();
      coord wrapped_pos = pos;
      coord_periodic_boundary (wrapped_pos);
#if TREE
      Point failed_point =
        locate_level (wrapped_pos.x, wrapped_pos.y, wrapped_pos.z, node->depth);
#else
      Point failed_point =
        locate_nonlocal (wrapped_pos.x, wrapped_pos.y, wrapped_pos.z);
#endif
      fprintf (stderr,
               "[rank %d] ERROR: unresolved IB node owner in "
               "ibmeshmanager_update_pid: node_id=%zu old_pid=%d "
               "new_pid=%d node_depth=%d raw_pos=(%g,%g,%g) "
               "wrapped_pos=(%g,%g,%g) locate_level=%d "
               "point=(%d,%d,%d) grid_depth=%d cells=%ld L0=%g "
               "origin=(%g,%g,%g)\n",
               pid (), node_id, old_pid, new_pid, node->depth, pos.x, pos.y,
               pos.z, wrapped_pos.x, wrapped_pos.y, wrapped_pos.z,
               failed_point.level, failed_point.i, failed_point.j,
               failed_point.k, depth (), grid->tn, L0, X0, Y0, Z0);
      MPI_Abort (MPI_COMM_WORLD, 3);
    }
#endif
    node->pid = new_pid;

    // If the new or old owner was us, we must exchange
    if (old_pid != -1 && old_pid != new_pid) {
      if (old_pid == pid ()) {
        ibexchangelist_push (&ibmm.snd_migrate[new_pid], node);
      }
      if (new_pid == pid ()) {
        ibexchangelist_push (&ibmm.rcv_migrate[old_pid], node);
      }
    }

    const int support_word = pid () / 64;
    const int support_bit = pid () % 64;
    const bool has_local_support =
      node_support_masks[node_id * support_mask_words + support_word] &
      (1ull << support_bit);
    if (node->pid == pid ()) { // local node
      ibnodelist_push (&ibmm.local, node);
      for (int peer = 0; peer < npe (); peer++) {
        const int peer_word = peer / 64;
        const int peer_bit = peer % 64;
        if (peer != pid () &&
            (node_support_masks[node_id * support_mask_words + peer_word] &
             (1ull << peer_bit)))
          ibexchangelist_push_unique (&ibmm.snd_boundary[peer], node);
      }
    } // local node
    else { // remote node
      if (has_local_support)
        ibexchangelist_push_unique (&ibmm.rcv_boundary[node->pid], node);
    } // remote node
  }

  free (node_support_masks);
  free (node_pids);

  // Migrate nodes
  IBscalar* slist = iball;
  size_t nscalars = iblist_len (slist);

  //ibmeshmanager_check_exchange_counts ("migration", ibmm.snd_migrate,
                                       ibmm.rcv_migrate);

  for (int peer = 0; peer < npe (); peer++) {
    if (peer != pid ()) {
      ibexchangelist_init_buffer (&ibmm.snd_migrate[peer], nscalars);
      ibexchangelist_init_buffer (&ibmm.rcv_migrate[peer], nscalars);
    }
  }

  // Packing
  for (int peer = 0; peer < npe (); peer++) {
    if (peer != pid ()) {
      int si = 0;
      foreach_ibscalar (slist) {
        int nn = ibmm.snd_migrate[peer].nodes.size;
        for (int ni = 0; ni < nn; ni++) {
          // bool ib_set_dirty = false;
          IBNode* node = ibmm.snd_migrate[peer].nodes.ptrs[ni];
          ibmm.snd_migrate[peer].buffs[si * nn + ni] = ibval (s);
        }
        si++;
      }
    }
  }

  // Exchange
  int maxreq = 2 * (npe () - 1);
  MPI_Request* requests =
    maxreq > 0 ? (MPI_Request*) malloc (maxreq * sizeof (MPI_Request)) : NULL;
  int nreq = 0;

  for (int peer = 0; peer < npe (); peer++) {
    if (peer == pid ())
      continue;
    double* buf = ibmm.rcv_migrate[peer].buffs;
    int count = (int) (ibmm.rcv_migrate[peer].nodes.size * nscalars);
    MPI_Irecv (
      buf, count, MPI_DOUBLE, peer, 0, MPI_COMM_WORLD, &requests[nreq++]);
  }

  for (int peer = 0; peer < npe (); peer++) {
    if (peer == pid ())
      continue;
    double* buf = ibmm.snd_migrate[peer].buffs;
    int count = (int) (ibmm.snd_migrate[peer].nodes.size * nscalars);
    MPI_Isend (
      buf, count, MPI_DOUBLE, peer, 0, MPI_COMM_WORLD, &requests[nreq++]);
  }

  if (nreq > 0)
    MPI_Waitall (nreq, requests, MPI_STATUSES_IGNORE);
  free (requests);

  // Unpacking
  for (int peer = 0; peer < npe (); peer++) {
    if (peer != pid ()) {
      int si = 0;
      foreach_ibscalar (slist) {
        int nn = ibmm.rcv_migrate[peer].nodes.size;
        for (int ni = 0; ni < nn; ni++) {
          // bool ib_set_dirty = false;
          IBNode* node = ibmm.rcv_migrate[peer].nodes.ptrs[ni];
          ibval (s) = ibmm.rcv_migrate[peer].buffs[si * nn + ni];
        }
        si++;
      }
    }
  }

  for (int peer = 0; peer < npe (); peer++) {
    if (peer != pid ()) {
      ibexchangelist_free_buffer (&ibmm.snd_migrate[peer]);
      ibexchangelist_free_buffer (&ibmm.rcv_migrate[peer]);
    }
  }

  ibmm.dirty = false;
}

/**
 * @brief Updates the list of local IBNodes as well as those with kernel support
 * into an ajacent process
 */
trace void ibmeshmanager_boundary (IBscalar* slist = iball) {
  ibmeshmanager_update_pid ();

  // IBscalar* slist = NULL;

  // foreach_ibscalar (list) {
  //   if (ibdirty (s))
  //     slist = iblist_add (slist, s);
  // }

  size_t nscalars = iblist_len (slist);

  if (!nscalars)
    return;

  //ibmeshmanager_check_exchange_counts ("boundary", ibmm.snd_boundary,
                                       ibmm.rcv_boundary);

  for (int peer = 0; peer < npe (); peer++) {
    if (peer != pid ()) {
      ibexchangelist_init_buffer (&ibmm.snd_boundary[peer], nscalars);
      ibexchangelist_init_buffer (&ibmm.rcv_boundary[peer], nscalars);
    }
  }

  // Packing
  for (int peer = 0; peer < npe (); peer++) {
    if (peer != pid ()) {
      int si = 0;
      foreach_ibscalar (slist) {
        int nn = ibmm.snd_boundary[peer].nodes.size;
        for (int ni = 0; ni < nn; ni++) {
          // bool ib_set_dirty = false;
          IBNode* node = ibmm.snd_boundary[peer].nodes.ptrs[ni];
          ibmm.snd_boundary[peer].buffs[si * nn + ni] = ibval (s);
        }
        si++;
      }
    }
  }

  // Exchange
  int maxreq = 2 * (npe () - 1);
  MPI_Request* requests =
    maxreq > 0 ? (MPI_Request*) malloc (maxreq * sizeof (MPI_Request)) : NULL;
  int nreq = 0;

  for (int peer = 0; peer < npe (); peer++) {
    if (peer == pid ())
      continue;
    double* buf = ibmm.rcv_boundary[peer].buffs;
    int count = (int) (ibmm.rcv_boundary[peer].nodes.size * nscalars);
    MPI_Irecv (
      buf, count, MPI_DOUBLE, peer, 0, MPI_COMM_WORLD, &requests[nreq++]);
  }

  for (int peer = 0; peer < npe (); peer++) {
    if (peer == pid ())
      continue;
    double* buf = ibmm.snd_boundary[peer].buffs;
    int count = (int) (ibmm.snd_boundary[peer].nodes.size * nscalars);
    MPI_Isend (
      buf, count, MPI_DOUBLE, peer, 0, MPI_COMM_WORLD, &requests[nreq++]);
  }

  if (nreq > 0)
    MPI_Waitall (nreq, requests, MPI_STATUSES_IGNORE);
  free (requests);

  // Unpacking
  for (int peer = 0; peer < npe (); peer++) {
    if (peer != pid ()) {
      int si = 0;
      foreach_ibscalar (slist) {
        int nn = ibmm.rcv_boundary[peer].nodes.size;
        for (int ni = 0; ni < nn; ni++) {
          // bool ib_set_dirty = false;
          IBNode* node = ibmm.rcv_boundary[peer].nodes.ptrs[ni];
          ibval (s) = ibmm.rcv_boundary[peer].buffs[si * nn + ni];
        }
        si++;
      }
    }
  }

  for (int peer = 0; peer < npe (); peer++) {
    if (peer != pid ()) {
      ibexchangelist_free_buffer (&ibmm.snd_boundary[peer]);
      ibexchangelist_free_buffer (&ibmm.rcv_boundary[peer]);
    }
  }
}

#endif // _MPI
