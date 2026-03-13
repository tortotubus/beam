#include "library/ibm/IBNode.h"
#include "library/ibm/IBMempool.h"
#include "library/ibm/IBMesh.h"
#include "library/ibm/IBMeshManagerMPI.h"

#ifndef MPI_AUTO_BOUNDARY
#define MPI_AUTO_BOUNDARY 0
#endif

/**
 * Type definitions
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

/**
 * Globals
 */
IBMeshManager ibmm = {0};

/* Function declarations */

void ibmeshmanager_init (int mesh_count);
void ibmeshmanager_free ();
int ibmeshmanager_add_mesh ();
void ibmeshmanager_delete_mesh (int mesh_id);
void ibmeshmanager_add_nodes (int mesh_id, int count);
void ibmeshmanager_delete_all_nodes (int mesh_id);

#if _MPI
int _ibmeshmanager_get_pid (Point p);
void ibmeshmanager_update_pid ();
void ibmeshmanager_boundary (IBscalar* list = iball);
#endif

/**
 * Macros
 */

#define IBMESHMANAGER_POOL_SIZE_BYTES (1 << 19)

/**
 * @brief Loops through all meshes
 *
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
 * @brief Loops through all nodes in the node pool
 *
 * @relates IBMeshManager
 */
macro foreach_ibnode (bool local_only = false) {
#if _MPI
  {
    IBNodeList nlist = local_only ? ibmm.local : ibmm.pool.active;
    for (size_t node_id = 0; node_id < nlist.size; node_id++) {
      IBNode* node = nlist.ptrs[node_id];
      NOT_UNUSED (node);
      bool ib_set_dirty = false;
      NOT_UNUSED (ib_set_dirty);
      // clang-format off
      {...}
      // clang-format on
    }
    if (MPI_AUTO_BOUNDARY)
      ibmeshmanager_boundary ();
  }
#else
  for (size_t node_id = 0; node_id < ibmm.pool.active.size; node_id++) {
    IBNode* node = ibmm.pool.active.ptrs[node_id];
    NOT_UNUSED (node);
    bool ib_set_dirty = false;
    NOT_UNUSED (ib_set_dirty);
    // clang-format off
    {...}
    // clang-format on
  }
#endif
}

/**
 * @brief Loops through all nodes of all meshes
 *
 * @relates IBMeshManager
 */
macro foreach_ibnode_per_ibmesh (bool boundary = MPI_AUTO_BOUNDARY) {
  for (size_t mesh_id = 0; mesh_id < ibmm.nm; mesh_id++) {
    IBMesh* mesh = &ibmm.meshes[mesh_id];
    NOT_UNUSED (mesh);
    for (size_t node_id = 0; node_id < mesh->nodes.size; node_id++) {
      IBNode* node = mesh->nodes.ptrs[node_id];
      NOT_UNUSED (node);
      bool ib_set_dirty = false;
      NOT_UNUSED (ib_set_dirty);
      // clang-format off
      {...}
      // clang-format on
    }
  }
#if _MPI
  // if (MPI_AUTO_BOUNDARY)
  // ibmeshmanager_boundary ();
#endif
}

/* Function definitions */

/**
 * @brief Initialize the immersed boundary mesh manager
 *
 * @param nm The number of meshes you plan to have
 *
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

#endif
#endif

  foreach_ibmesh () {
    ibmesh_init (mesh);
  }
}

/**
 * @brief Free all members in the immersed boundary mesh manager.
 *
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
 *
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
 *
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
 * @brief Deletes the mesh from the manager, freeing the mesh object, and
 * marking as free its nodes in the pool.
 *
 * @relates IBMeshManager
 */
// void ibmeshmanager_delete_all_meshes () {
//   IBMempool* pool = &ibmm.pool;
//   foreach_ibmesh (ibmm) {
//     ibmesh_free (mesh);
//   }
//   free (ibmm.meshes);
//   ibmempool_free (&ibmm.pool);
//   ibmm.meshes = NULL;
//   ibmm.nm = 0;
// }

/**
 * @brief Bulk adds nodes to a given mesh
 *
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
 *
 * @relates IBMeshManager
 */
void ibmeshmanager_delete_all_nodes (int mesh_id) {
  IBMesh* mesh = &ibmm.meshes[mesh_id];
  IBMempool* pool = &ibmm.pool;
  ibmesh_delete_all_nodes (mesh, pool);
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
  int* node_is_local = (int*) calloc ((size_t) nnode, sizeof (int));
  assert (node_pids && node_is_local);
  for (int i = 0; i < nnode; i++)
    node_pids[i] = -1;

  /* Pass 1: local owner candidates only. */
  foreach_ibnode () {
    coord d = node->pos;
    coord_periodic_boundary (d);
#if TREE
    Point point = locate_level (d.x, d.y, d.z, node->depth);
    int ig = 0, jg = 0, kg = 0;
    NOT_UNUSED (ig);
    NOT_UNUSED (jg);
    NOT_UNUSED (kg);
    POINT_VARIABLES ();

    if (point.level >= 0) {
      if (allocated (0)) {
        if (is_local (cell)) {
          node_pids[node_id] = pid ();
          node_is_local[node_id] = 1;
        }
      }
    }
#else
    Point point = locate (d.x, d.y, d.z);
    if (point.level >= 0) {
      if (!is_boundary (point)) {
        node_pids[node_id] = pid ();
        node_is_local[node_id] = 1;
      }
    }
#endif
  }

  // return;

  /* Synchronize globally: each node gets one agreed owner pid. */
  mpi_all_reduce_array (node_pids, MPI_INT, MPI_MAX, nnode);

  // return;

  /* Pass 2: build local/snd_boundary/rcv_boundary from resolved owners. */
  foreach_ibnode () {

    // Check if the owner has actually changed
    int old_pid = node->pid;
    int new_pid = node_pids[node_id];
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

#if TREE
    int ig = 0, jg = 0, kg = 0;
    NOT_UNUSED (ig);
    NOT_UNUSED (jg);
    NOT_UNUSED (kg);
    coord d = node->pos;
    coord_periodic_boundary (d);
    Point point = locate_level (d.x, d.y, d.z, node->depth);
    POINT_VARIABLES ();

    if (node->pid == pid ()) { // local node
      ibnodelist_push (&ibmm.local, node);
      if (point.level >= 0 && allocated (0)) {
        foreach_neighbor (PESKIN_SUPPORT_RADIUS) {
          if (allocated (0) && !is_local (cell))
            ibexchangelist_push_unique (&ibmm.snd_boundary[cell.pid], node);
        }
      }
    } // local node
    else { // remote node
      bool has_local_support = false;
      if (point.level >= 0 && allocated (0)) {
        foreach_neighbor (PESKIN_SUPPORT_RADIUS) {
          if (!has_local_support && allocated (0) && is_local (cell)) {
            has_local_support = true;
          }
        }
      }
      if (has_local_support)
        ibexchangelist_push_unique (&ibmm.rcv_boundary[node->pid], node);
    } // remote node
#else // !TREE
    coord d = node->pos;
    coord_periodic_boundary (d);
    Point point = locate_nonlocal (d.x, d.y, d.z);


    int ig = 0, jg = 0, kg = 0;
    NOT_UNUSED (ig); NOT_UNUSED (jg); NOT_UNUSED (kg);

    if (point.level >= 0) {
      if (node->pid == pid ()) { // local node
        foreach_neighbor(PESKIN_SUPPORT_RADIUS) {
          if (is_boundary (point)) {
            int point_pid = _ibmeshmanager_get_pid (point);
            // printf("%d\n", point_pid);
            if (point_pid >= 0) // MPI boundary, not domain boundary
              ibexchangelist_push_unique (&ibmm.snd_boundary[point_pid], node);
          }
        }
      } else { // remote node
        ibexchangelist_push_unique (&ibmm.rcv_boundary[node->pid], node);
      }
    } // remote node
#endif
  }

  free (node_is_local);
  free (node_pids);

  // Migrate nodes
  IBscalar* slist = iball;
  size_t nscalars = iblist_len (slist);

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
          bool ib_set_dirty = false;
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
          bool ib_set_dirty = false;
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
void ibmeshmanager_boundary (IBscalar* slist = iball) {
  if (ibmm.dirty)
    ibmeshmanager_update_pid ();

  // IBscalar* slist = NULL;

  // foreach_ibscalar (list) {
  //   if (ibdirty (s))
  //     slist = iblist_add (slist, s);
  // }

  size_t nscalars = iblist_len (slist);

  if (!nscalars)
    return;

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
          bool ib_set_dirty = false;
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
          bool ib_set_dirty = false;
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
