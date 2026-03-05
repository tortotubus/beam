#include "library/ibm/IBNode.h"
#include "library/ibm/IBMempool.h"
#include "library/ibm/IBMesh.h"

/**
 * Type definitions
 */

typedef struct {
  IBMesh* meshes;
  IBMempool pool;
  int nm;
#if _MPI
  bool dirty;
  IBNodeList local;
  IBNodeList send;
  IBNodeList receive;
#endif
} IBMeshManager;

/**
 * Globals
 */
IBMeshManager ibmm = {0};

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
macro foreach_ibnode () {
  for (size_t node_id = 0; node_id < ibmm.pool.active.size; node_id++) {
    IBNode* node = ibmm.pool.active.ptrs[node_id];
    NOT_UNUSED (node);
    // clang-format off
    {...}
    // clang-format on
  }
}

/**
 * @brief Loops through all nodes of all meshes
 *
 * @relates IBMeshManager
 */
macro foreach_ibnode_per_ibmesh () {
  for (size_t mesh_id = 0; mesh_id < ibmm.nm; mesh_id++) {
    IBMesh* mesh = &ibmm.meshes[mesh_id];
    NOT_UNUSED (mesh);
    for (size_t node_id = 0; node_id < mesh->nodes.size; node_id++) {
      IBNode* node = mesh->nodes.ptrs[node_id];
      NOT_UNUSED (node);
      // clang-format off
      {...}
      // clang-format on
    }
  }
}

/* Function declarations */

void ibmeshmanager_init (int mesh_count);
void ibmeshmanager_free ();
int ibmeshmanager_add_mesh ();
void ibmeshmanager_delete_mesh (int mesh_id);
void ibmeshmanager_add_nodes (int mesh_id, int count);
void ibmeshmanager_delete_all_nodes (int mesh_id);

#if TREE
void ibmeshmanager_init_stencil_caches ();
void ibmeshmanager_update_stencil_caches ();
#endif

#if _MPI
void ibmeshmanager_update_pid ();
#endif

/* Function definitions */

/**
 * @brief Initialize the immersed boundary mesh manager
 *
 * @param nm The number of meshes you plan to have
 *
 * @relates IBMeshManager
 */
void ibmeshmanager_init (int mesh_count) {
  if (!ibmm.meshes) {
    foreach_ibmesh () ibmesh_free (mesh);
    free (ibmm.meshes);
  }

  ibmm.meshes = NULL;
  ibmm.nm = mesh_count;
  ibmm.meshes = (IBMesh*) calloc (mesh_count, sizeof (IBMesh));
  ibmm.pool = ibmempool_init (IBMESHMANAGER_POOL_SIZE_BYTES);

#if _MPI
  ibnodelist_init (&ibmm.local, 0);
  ibnodelist_init (&ibmm.send, 0);
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

#if TREE
/**
 * @brief Initialize all node stencil caches
 *
 * @relates IBMeshManager
 */
void ibmeshmanager_init_stencil_caches () {
  foreach_ibnode_per_ibmesh () {
    foreach_neighborhood_coord_level (
      node->lagpos, PESKIN_SUPPORT_RADIUS, mesh->refinement_level) {
      cache_append (&node->stencil, point, 0);
    }
  }
}

/**
 * @brief Update all node stencil caches
 *
 * @relates IBMeshManager
 */
void ibmeshmanager_update_stencil_caches () {
  foreach_ibnode_per_ibmesh () {
    ibnode_free_stencil (node);
    ibnode_init_stencil (node);
  }
  ibmeshmanager_init_stencil_caches ();
}

#endif

#if _MPI
/**
 * @brief Update all IBNode pids
 *
 * Locally queries the eulerian grid cell enclosed by a box
 * and sets the pid of the node to the pid of that cell
 *
 * @relates IBMeshManager
 */
void ibmeshmanager_update_pid () {
  foreach_ibnode_per_ibmesh () {
    coord centre_coord = node->lagpos;
    Point centre_point = locate_level (
      centre_coord.x, centre_coord.y, centre_coord.z, mesh->refinement_level);
    if (centre_point.level > 0) {
      Point point = {0};
      point = centre_point;
#if dimension == 1
      int ig = point.i;
#elif dimension == 2
      int ig = point.i;
      int jg = point.j;
#else
      int ig = point.i;
      int jg = point.j;
      int kg = point.j;
#endif
      POINT_VARIABLES ();
      node->pid = cell.pid;
    }
  }
}

/**
 * @brief Updates the list of local IBNodes as well as those with kernel support
 * into an ajacent process
 */
void ibmeshmanager_update_mpi_boundaries () {
  if (!ibmm.dirty) 
    return;

  ibmeshmanager_update_pid();

  ibnodelist_clear(&ibmm.local);
  ibnodelist_clear(&ibmm.receive);
  ibnodelist_clear(&ibmm.send);

  foreach_ibnode_per_ibmesh() {
    if (node->pid == pid()) {
      ibnodelist_push(&ibmm.local, node);
      foreach_cache(node->stencil) {
        if (cell.pid != pid()) {
          ibnodelist_push(&ibmm.send, node);
        }
      }
    } else {
      foreach_cache(node->stencil) {
        if (cell.pid == pid()) {
          ibnodelist_push(&ibmm.receive, node);
        }
      }
    }
  }
}

#endif

// Print Basilisk MPI neighbor sets (send/receive pid lists).
// Works with grid/tree-mpi.h (MPI tree) where mpi_boundary is the active
// boundary.
//
// Usage examples:
//   event init (t = 0) { mpi_boundary_update_buffers(); mpi_print_neighbors(-1,
//   stderr); } event adapt (i++)  { /* after adapt/balance */
//   mpi_boundary_update_buffers(); mpi_print_neighbors(pid(), stderr); }
