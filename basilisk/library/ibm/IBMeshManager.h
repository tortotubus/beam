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
  bool dirty;
} IBMeshManager;

/**
 * Globals
 */
IBMeshManager ibmm = {
  .meshes = NULL,
  .pool = {
    .pool = NULL,
    .active = {0},
    .len = 0, 
  }, 
  .nm = 0,
  .dirty = true
};


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
    IBMesh *mesh = &ibmm.meshes[mesh_id];
    NOT_UNUSED(mesh);
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
int ibmeshmanager_add_mesh();
void ibmeshmanager_delete_mesh (int mesh_id);
void ibmeshmanager_add_nodes (int mesh_id, int count);
void ibmeshmanager_delete_all_nodes (int mesh_id);
void ibmeshmanager_init_stencil_caches();
void ibmeshmanager_update_stencil_caches();

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
    foreach_ibmesh()
      ibmesh_free(mesh);
    free(ibmm.meshes);
  }

  ibmm.meshes = NULL;
  ibmm.nm = mesh_count;
  ibmm.meshes = (IBMesh*) calloc (mesh_count, sizeof (IBMesh));
  ibmm.pool = ibmempool_init (IBMESHMANAGER_POOL_SIZE_BYTES);

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

/**
 * @brief Initialize all node stencil caches
 * 
 * @relates IBMeshManager
 */
void ibmeshmanager_init_stencil_caches() {
  foreach_ibnode_per_ibmesh() {
    foreach_neighborhood_coord_level(node->lagpos, PESKIN_SUPPORT_RADIUS, mesh->refinement_level) {
      cache_append(&node->stencil, point, 0);
    }
  }
}

/**
 * @brief Update all node stencil caches
 * 
 * @relates IBMeshManages
 */
void ibmeshmanager_update_stencil_caches() {
  foreach_ibnode_per_ibmesh() {
    ibnode_free_stencil(node);
    ibnode_init_stencil(node);
  }
  ibmeshmanager_init_stencil_caches();
}

void ibmeshmanager_mpi_reduction() {

}