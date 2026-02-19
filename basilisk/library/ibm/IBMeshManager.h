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
  .nm = 0
};


/**
 * Macros
 */

#define IBMESHMANAGER_POOL_SIZE_BYTES (1 << 19)

/**
 * @brief Loops through all meshes
 *
 * @memberof IBMeshManager
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
 * @memberof IBMeshManager
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
 * @memberof IBMeshManager
 */
macro foreach_ibnode_per_ibmesh () {
  for (size_t mesh_id = 0; mesh_id < ibmm.nm; mesh_id++) {
    IBMesh *mesh = &ibmm.meshes[mesh_id];
    NOT_UNUSED(mesh);
    for (size_t node_id = 0; node_id < ibmm.pool.active.size; node_id++) {
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
// void ibmeshmanager_delete_all_meshes();

void ibmeshmanager_add_nodes (int mesh_id, int count);
void ibmeshmanager_delete_all_nodes (int mesh_id);

/* Function definitions */

/**
 * @brief Initialize the immersed boundary mesh manager
 *
 * @param nm The number of meshes you plan to have
 *
 * @memberof IBMeshManager
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
 * @memberof IBMeshManager
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
 * @memberof IBMeshManager
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
 * @memberof IBMeshManager
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
 * @memberof IBMeshManager
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
 * @memberof IBMeshManager
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
 * @memberof IBMeshManager
 */
void ibmeshmanager_delete_all_nodes (int mesh_id) {
  IBMesh* mesh = &ibmm.meshes[mesh_id];
  IBMempool* pool = &ibmm.pool;
  ibmesh_delete_all_nodes (mesh, pool);
}
