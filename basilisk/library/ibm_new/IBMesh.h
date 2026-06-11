

#include "ibm/IBMempool.h"
#include "ibm/IBNodeList.h"
#include "ibm/IBNode.h"
#include "ibm/IBMeshModel.h"

/*
 * Type Definitions
 */
typedef struct {
  int pid;
  int gid;
  int depth;
  int nnode_global;
  bool sparse_managed;
  IBNodeList nodes;
  IBMeshModel model;
} IBMesh;

/**
 * Function Declarations
 */

void ibmesh_init (IBMesh* mesh);
void ibmesh_free (IBMesh* mesh);
void ibmesh_add_nodes (IBMesh* mesh, IBMempool* pool, int count);
void ibmesh_delete_all_nodes (IBMesh* mesh, IBMempool* pool);
IBNode* ibmesh_find_node_by_lid (IBMesh* mesh, int node_lid);
IBNode* ibmesh_get_or_create_node_by_lid (IBMesh* mesh,
                                          IBMempool* pool,
                                          int node_lid);
void ibmesh_prune_untouched_nodes (IBMesh* mesh,
                                   IBMempool* pool,
                                   int epoch);

/**
 * Functions Definitions
 */

/**
 * @brief Initialize an IBMesh structure.
 *
 * @param mesh Pointer to the IBMesh to initialize. If NULL, does nothing.
 *
 * @details Initializes the internal IBNodeList with zero initial capacity.
 * After this, nodes may be manually added by the user, or set a model on the
 * mesh.
 *
 * @memberof IBMesh
 */
void ibmesh_init (IBMesh* mesh) {
  if (!mesh)
    return;

  mesh->pid = 0;
  mesh->gid = -1;
  mesh->depth = 0;
  mesh->nnode_global = 0;
  mesh->sparse_managed = false;
  mesh->model = ibmeshmodel_init();
  ibnodelist_init (&mesh->nodes, 0);
}

/**
 * @brief Free all resources associated with an IBMesh.
 *
 * @param mesh Pointer to the IBMesh to free. If NULL, does nothing.
 *
 * @details Deallocates the internal IBNodeList. Does not free the IBNode
 * objects themselves (they should be freed through the IBMempool).
 *
 * @memberof IBMesh
 */
void ibmesh_free (IBMesh* mesh) {
  if (!mesh)
    return;

  ibnodelist_free (&mesh->nodes);
  ibmeshmodel_destroy (&mesh->model);
}

/**
 * @brief Add multiple new nodes to the mesh.
 *
 * @param mesh Pointer to the IBMesh. Must not be NULL.
 * @param pool Pointer to the IBMempool to allocate nodes from. Must not be
 * NULL.
 * @param count Number of nodes to add. Must be > 0.
 *
 * @details Allocates the specified number of IBNode objects from the pool,
 * initializes each one, and adds them to the mesh's node list. If any
 * allocation or initialization fails for a particular node, that node is
 * skipped and the operation continues. Reserves capacity in the node list
 * before adding.
 *
 * @memberof IBMesh
 */
void ibmesh_add_nodes (IBMesh* mesh, IBMempool* pool, int count) {
  if (!mesh || !pool || count <= 0)
    return;

  ibnodelist_reserve (&mesh->nodes, mesh->nodes.size + (size_t) count);

  for (int node_idx = 0; node_idx < count; node_idx++) {
    IBNode* node = ibmempool_alloc_node (pool);
    if (!node)
      continue;
    if (ibnode_init (node) != 0) {
      ibmempool_free_node_ptr (pool, node);
      continue;
    }
    node->mesh_gid = mesh->gid;
    node->node_lid = node_idx;
    node->depth = mesh->depth;
    ibnodelist_push (&mesh->nodes, node);
  }
}

/**
 * @brief Delete all nodes from the mesh and return them to the pool.
 *
 * @param mesh Pointer to the IBMesh. Must not be NULL.
 * @param pool Pointer to the IBMempool to return freed nodes to. Must not be
 * NULL.
 *
 * @details Removes all nodes from the mesh using swap-with-last removal and
 * returns each node's memory to the pool. The mesh will be empty after this
 * call.
 *
 * @memberof IBMesh
 */
void ibmesh_delete_all_nodes (IBMesh* mesh, IBMempool* pool) {
  if (!mesh || !pool)
    return;

  while (mesh->nodes.size > 0) {
    IBNode* node = ibnodelist_remove_swap (&mesh->nodes, mesh->nodes.size - 1);
    if (node)
      ibmempool_free_node_ptr (pool, node);
  }
}

IBNode* ibmesh_find_node_by_lid (IBMesh* mesh, int node_lid) {
  if (!mesh || node_lid < 0)
    return NULL;

  for (size_t i = 0; i < mesh->nodes.size; i++) {
    IBNode* node = mesh->nodes.ptrs[i];
    if (node && node->node_lid == node_lid)
      return node;
  }
  return NULL;
}

IBNode* ibmesh_get_or_create_node_by_lid (IBMesh* mesh,
                                          IBMempool* pool,
                                          int node_lid) {
  if (!mesh || !pool || node_lid < 0)
    return NULL;

  IBNode* node = ibmesh_find_node_by_lid (mesh, node_lid);
  if (node)
    return node;

  node = ibmempool_alloc_node (pool);
  if (!node)
    return NULL;
  if (ibnode_init (node) != 0) {
    ibmempool_free_node_ptr (pool, node);
    return NULL;
  }

  node->mesh_gid = mesh->gid;
  node->node_lid = node_lid;
  node->depth = mesh->depth;
  ibnodelist_push (&mesh->nodes, node);
  return node;
}

void ibmesh_prune_untouched_nodes (IBMesh* mesh,
                                   IBMempool* pool,
                                   int epoch) {
  if (!mesh || !pool || !mesh->sparse_managed)
    return;

  for (size_t i = 0; i < mesh->nodes.size;) {
    IBNode* node = mesh->nodes.ptrs[i];
    if (node && node->touched_epoch == epoch) {
      i++;
      continue;
    }

    IBNode* dead = ibnodelist_remove_swap (&mesh->nodes, i);
    if (dead)
      ibmempool_free_node_ptr (pool, dead);
  }
}

/**
 * @brief Set the mesh model
 */
void ibmesh_set_model (IBMesh* mesh, IBMempool* pool, IBMeshModel model) {
  
  switch (model.type) {
  case IB_MODEL_VELOCITY_COUPLED: {
    mesh->model = model;
    ibmesh_delete_all_nodes (mesh, pool);
    int node_count = mesh->model.velocity_ops->node_count (model.ctx);
    mesh->nnode_global = node_count;
    if (!mesh->sparse_managed)
      ibmesh_add_nodes (mesh, pool, node_count);
    mesh->model.velocity_ops->sync (model.ctx, mesh);
    break;
  }
  case IB_MODEL_FORCE_COUPLED: {
    mesh->model = model;
    ibmesh_delete_all_nodes (mesh, pool);
    int node_count = mesh->model.force_ops->node_count (model.ctx);
    mesh->nnode_global = node_count;
    if (!mesh->sparse_managed)
      ibmesh_add_nodes (mesh, pool, node_count);
    mesh->model.force_ops->sync (model.ctx, mesh);
    break;
  }
  default:
    break;
  }
}
