

#include "ibm/IBNode.h"
#include "ibm/IBMempool.h"
#include "ibm/IBMesh.h"
#include "ibm/IBMeshManagerMPI.h"

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
static bool ibmeshmanager_initialized = false;
static int ibmeshmanager_sparse_epoch = 0;

// ============================================================================
// Function declarations
// ============================================================================

void ibmeshmanager_init (int mesh_count);
void ibmeshmanager_free ();
int ibmeshmanager_add_mesh ();
void ibmeshmanager_delete_mesh (int mesh_id);
void ibmeshmanager_add_nodes (int mesh_id, int count);
void ibmeshmanager_delete_all_nodes (int mesh_id);
void ibmeshmanager_set_model (int mesh_id, IBMeshModel model, int owner = 0);
void ibmeshmanager_set_model_sparse (int mesh_id,
                                     int mesh_gid,
                                     IBMeshModel model,
                                     int owner = 0);
void ibmeshmanager_advance_positions (double dt);
void ibmeshmanager_advance_force_coupled_positions (double dt);
void ibmeshmanager_evaluate_velocity_coupled_midpoints (double dt);
void ibmeshmanager_advance_velocity_coupled_positions (double dt);
void ibmeshmanager_sync_model_outputs ();
void ibmeshmanager_sync_velocity_coupled_model_outputs ();
void ibmeshmanager_sync_force_coupled_model_outputs ();

#if _MPI
int _ibmeshmanager_get_pid (Point p);
void ibmeshmanager_update_pid ();
void ibmeshmanager_boundary (IBscalar* list = iball);
#endif

// ============================================================================
// Macros
// ============================================================================

#define IBMESHMANAGER_POOL_SIZE_BYTES (1 << 19)

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
  if (ibmeshmanager_initialized)
    return;

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
  ibmeshmanager_initialized = true;
}

/**
 * @brief Free all members in the immersed boundary mesh manager.
 * @relates IBMeshManager
 */
void ibmeshmanager_free () {
  if (!ibmeshmanager_initialized)
    return;

  IBMempool* pool = &ibmm.pool;
  foreach_ibmesh () {
    ibmesh_free (mesh);
  }
  free (ibmm.meshes);
  ibmempool_free (&ibmm.pool);
  ibmm.meshes = NULL;
  ibmm.nm = 0;
  ibmeshmanager_initialized = false;

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
static void ibmeshmanager_set_model_internal (int mesh_id,
                                              IBMeshModel model,
                                              int owner,
                                              bool sparse_managed_enabled,
                                              int mesh_gid) {
  IBMesh* mesh = &ibmm.meshes[mesh_id];
  IBMempool* pool = &ibmm.pool;
  mesh->pid = owner;
  mesh->gid = mesh_gid;
  mesh->sparse_managed = sparse_managed_enabled;
  int node_count = 0;
  IBscalar* input_fields = NULL;
  IBscalar* output_fields = NULL;

  switch (model.type) {
  case IB_MODEL_VELOCITY_COUPLED:
    assert (model.velocity_ops && model.velocity_ops->node_count);
    node_count = model.velocity_ops->node_count (model.ctx);
    if (model.velocity_ops->input_fields)
      input_fields = model.velocity_ops->input_fields (model.ctx);
    if (model.velocity_ops->output_fields)
      output_fields = model.velocity_ops->output_fields (model.ctx);
    break;
  case IB_MODEL_FORCE_COUPLED:
    assert (model.force_ops && model.force_ops->node_count);
    node_count = model.force_ops->node_count (model.ctx);
    if (model.force_ops->input_fields)
      input_fields = model.force_ops->input_fields (model.ctx);
    if (model.force_ops->output_fields)
      output_fields = model.force_ops->output_fields (model.ctx);
    break;
  default:
    break;
  }

#if _MPI
  if (owner != pid ()) {
    IBMeshModel passive = ibmeshmodel_passive_from (
      model.type, node_count, input_fields, output_fields);
    ibmeshmodel_destroy (&model);
    model = passive;
  }
#endif

  mesh->nnode_global = node_count;
  free (input_fields);
  free (output_fields);
  ibmesh_set_model (mesh, pool, model);
  mesh->pid = owner;
  mesh->gid = mesh_gid;
  mesh->nnode_global = node_count;
}

void ibmeshmanager_set_model (int mesh_id, IBMeshModel model, int owner = 0) {
  ibmeshmanager_set_model_internal (mesh_id, model, owner, false, mesh_id);
}

void ibmeshmanager_set_model_sparse (int mesh_id,
                                     int mesh_gid,
                                     IBMeshModel model,
                                     int owner = 0) {
#if _MPI
  ibmeshmanager_set_model_internal (mesh_id, model, owner, true, mesh_gid);
#else
  ibmeshmanager_set_model_internal (mesh_id, model, owner, false, mesh_gid);
#endif
}

#if _MPI
static IBscalar* ibmeshmanager_add_vector_fields (IBscalar* fields,
                                                  IBvector v) {
  foreach_dimension ()
    fields = iblist_add (fields, v.x);
  return fields;
}

static IBscalar* ibmeshmanager_default_input_fields (int model_type) {
  IBscalar* fields = NULL;
  switch (model_type) {
  case IB_MODEL_VELOCITY_COUPLED:
    fields = ibmeshmanager_add_vector_fields (fields, nvel);
    fields = ibmeshmanager_add_vector_fields (fields, nforce);
    break;
  default:
    break;
  }
  return fields;
}

static IBscalar* ibmeshmanager_default_output_fields (int model_type) {
  IBscalar* fields = NULL;
  switch (model_type) {
  case IB_MODEL_VELOCITY_COUPLED:
    fields = ibmeshmanager_add_vector_fields (fields, npos);
    fields = ibmeshmanager_add_vector_fields (fields, nvel);
    fields = ibmeshmanager_add_vector_fields (fields, nforce);
    fields = iblist_add (fields, nweight);
    break;
  case IB_MODEL_FORCE_COUPLED:
    fields = ibmeshmanager_add_vector_fields (fields, npos);
    fields = ibmeshmanager_add_vector_fields (fields, nvel);
    break;
  default:
    break;
  }
  return fields;
}

static IBscalar* ibmeshmanager_model_input_fields (IBMesh* mesh) {
  if (!mesh)
    return NULL;

  IBscalar* fields = NULL;
  switch (mesh->model.type) {
  case IB_MODEL_VELOCITY_COUPLED:
    if (mesh->model.velocity_ops && mesh->model.velocity_ops->input_fields)
      fields = mesh->model.velocity_ops->input_fields (mesh->model.ctx);
    break;
  case IB_MODEL_FORCE_COUPLED:
    if (mesh->model.force_ops && mesh->model.force_ops->input_fields)
      fields = mesh->model.force_ops->input_fields (mesh->model.ctx);
    break;
  default:
    break;
  }

  return fields ? fields : ibmeshmanager_default_input_fields (mesh->model.type);
}

static IBscalar* ibmeshmanager_model_output_fields (IBMesh* mesh) {
  if (!mesh)
    return NULL;

  IBscalar* fields = NULL;
  switch (mesh->model.type) {
  case IB_MODEL_VELOCITY_COUPLED:
    if (mesh->model.velocity_ops && mesh->model.velocity_ops->output_fields)
      fields = mesh->model.velocity_ops->output_fields (mesh->model.ctx);
    break;
  case IB_MODEL_FORCE_COUPLED:
    if (mesh->model.force_ops && mesh->model.force_ops->output_fields)
      fields = mesh->model.force_ops->output_fields (mesh->model.ctx);
    break;
  default:
    break;
  }

  return fields ? fields : ibmeshmanager_default_output_fields (mesh->model.type);
}

static int ibmeshmanager_field_index (IBscalar* fields, IBscalar target) {
  int index = 0;
  foreach_ibscalar (fields) {
    if (s.i == target.i)
      return index;
    index++;
  }
  return -1;
}

static bool ibmeshmanager_has_local_kernel_support (coord pos, int depth) {
  bool has_support = false;
  foreach_neighbor_coord_level (PESKIN_SUPPORT_RADIUS, depth, pos) {
#if TREE
    if (!has_support && allocated (0) && is_local (cell))
      has_support = true;
#else
    has_support = true;
#endif
  }
  return has_support;
}

static bool ibmeshmanager_node_center_is_local (coord pos, int node_depth) {
  coord d = pos;
  coord_periodic_boundary (d);
#if TREE
  int effective_depth = node_depth > depth () ? depth () : node_depth;
  Point point = locate_level (d.x, d.y, d.z, effective_depth);
  int ig = 0, jg = 0, kg = 0;
  NOT_UNUSED (ig);
  NOT_UNUSED (jg);
  NOT_UNUSED (kg);
  POINT_VARIABLES ();
  return point.level >= 0 && allocated (0) && is_local (cell);
#else
  Point point = locate (d.x, d.y, d.z);
  return point.level >= 0 && !is_boundary (point);
#endif
}

#if _MPI
static IBNode* ibmeshmanager_temp_node_alloc (void) {
  IBNode* node = (IBNode*) calloc (1, ibmempool_stride (&ibmm.pool));
  assert (node);
  ibnode_init (node);
  return node;
}

static void ibmeshmanager_temp_node_free (IBNode* node) {
  free (node);
}

static void ibmeshmanager_reduce_fields_to_mesh_owner_sparse (IBMesh* mesh,
                                                              IBscalar* fields) {
  if (!mesh || !fields)
    return;

  const int nscalars = iblist_len (fields);
  const int nn = mesh->nnode_global;
  if (nscalars <= 0 || nn <= 0)
    return;

  const int owner = mesh->pid >= 0 ? mesh->pid : 0;
  double* send_values =
    (double*) calloc ((size_t) nn * (size_t) nscalars, sizeof (double));
  double* recv_values =
    owner == pid ()
      ? (double*) calloc ((size_t) nn * (size_t) nscalars, sizeof (double))
      : NULL;
  assert (send_values && (owner != pid () || recv_values));

  for (size_t ni = 0; ni < mesh->nodes.size; ni++) {
    IBNode* node = mesh->nodes.ptrs[ni];
    if (!node || node->node_lid < 0 || node->node_lid >= nn)
      continue;

    int si = 0;
    foreach_ibscalar (fields)
      send_values[(size_t) node->node_lid * (size_t) nscalars +
                  (size_t) si++] += ibval (s);
  }

  MPI_Reduce (send_values,
              recv_values,
              nn * nscalars,
              MPI_DOUBLE,
              MPI_SUM,
              owner,
              MPI_COMM_WORLD);

  if (owner == pid ()) {
    IBVelocityCoupledModelOps* ops = mesh->model.velocity_ops;
    bool handled_by_model =
      mesh->model.type == IB_MODEL_VELOCITY_COUPLED && ops && ops->input_node;

    for (int node_lid = 0; node_lid < nn; node_lid++) {
      IBNode* node = handled_by_model ? ibmeshmanager_temp_node_alloc ()
                                      : ibmesh_find_node_by_lid (mesh, node_lid);
      if (!node)
        continue;

      node->mesh_gid = mesh->gid;
      node->node_lid = node_lid;
      node->depth = mesh->depth;

      int si = 0;
      foreach_ibscalar (fields)
        ibval (s) =
          recv_values[(size_t) node_lid * (size_t) nscalars + (size_t) si++];

      if (handled_by_model) {
        ops->input_node (mesh->model.ctx, mesh, node_lid, node);
        ibmeshmanager_temp_node_free (node);
      }
    }
  }

  free (send_values);
  free (recv_values);
}

static void ibmeshmanager_broadcast_fields_from_mesh_owner_sparse (
  IBMesh* mesh, IBscalar* fields) {
  if (!mesh || !fields)
    return;

  const int nscalars = iblist_len (fields);
  const int nn = mesh->nnode_global;
  if (nscalars <= 0 || nn <= 0)
    return;

  const int owner = mesh->pid >= 0 ? mesh->pid : 0;
  int pos_index[3] = {-1, -1, -1};
  pos_index[0] = ibmeshmanager_field_index (fields, npos.x);
#if dimension >= 2
  pos_index[1] = ibmeshmanager_field_index (fields, npos.y);
#endif
#if dimension >= 3
  pos_index[2] = ibmeshmanager_field_index (fields, npos.z);
#endif

  bool missing_position = pos_index[0] < 0
#if dimension >= 2
                          || pos_index[1] < 0
#endif
#if dimension >= 3
                          || pos_index[2] < 0
#endif
    ;
  if (missing_position) {
    fprintf (stderr,
             "ibmeshmanager: sparse mesh %d output fields must include npos.\n",
             mesh->gid);
    assert (false);
  }

  int* depths = (int*) malloc ((size_t) nn * sizeof (int));
  double* values =
    (double*) calloc ((size_t) nn * (size_t) nscalars, sizeof (double));
  assert (depths && values);
  for (int i = 0; i < nn; i++)
    depths[i] = mesh->depth;

  if (owner == pid ()) {
    IBVelocityCoupledModelOps* ops = mesh->model.velocity_ops;
    bool handled_by_model =
      mesh->model.type == IB_MODEL_VELOCITY_COUPLED && ops && ops->sync_node;

    if (handled_by_model) {
      for (int node_lid = 0; node_lid < nn; node_lid++) {
        IBNode* node = ibmeshmanager_temp_node_alloc ();
        node->mesh_gid = mesh->gid;
        node->node_lid = node_lid;
        node->depth = mesh->depth;

        ops->sync_node (mesh->model.ctx, mesh, node_lid, node);

        depths[node_lid] = node->depth;
        int si = 0;
        foreach_ibscalar (fields)
          values[(size_t) node_lid * (size_t) nscalars + (size_t) si++] =
            ibval (s);

        ibmeshmanager_temp_node_free (node);
      }
    } else {
      for (size_t ni = 0; ni < mesh->nodes.size; ni++) {
        IBNode* node = mesh->nodes.ptrs[ni];
        if (!node || node->node_lid < 0 || node->node_lid >= nn)
          continue;

        depths[node->node_lid] = node->depth;
        int si = 0;
        foreach_ibscalar (fields)
          values[(size_t) node->node_lid * (size_t) nscalars + (size_t) si++] =
            ibval (s);
      }
    }
  }

  MPI_Bcast (depths, nn, MPI_INT, owner, MPI_COMM_WORLD);
  MPI_Bcast (values, nn * nscalars, MPI_DOUBLE, owner, MPI_COMM_WORLD);

  if (mesh->sparse_managed) {
    const int epoch = ++ibmeshmanager_sparse_epoch;
    for (int node_lid = 0; node_lid < nn; node_lid++) {
      coord pos = {0};
      pos.x =
        values[(size_t) node_lid * (size_t) nscalars + (size_t) pos_index[0]];
#if dimension >= 2
      pos.y =
        values[(size_t) node_lid * (size_t) nscalars + (size_t) pos_index[1]];
#endif
#if dimension >= 3
      pos.z =
        values[(size_t) node_lid * (size_t) nscalars + (size_t) pos_index[2]];
#endif
      bool has_support =
        ibmeshmanager_has_local_kernel_support (pos, depths[node_lid]);
      bool owns_center =
        ibmeshmanager_node_center_is_local (pos, depths[node_lid]);
      if (!has_support && !owns_center)
        continue;

      IBNode* node = ibmesh_get_or_create_node_by_lid (
        mesh, &ibmm.pool, node_lid);
      assert (node);
      node->mesh_gid = mesh->gid;
      node->node_lid = node_lid;
      node->depth = depths[node_lid];
      node->touched_epoch = epoch;

      int si = 0;
      foreach_ibscalar (fields)
        ibval (s) =
          values[(size_t) node_lid * (size_t) nscalars + (size_t) si++];
    }

    ibmesh_prune_untouched_nodes (mesh, &ibmm.pool, epoch);
  }

  free (depths);
  free (values);
}
#endif

static void ibmeshmanager_reduce_fields_to_mesh_owner (IBMesh* mesh,
                                                       IBscalar* fields) {
  if (!mesh || !fields)
    return;

  if (mesh->sparse_managed) {
    ibmeshmanager_reduce_fields_to_mesh_owner_sparse (mesh, fields);
    return;
  }

  const int nscalars = iblist_len (fields);
  const int nn = (int) mesh->nodes.size;
  if (nscalars <= 0 || nn <= 0)
    return;

  double* values = (double*) calloc ((size_t) nn*nscalars, sizeof (double));
  assert (values);

  for (int ni = 0; ni < nn; ni++) {
    IBNode* node = mesh->nodes.ptrs[ni];
    if (node->pid == pid ()) {
      int si = 0;
      foreach_ibscalar (fields)
        values[ni*nscalars + si++] = ibval (s);
    }
  }

  MPI_Allreduce (
    MPI_IN_PLACE, values, nn*nscalars, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  const int owner = mesh->pid >= 0 ? mesh->pid : 0;
  if (owner == pid ()) {
    for (int ni = 0; ni < nn; ni++) {
      IBNode* node = mesh->nodes.ptrs[ni];
      int si = 0;
      foreach_ibscalar (fields)
        ibval (s) = values[ni*nscalars + si++];
    }
  }

  free (values);
}

static void ibmeshmanager_broadcast_fields_from_mesh_owner (IBMesh* mesh,
                                                            IBscalar* fields) {
  if (!mesh || !fields)
    return;

  if (mesh->sparse_managed) {
    ibmeshmanager_broadcast_fields_from_mesh_owner_sparse (mesh, fields);
    return;
  }

  const int nscalars = iblist_len (fields);
  const int nn = (int) mesh->nodes.size;
  if (nscalars <= 0 || nn <= 0)
    return;

  const int owner = mesh->pid >= 0 ? mesh->pid : 0;
  double* values = (double*) malloc ((size_t) nn*nscalars*sizeof (double));
  assert (values);

  if (owner == pid ()) {
    for (int ni = 0; ni < nn; ni++) {
      IBNode* node = mesh->nodes.ptrs[ni];
      int si = 0;
      foreach_ibscalar (fields)
        values[ni*nscalars + si++] = ibval (s);
    }
  }

  MPI_Bcast (values, nn*nscalars, MPI_DOUBLE, owner, MPI_COMM_WORLD);

  if (owner != pid ()) {
    for (int ni = 0; ni < nn; ni++) {
      IBNode* node = mesh->nodes.ptrs[ni];
      int si = 0;
      foreach_ibscalar (fields)
        ibval (s) = values[ni*nscalars + si++];
    }
  }

  free (values);
}
#endif

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
  foreach_ibmesh () {
    if (model_type != IB_MODEL_INVALID && mesh->model.type != model_type)
      continue;
    IBscalar* fields = ibmeshmanager_model_input_fields (mesh);
    ibmeshmanager_reduce_fields_to_mesh_owner (mesh, fields);
    free (fields);
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
    IBscalar* fields = ibmeshmanager_model_output_fields (mesh);
    ibmeshmanager_broadcast_fields_from_mesh_owner (mesh, fields);
    free (fields);
  }

  ibmm.dirty = true;
  ibmeshmanager_update_pid ();
#endif
}

trace void ibmeshmanager_sync_model_outputs_filtered (int model_type) {
#if _MPI
  ibmeshmanager_update_pid ();
#endif

  foreach_ibmesh () {
    if (model_type != IB_MODEL_INVALID && mesh->model.type != model_type)
      continue;

    const int owner = mesh->pid >= 0 ? mesh->pid : 0;
    if (owner != pid ())
      continue;

    switch (mesh->model.type) {
    case IB_MODEL_VELOCITY_COUPLED:
      if (mesh->model.velocity_ops && mesh->model.velocity_ops->sync)
        mesh->model.velocity_ops->sync (mesh->model.ctx, mesh);
      break;
    case IB_MODEL_FORCE_COUPLED:
      if (mesh->model.force_ops && mesh->model.force_ops->sync)
        mesh->model.force_ops->sync (mesh->model.ctx, mesh);
      break;
    default:
      break;
    }
  }

#if _MPI
  foreach_ibmesh () {
    if (model_type != IB_MODEL_INVALID && mesh->model.type != model_type)
      continue;

    IBscalar* fields = ibmeshmanager_model_output_fields (mesh);
    ibmeshmanager_broadcast_fields_from_mesh_owner (mesh, fields);
    free (fields);
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

trace void ibmeshmanager_sync_model_outputs () {
  ibmeshmanager_sync_model_outputs_filtered (IB_MODEL_INVALID);
}

trace void ibmeshmanager_sync_velocity_coupled_model_outputs () {
  ibmeshmanager_sync_model_outputs_filtered (IB_MODEL_VELOCITY_COUPLED);
}

trace void ibmeshmanager_sync_force_coupled_model_outputs () {
  ibmeshmanager_sync_model_outputs_filtered (IB_MODEL_FORCE_COUPLED);
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

  const int np = npe ();
  ibnodelist_clear (&ibmm.local);

  for (int peer = 0; peer < np; peer++) {
    ibexchangelist_clear (&ibmm.rcv_boundary[peer]);
    ibexchangelist_clear (&ibmm.snd_boundary[peer]);
    ibexchangelist_clear (&ibmm.rcv_migrate[peer]);
    ibexchangelist_clear (&ibmm.snd_migrate[peer]);
  }

  foreach_ibmesh () {
    if (mesh->sparse_managed) {
      for (size_t node_id = 0; node_id < mesh->nodes.size; node_id++) {
        IBNode* node = mesh->nodes.ptrs[node_id];
        if (!node)
          continue;

        coord d = {0};
        foreach_dimension()
          d.x = ibval(npos.x);
        coord_periodic_boundary(d);

        node->pid = -1;
#if TREE
        int effective_depth = node->depth > depth() ? depth() : node->depth;
        Point point = locate_level(d.x, d.y, d.z, effective_depth);
        int ig = 0, jg = 0, kg = 0;
        NOT_UNUSED(ig);
        NOT_UNUSED(jg);
        NOT_UNUSED(kg);
        POINT_VARIABLES();

        if (point.level >= 0 && allocated(0) && is_local(cell))
          node->pid = pid();
#else
        Point point = locate(d.x, d.y, d.z);
        if (point.level >= 0 && !is_boundary(point))
          node->pid = pid();
#endif
        if (ibmeshmanager_has_local_kernel_support(d, node->depth))
          ibnodelist_push(&ibmm.local, node);
      }
      continue;
    }

    const int nnode = (int) mesh->nodes.size;
    if (nnode <= 0)
      continue;

    int* node_pids = (int*) malloc((size_t) nnode*sizeof(int));
    assert(node_pids);
    for (int i = 0; i < nnode; i++)
      node_pids[i] = -1;

    for (int ni = 0; ni < nnode; ni++) {
      IBNode* node = mesh->nodes.ptrs[ni];
      if (!node)
        continue;

      coord d = {0};
      foreach_dimension()
        d.x = ibval(npos.x);
      coord_periodic_boundary(d);
#if TREE
      int effective_depth = node->depth > depth() ? depth() : node->depth;
      Point point = locate_level(d.x, d.y, d.z, effective_depth);
      int ig = 0, jg = 0, kg = 0;
      NOT_UNUSED(ig);
      NOT_UNUSED(jg);
      NOT_UNUSED(kg);
      POINT_VARIABLES();

      if (point.level >= 0 && allocated(0) && is_local(cell))
        node_pids[ni] = pid();
#else
      Point point = locate(d.x, d.y, d.z);
      if (point.level >= 0 && !is_boundary(point))
        node_pids[ni] = pid();
#endif
    }

    mpi_all_reduce_array(node_pids, MPI_INT, MPI_MAX, nnode);

    for (int ni = 0; ni < nnode; ni++) {
      IBNode* node = mesh->nodes.ptrs[ni];
      if (!node)
        continue;

      int old_pid = node->pid;
      int new_pid = node_pids[ni];
      node->pid = new_pid;

      if (old_pid != -1 && old_pid != new_pid) {
        if (old_pid == pid() && new_pid >= 0 && new_pid < np)
          ibexchangelist_push(&ibmm.snd_migrate[new_pid], node);
        if (new_pid == pid() && old_pid >= 0 && old_pid < np)
          ibexchangelist_push(&ibmm.rcv_migrate[old_pid], node);
      }

#if TREE
      int ig = 0, jg = 0, kg = 0;
      NOT_UNUSED(ig);
      NOT_UNUSED(jg);
      NOT_UNUSED(kg);
      coord d = {0};
      foreach_dimension()
        d.x = ibval(npos.x);
      coord_periodic_boundary(d);
      int effective_depth = node->depth > depth() ? depth() : node->depth;
      Point point = locate_level(d.x, d.y, d.z, effective_depth);
      POINT_VARIABLES();

      if (node->pid == pid()) {
        ibnodelist_push(&ibmm.local, node);
        if (point.level >= 0 && allocated(0)) {
          foreach_neighbor(PESKIN_SUPPORT_RADIUS) {
            if (allocated(0) && !is_local(cell) &&
                cell.pid >= 0 && cell.pid < np)
              ibexchangelist_push_unique(&ibmm.snd_boundary[cell.pid], node);
          }
        }
      } else {
        bool has_local_support = false;
        if (point.level >= 0 && allocated(0)) {
          foreach_neighbor(PESKIN_SUPPORT_RADIUS) {
            if (!has_local_support && allocated(0) && is_local(cell))
              has_local_support = true;
          }
        }
        if (has_local_support && node->pid >= 0 && node->pid < np)
          ibexchangelist_push_unique(&ibmm.rcv_boundary[node->pid], node);
      }
#else // !TREE
      coord d = {0};
      foreach_dimension()
        d.x = ibval(npos.x);
      coord_periodic_boundary(d);
      Point point = locate_nonlocal(d.x, d.y, d.z);

      int ig = 0, jg = 0, kg = 0;
      NOT_UNUSED(ig);
      NOT_UNUSED(jg);
      NOT_UNUSED(kg);

      if (point.level >= 0) {
        if (node->pid == pid()) {
          ibnodelist_push(&ibmm.local, node);
          foreach_neighbor(PESKIN_SUPPORT_RADIUS) {
            if (is_boundary(point)) {
              int point_pid = _ibmeshmanager_get_pid(point);
              if (point_pid >= 0 && point_pid < np)
                ibexchangelist_push_unique(&ibmm.snd_boundary[point_pid], node);
            }
          }
        } else if (node->pid >= 0 && node->pid < np) {
          ibexchangelist_push_unique(&ibmm.rcv_boundary[node->pid], node);
        }
      }
#endif
    }

    free(node_pids);
  }

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
