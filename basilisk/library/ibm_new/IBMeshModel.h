

// ============================================================================
// Type Definitions
// ============================================================================

enum {
  IB_MODEL_INVALID = -1,
  IB_MODEL_NONE = 0,
  IB_MODEL_VELOCITY_COUPLED = 1,
  IB_MODEL_FORCE_COUPLED = 2,
  IB_MODEL_OTHER = 3,
};

/**
 * @struct IBVelocityCoupledModelOps
 *
 * @brief Callback table for immersed-boundary models whose positions are
 * advanced from gathered Eulerian velocity.
 *
 * @details The mesh manager owns the gather/spread/MPI lifecycle. A velocity
 * coupled model owns only model-specific state and implements the callbacks
 * needed to copy that state into IB node fields and advance it.
 *
 * @param node_count Returns the number of IB nodes required by the model. This
 * is queried during @ref ibmeshmanager_set_model() before nodes are allocated.
 * All ranks must return the same value for a logical mesh.
 * @param sync Copies model-owned state into IB node fields without advancing
 * time. This is called during model registration and by explicit manager sync
 * hooks. In MPI, only the mesh owner is expected to run the authoritative sync.
 * @param sync_node Optional keyed sparse sync. Copies model-owned state for one
 * logical node id into a supplied temporary/local IBNode. Sparse-managed meshes
 * use this to avoid requiring a dense owner-side IBNode array.
 * @param input_node Optional keyed sparse input. Copies reduced node input
 * fields for one logical node id back into model-owned state on the mesh owner.
 * @param input_fields Optional heap-allocated IBscalar list reduced from node
 * owners to the mesh owner before owner-only callbacks. Return NULL to use
 * manager defaults. The manager frees the returned list.
 * @param output_fields Optional heap-allocated IBscalar list broadcast from the
 * mesh owner after sync/advance/midpoint. Return NULL to use manager defaults.
 * The manager frees the returned list.
 * @param midpoint Optional velocity-coupled midpoint update. Called only on
 * the mesh owner in MPI.
 * @param advance Full position/model update. Called only on the mesh owner in
 * MPI.
 * @param destroy Frees model-owned context. Non-owner ranks may have this
 * called immediately during distributed registration after node counts and
 * field lists have been sampled.
 *
 * @note @c mesh arguments are really supposed to be IBMesh type. However, qcc
 * does not play nice with forward declaration of struct. In the function it
 * points to, implementers should treat and cast mesh to IBMesh.
 */
typedef struct {
  int (*node_count)(void *ctx);
  int (*sync)(void *ctx, void *mesh);
  int (*sync_node)(void *ctx, void *mesh, int node_lid, void *node);
  int (*input_node)(void *ctx, void *mesh, int node_lid, void *node);
  IBscalar *(*input_fields)(void *ctx);
  IBscalar *(*output_fields)(void *ctx);
  int (*midpoint)(void *ctx, void *mesh, double dt);
  int (*advance)(void *ctx, void *mesh, double dt);
  void (*destroy)(void *ctx);
} IBVelocityCoupledModelOps;

/**
 * @struct IBPassiveModelCtx
 *
 * @brief Minimal context for a non-owner replica of an IB model.
 *
 * @details Passive models allocate the same IB node slots as the owner model
 * and preserve the same input/output field synchronization contract, but all
 * model callbacks are no-ops. They are used by the manager on MPI ranks that do
 * not own the model state. These ranks still need IB node fields for local grid
 * gather/spread/adaptivity, but should not retain expensive model-specific
 * state.
 *
 * @field node_count Number of IB nodes in the logical mesh.
 * @field input_fields Copied IBscalar list returned by the original model's
 * input_fields callback, or NULL to preserve default manager behavior.
 * @field output_fields Copied IBscalar list returned by the original model's
 * output_fields callback, or NULL to preserve default manager behavior.
 */
typedef struct {
  int node_count;
  IBscalar *input_fields;
  IBscalar *output_fields;
} IBPassiveModelCtx;

/**
 * @struct IBForceCoupledModelOps
 *
 * @brief Callback table for immersed-boundary models whose forces are coupled
 * to the Eulerian solver without velocity-driven position ownership.
 *
 * @details The ownership and field-list rules match
 * @ref IBVelocityCoupledModelOps. In MPI, sync/advance are owner-only and field
 * lists define the data reduced to, or broadcast from, that owner.
 */
typedef struct {
  int (*node_count)(void *ctx);
  void (*sync)(void *ctx, void *mesh);
  IBscalar *(*input_fields)(void *ctx);
  IBscalar *(*output_fields)(void *ctx);
  void (*advance)(void *ctx, void *mesh, double dt);
  void (*destroy)(void *ctx);
} IBForceCoupledModelOps;

/**
 * @struct IBMeshModel
 *
 * @brief Runtime model object attached to one IBMesh.
 *
 * @details Exactly one ops table should be active for the selected @c type.
 * @c ctx is owned by the model and released through the ops @c destroy callback
 * via @ref ibmeshmodel_destroy(). The mesh manager may replace a non-owner
 * model with a passive equivalent during registration.
 */
typedef struct {
  int type;
  IBForceCoupledModelOps *force_ops;
  IBVelocityCoupledModelOps *velocity_ops;
  void *ctx;
} IBMeshModel;

// ============================================================================
// Function Declarations
// ============================================================================

// ============================================================================
// Function Definitions
// ============================================================================

/**
 * @brief Return an empty model with type @ref IB_MODEL_NONE.
 *
 * @relates IBMeshModel
 */
IBMeshModel ibmeshmodel_init() {
  IBMeshModel model = {0};
  return model;
}

/**
 * @brief Allocate an empty force-coupled model object.
 *
 * @details The caller must fill the required callbacks before registering the
 * model with a mesh manager.
 *
 * @relates IBMeshModel
 */
IBMeshModel ibmeshmodel_force_coupled_init() {
  IBMeshModel model = {0};
  model.type = IB_MODEL_FORCE_COUPLED;
  model.force_ops =
      (IBForceCoupledModelOps *)calloc(1, sizeof(IBForceCoupledModelOps));
  return model;
}

/**
 * @brief Allocate an empty velocity-coupled model object.
 *
 * @details The caller must fill the required callbacks before registering the
 * model with a mesh manager.
 *
 * @relates IBMeshModel
 */
IBMeshModel ibmeshmodel_velocity_coupled_init() {
  IBMeshModel model = {0};
  model.type = IB_MODEL_VELOCITY_COUPLED;
  model.velocity_ops =
      (IBVelocityCoupledModelOps *)calloc(1, sizeof(IBVelocityCoupledModelOps));
  return model;
}

/**
 * @brief Destroy a model and reset it to @ref IB_MODEL_NONE.
 *
 * @details Calls the model-specific @c destroy callback on @c ctx if one is
 * provided, frees the ops table, and clears all pointers. Safe to call on empty
 * models.
 *
 * @relates IBMeshModel
 */
void ibmeshmodel_destroy(IBMeshModel *model) {
  if (!model)
    return;

  switch (model->type) {
  case IB_MODEL_VELOCITY_COUPLED: {
    if (model->velocity_ops) {
      if (model->velocity_ops->destroy)
        model->velocity_ops->destroy(model->ctx);
      free(model->velocity_ops);
    }
    break;
  }
  case IB_MODEL_FORCE_COUPLED: {
    if (model->force_ops) {
      if (model->force_ops->destroy)
        model->force_ops->destroy(model->ctx);
      free(model->force_ops);
    }
    break;
  }
  default: {
    break;
  }
  }

  model->force_ops = NULL;
  model->velocity_ops = NULL;
  model->ctx = NULL;
  model->type = IB_MODEL_NONE;
}

/** @brief Passive callback returning the replica node count. */
static int ibmeshmodel_passive_node_count(void *ctx) {
  IBPassiveModelCtx *passive = (IBPassiveModelCtx *)ctx;
  return passive ? passive->node_count : 0;
}

/** @brief Passive velocity sync callback; intentionally performs no work. */
static int ibmeshmodel_passive_velocity_sync(void *ctx, void *mesh) {
  (void)ctx;
  (void)mesh;
  return 0;
}

/** @brief Passive keyed sync callback; intentionally performs no work. */
static int ibmeshmodel_passive_velocity_sync_node(void *ctx, void *mesh,
                                                  int node_lid, void *node) {
  (void)ctx;
  (void)mesh;
  (void)node_lid;
  (void)node;
  return 0;
}

/** @brief Passive keyed input callback; intentionally performs no work. */
static int ibmeshmodel_passive_velocity_input_node(void *ctx, void *mesh,
                                                   int node_lid, void *node) {
  (void)ctx;
  (void)mesh;
  (void)node_lid;
  (void)node;
  return 0;
}

/** @brief Passive force sync callback; intentionally performs no work. */
static void ibmeshmodel_passive_force_sync(void *ctx, void *mesh) {
  (void)ctx;
  (void)mesh;
}

/**
 * @brief Return a copy of the preserved passive input field list.
 *
 * @details The caller owns and must free the returned list.
 */
static IBscalar *ibmeshmodel_passive_input_fields(void *ctx) {
  IBPassiveModelCtx *passive = (IBPassiveModelCtx *)ctx;
  return passive ? iblist_copy(passive->input_fields) : NULL;
}

/**
 * @brief Return a copy of the preserved passive output field list.
 *
 * @details The caller owns and must free the returned list.
 */
static IBscalar *ibmeshmodel_passive_output_fields(void *ctx) {
  IBPassiveModelCtx *passive = (IBPassiveModelCtx *)ctx;
  return passive ? iblist_copy(passive->output_fields) : NULL;
}

/** @brief Passive velocity midpoint callback; intentionally performs no work.
 */
static int ibmeshmodel_passive_velocity_midpoint(void *ctx, void *mesh,
                                                 double dt) {
  (void)ctx;
  (void)mesh;
  (void)dt;
  return 0;
}

/** @brief Passive velocity advance callback; intentionally performs no work. */
static int ibmeshmodel_passive_velocity_advance(void *ctx, void *mesh,
                                                double dt) {
  (void)ctx;
  (void)mesh;
  (void)dt;
  return 0;
}

/** @brief Passive force advance callback; intentionally performs no work. */
static void ibmeshmodel_passive_force_advance(void *ctx, void *mesh,
                                              double dt) {
  (void)ctx;
  (void)mesh;
  (void)dt;
}

/** @brief Free a passive model context and its preserved field lists. */
static void ibmeshmodel_passive_destroy(void *ctx) {
  IBPassiveModelCtx *passive = (IBPassiveModelCtx *)ctx;
  if (!passive)
    return;
  free(passive->input_fields);
  free(passive->output_fields);
  free(passive);
}

/**
 * @brief Build a passive replica of a velocity- or force-coupled model.
 *
 * @param model_type Original model type. Only velocity-coupled and
 * force-coupled models are supported.
 * @param node_count Number of IB nodes to allocate for the passive replica.
 * @param input_fields Optional field list sampled from the original model.
 * The list is copied; caller retains ownership of the input pointer.
 * @param output_fields Optional field list sampled from the original model.
 * The list is copied; caller retains ownership of the input pointer.
 *
 * @return A passive IBMeshModel with no-op sync/advance callbacks and preserved
 * field-list callbacks, or an empty model for unsupported model types.
 *
 * @details The manager uses this on non-owner MPI ranks after sampling the real
 * model's node count and field synchronization contract. This keeps proxy IB
 * node fields available locally while immediately releasing expensive
 * model-owned state on ranks that do not own it.
 */
static IBMeshModel ibmeshmodel_passive_from(int model_type, int node_count,
                                            IBscalar *input_fields,
                                            IBscalar *output_fields) {
  IBPassiveModelCtx *passive =
      (IBPassiveModelCtx *)calloc(1, sizeof(IBPassiveModelCtx));
  assert(passive);
  passive->node_count = node_count;
  passive->input_fields = iblist_copy(input_fields);
  passive->output_fields = iblist_copy(output_fields);

  if (model_type == IB_MODEL_VELOCITY_COUPLED) {
    IBMeshModel model = ibmeshmodel_velocity_coupled_init();
    model.ctx = passive;
    model.velocity_ops->node_count = ibmeshmodel_passive_node_count;
    model.velocity_ops->sync = ibmeshmodel_passive_velocity_sync;
    model.velocity_ops->sync_node = ibmeshmodel_passive_velocity_sync_node;
    model.velocity_ops->input_node = ibmeshmodel_passive_velocity_input_node;
    model.velocity_ops->input_fields = ibmeshmodel_passive_input_fields;
    model.velocity_ops->output_fields = ibmeshmodel_passive_output_fields;
    model.velocity_ops->midpoint = ibmeshmodel_passive_velocity_midpoint;
    model.velocity_ops->advance = ibmeshmodel_passive_velocity_advance;
    model.velocity_ops->destroy = ibmeshmodel_passive_destroy;
    return model;
  }

  if (model_type == IB_MODEL_FORCE_COUPLED) {
    IBMeshModel model = ibmeshmodel_force_coupled_init();
    model.ctx = passive;
    model.force_ops->node_count = ibmeshmodel_passive_node_count;
    model.force_ops->sync = ibmeshmodel_passive_force_sync;
    model.force_ops->input_fields = ibmeshmodel_passive_input_fields;
    model.force_ops->output_fields = ibmeshmodel_passive_output_fields;
    model.force_ops->advance = ibmeshmodel_passive_force_advance;
    model.force_ops->destroy = ibmeshmodel_passive_destroy;
    return model;
  }

  ibmeshmodel_passive_destroy(passive);
  return ibmeshmodel_init();
}
