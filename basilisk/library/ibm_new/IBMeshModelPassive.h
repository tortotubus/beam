#pragma once

#include "IBMeshModel.h"

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
    model.force_ops->sync_node = ibmeshmodel_passive_velocity_sync_node;
    model.force_ops->input_node = ibmeshmodel_passive_velocity_input_node;
    model.force_ops->input_fields = ibmeshmodel_passive_input_fields;
    model.force_ops->output_fields = ibmeshmodel_passive_output_fields;
    model.force_ops->advance = ibmeshmodel_passive_force_advance;
    model.force_ops->destroy = ibmeshmodel_passive_destroy;
    return model;
  }

  ibmeshmodel_passive_destroy(passive);
  return ibmeshmodel_init();
}
