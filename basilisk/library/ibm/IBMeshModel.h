#pragma once

// #include "library/ibm/IBMesh.h"

/**
 * @enum IBMeshModelType;
 */
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
 * @note \c mesh arguments are really supposed to be IBMesh type. However, qcc
 * does not play nice with forward declaration of struct. In the function it
 * points to, implememters should treat and cast mesh to IBMesh
 */
typedef struct {
  int (*node_count) (void* ctx);
  int (*sync) (void* ctx, void* mesh);
  int (*advance) (void* ctx, void* mesh, double dt);
  void (*destroy) (void* ctx);
} IBVelocityCoupledModelOps;

/**
 * @struct IBForceCoupledModelOps
 */
typedef struct {
  int (*node_count) (void* ctx);
  void (*sync) (void* ctx, void* mesh);
  void (*advance) (void* ctx, void* mesh, double dt);
  void (*destroy) (void* ctx);
} IBForceCoupledModelOps;

/**
 * @struct IBMeshModel
 */
typedef struct {
  int type;
  IBForceCoupledModelOps* force_ops;
  IBVelocityCoupledModelOps* velocity_ops;
  void* ctx;
} IBMeshModel;

/**
 * @relates IBMeshModel
 */
IBMeshModel ibmeshmodel_init () {
  IBMeshModel model = {0};
  return model;
}

/**
 * @relates IBMeshModel
 */
IBMeshModel ibmeshmodel_force_coupled_init () {
  IBMeshModel model = {0};
  model.type = IB_MODEL_FORCE_COUPLED;
  model.force_ops = (IBForceCoupledModelOps*) calloc (1, sizeof (IBForceCoupledModelOps));
  return model;
}

/**
 * @relates IBMeshModel
 */
void ibmeshmodel_destroy (IBMeshModel* model) {
  if (!model)
    return;

  switch (model->type) {
  case IB_MODEL_FORCE_COUPLED: {
    if (model->force_ops) {
      if (model->force_ops->destroy)
        model->force_ops->destroy (model->ctx);
      free (model->force_ops);
    }
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
