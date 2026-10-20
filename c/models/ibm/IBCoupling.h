#pragma once
#include "elff/c/models/ModelRuntime.h"
#ifdef __cplusplus
extern "C" {
#endif
/** @brief Borrowed runtime-owned IB component; never delete independently. */
typedef void* ib_coupling_t;
typedef enum {
  IB_MODEL_VIEW_INVALID = 0,
  IB_MODEL_VIEW_FORCE_COUPLED,
  IB_MODEL_VIEW_VELOCITY_COUPLED
} ib_model_view_kind_t;

/** Stable metadata handle, borrowed from a sealed runtime on every rank. */
typedef const void* ib_model_descriptor_t;
typedef struct {
  int sparse;
  int owner_rank; /* Rank in the runtime communicator. */
  int level;
  int node_count;
} ib_model_info_t;
int ib_model_descriptor_get_info(ib_model_descriptor_t descriptor,
                                 ib_model_info_t* info);

typedef struct {
  ib_model_view_kind_t kind;
  void* handle; /* Borrowed local model; NULL on remote sparse ranks. */
  ib_model_descriptor_t descriptor; /* Present on every rank. */
} ib_model_view_t;

typedef struct {
  size_t count;
  const ib_model_view_t* data;
} ib_model_view_array_t;


/** @brief Attach or retrieve the IB component; a new attachment requires an unsealed runtime. */
ib_coupling_t ib_coupling_attach(model_runtime_t runtime);
/** @brief Retrieve an existing component without attaching one. */
ib_coupling_t ib_coupling_get(model_runtime_t runtime);
/** @brief Seal the owning runtime and expose stable ordered IB-only views. */
int ib_coupling_get_model_views(ib_coupling_t coupling, ib_model_view_array_t* views);
#ifdef __cplusplus
}
#endif
