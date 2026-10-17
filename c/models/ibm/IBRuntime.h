#pragma once

#include "elff/config/_config.hpp"

#include <stddef.h>

#ifdef ELFF_USE_MPI
#include <mpi.h>
#endif

#ifdef __cplusplus
extern "C"
{
#endif

/**
 * @name ib_runtime_t
 * @brief C handle and related functions for the immersed-boundary runtime.
 * @{
 */
typedef void* ib_runtime_t;

typedef enum {
  IB_MODEL_VIEW_INVALID = 0,
  IB_MODEL_VIEW_FORCE_COUPLED,
  IB_MODEL_VIEW_VELOCITY_COUPLED
} ib_model_view_kind_t;

typedef struct {
  ib_model_view_kind_t kind;
  void* handle; /* Borrowed; owned by the runtime. */
} ib_model_view_t;

typedef struct {
  size_t count;
  const ib_model_view_t* data;
} ib_model_view_array_t;

typedef struct {
  size_t model_index;
  long line;
  long column;
  char message[512];
} ib_runtime_error;

/** @brief Creates a new immersed-boundary runtime instance.
 *
 * @return Opaque runtime handle
 */
ib_runtime_t ib_runtime_new();

/** Build and seal a runtime from the ordered models in a case file. */
int ib_runtime_create_from_case_file(const char* path,
                                     ib_runtime_t* out_runtime,
                                     ib_runtime_error* error);

/**
 * Seal a runtime and return its ordered, borrowed coupling views.
 *
 * The array and handles remain stable until @ref ib_runtime_delete. Calling
 * this function prevents subsequent model adoption through C constructors.
 */
int ib_runtime_get_model_views(ib_runtime_t runtime,
                               ib_model_view_array_t* out_views);
#ifdef ELFF_USE_MPI
/** @brief Sets the MPI communicator used by the runtime.
 *
 * @param runtime Runtime handle
 * @param comm MPI communicator
 * @return Zero on success, nonzero on failure
 */
int ib_runtime_set_communicator(ib_runtime_t runtime, MPI_Comm comm);
#endif

/** @brief Saves the runtime state to disk.
 *
 * @param runtime Runtime handle
 * @param fname Output checkpoint filename
 * @return Zero on success, nonzero on failure
 */
int ib_runtime_checkpoint(ib_runtime_t runtime, const char *fname);

/** @brief Restores the runtime state from disk.
 *
 * @param runtime Runtime handle
 * @param fname Input checkpoint filename
 * @return Zero on success, nonzero on failure
 */
int ib_runtime_restore(ib_runtime_t runtime, const char *fname);

/** @brief Writes all runtime-exportable model PolyData to a static VTKHDF file.
 *
 * Models that do not override the IBModel PolyData hook are skipped.
 *
 * @param runtime Runtime handle
 * @param fname Output VTKHDF filename
 * @param overwrite Nonzero to overwrite an existing file
 * @return Zero on success, nonzero on failure
 */
int ib_runtime_write_polydata(ib_runtime_t runtime,
                              const char *fname,
                              int overwrite);

/** @brief Destroys a runtime instance.
 *
 * @param runtime Runtime handle to destroy
 */
void ib_runtime_delete(ib_runtime_t runtime);
/** @} */

#ifdef __cplusplus
}
#endif
