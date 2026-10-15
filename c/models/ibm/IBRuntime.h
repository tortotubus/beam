#pragma once

#ifdef ELFF_USE_MPI
#include <mpi.h>
#endif

#ifdef __cplusplus
extern "C"
{
#endif

/**
 * @name ib_model_t
 * @brief Opaque handle to a model registered with the IB runtime.
 * @{
 */
typedef void* ib_model_t;
/** @} */

/**
 * @name ib_runtime_t
 * @brief C handle and related functions for the immersed-boundary runtime.
 * @{
 */
typedef void* ib_runtime_t;

/** @brief Creates a new immersed-boundary runtime instance.
 *
 * @return Opaque runtime handle
 */
ib_runtime_t ib_runtime_new();

/** @brief Registers a model with the runtime.
 *
 * @param runtime Runtime handle
 * @param model Model handle to register
 * @return Zero on success, nonzero on failure
 */
int ib_runtime_register(ib_runtime_t runtime, ib_model_t model);

/** @brief Deregisters a model from the runtime.
 *
 * @param runtime Runtime handle
 * @param model Model handle to deregister
 * @return Zero on success, nonzero on failure
 */
int ib_runtime_deregister(ib_runtime_t runtime, ib_model_t model);

/** @brief Assigns a process id to a registered model.
 *
 * @param runtime Runtime handle
 * @param model Model handle
 * @param pid Process id to associate with the model
 * @return Zero on success, nonzero on failure
 */
int ib_runtime_set_pid(ib_runtime_t runtime, ib_model_t model, int pid);
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
