#pragma once

#include "elff/config/_config.hpp"
#include "elff/io/C/xmlCase.h"

#include <stddef.h>

#ifdef ELFF_USE_MPI
#include <mpi.h>
#endif

#ifdef __cplusplus
extern "C"
{
#endif

/**
 * @name model_runtime_t
 * @brief C handle and related functions for the model runtime.
 * @{
 */
typedef void* model_runtime_t;

typedef struct {
  size_t model_index;
  long line;
  long column;
  char message[512];
} model_runtime_error;

/** @brief Creates a new model runtime instance.
 *
 * @return Opaque runtime handle
 */
model_runtime_t model_runtime_new();

/** Build and seal a serial runtime from the ordered models in a case file.
 * Renders expressions in memory before validation and model construction; writes
 * a best-effort .rendered.xml diagnostic beside the source (rank zero for MPI).
 * Models distribution is rank0 by default, or round-robin in XML.
 * No MPI communicator is configured; MPI callers should use the MPI factory.
 */
int model_runtime_create_from_case_file(const char* path,
                                     model_runtime_t* out_runtime,
                                     model_runtime_error* error);

/** Build from an already rendered in-memory Case document. The document is
 * borrowed only during construction; no file is opened or written. */
int model_runtime_create_from_case_document(elff_xml_document_t document,
                                         model_runtime_t* out_runtime,
                                         model_runtime_error* error);

#ifdef ELFF_USE_MPI
/** Collectively build a runtime on an initialized, non-null communicator.
 * The runtime duplicates comm before any model construction. All members must
 * call this function; a local case/model error fails construction on all ranks.
 * On failure, out_runtime is null. Destroy the runtime before MPI_Finalize.
 */
int model_runtime_create_from_case_file_mpi(const char* path, MPI_Comm comm,
                                         model_runtime_t* out_runtime,
                                         model_runtime_error* error);

/** Collective document factory; duplicates comm before constructing models.
 * Each rank supplies its rendered document, which is borrowed for this call. */
int model_runtime_create_from_case_document_mpi(elff_xml_document_t document,
                                             MPI_Comm comm,
                                             model_runtime_t* out_runtime,
                                             model_runtime_error* error);

/** @brief Sets the MPI communicator used by the runtime.
 *
 * @param runtime Runtime handle
 * This duplicates a non-null communicator collectively. Call before adopting
 * models when construction depends on rank context. MPI must be initialized.
 *
 * @param comm MPI communicator
 * @return Zero on success, nonzero on failure
 */
int model_runtime_set_communicator(model_runtime_t runtime, MPI_Comm comm);
#endif

/** @brief Saves the runtime state to disk.
 * Collective for mixed model owners; rank zero writes the ordered state file.
 * All communicator members must participate.
 *
 * @param runtime Runtime handle
 * @param fname Output checkpoint filename
 * @return Zero on success, nonzero on failure
 */
int model_runtime_checkpoint(model_runtime_t runtime, const char *fname);

/** @brief Restores the runtime state from disk.
 * With mixed owners, the file must be readable by every runtime rank.
 *
 * @param runtime Runtime handle
 * @param fname Input checkpoint filename
 * @return Zero on success, nonzero on failure
 */
int model_runtime_restore(model_runtime_t runtime, const char *fname);

/** @brief Writes all runtime-exportable model PolyData to a static VTKHDF file.
 *
 * Models that do not override the Model PolyData hook are skipped.
 *
 * @param runtime Runtime handle
 * @param fname Output VTKHDF filename
 * @param overwrite Nonzero to overwrite an existing file
 * @return Zero on success, nonzero on failure
 */
int model_runtime_write_polydata(model_runtime_t runtime,
                              const char *fname,
                              int overwrite);

/** @brief Destroys a runtime instance.
 *
 * @param runtime Runtime handle to destroy
 */
void model_runtime_delete(model_runtime_t runtime);
/** @} */

#ifdef __cplusplus
}
#endif
