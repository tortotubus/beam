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

/** Build and seal a serial runtime from the ordered models in a case file.
 * Renders expressions in memory before validation and model construction; writes
 * a best-effort .rendered.xml diagnostic beside the source (rank zero for MPI).
 * Models distribution is rank0 by default, or round-robin in XML.
 * No MPI communicator is configured; MPI callers should use the MPI factory.
 */
int ib_runtime_create_from_case_file(const char* path,
                                     ib_runtime_t* out_runtime,
                                     ib_runtime_error* error);

/** Build from an already rendered in-memory Case document. The document is
 * borrowed only during construction; no file is opened or written. */
int ib_runtime_create_from_case_document(elff_xml_document_t document,
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
/** Collectively build a runtime on an initialized, non-null communicator.
 * The runtime duplicates comm before any model construction. All members must
 * call this function; a local case/model error fails construction on all ranks.
 * On failure, out_runtime is null. Destroy the runtime before MPI_Finalize.
 */
int ib_runtime_create_from_case_file_mpi(const char* path, MPI_Comm comm,
                                         ib_runtime_t* out_runtime,
                                         ib_runtime_error* error);

/** Collective document factory; duplicates comm before constructing models.
 * Each rank supplies its rendered document, which is borrowed for this call. */
int ib_runtime_create_from_case_document_mpi(elff_xml_document_t document,
                                             MPI_Comm comm,
                                             ib_runtime_t* out_runtime,
                                             ib_runtime_error* error);

/** @brief Sets the MPI communicator used by the runtime.
 *
 * @param runtime Runtime handle
 * This duplicates a non-null communicator collectively. Call before adopting
 * models when construction depends on rank context. MPI must be initialized.
 *
 * @param comm MPI communicator
 * @return Zero on success, nonzero on failure
 */
int ib_runtime_set_communicator(ib_runtime_t runtime, MPI_Comm comm);
#endif

/** @brief Saves the runtime state to disk.
 * Collective for mixed model owners; rank zero writes the ordered state file.
 * All communicator members must participate.
 *
 * @param runtime Runtime handle
 * @param fname Output checkpoint filename
 * @return Zero on success, nonzero on failure
 */
int ib_runtime_checkpoint(ib_runtime_t runtime, const char *fname);

/** @brief Restores the runtime state from disk.
 * With mixed owners, the file must be readable by every runtime rank.
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
