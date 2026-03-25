#pragma once

#if _MPI
#define ELFF_USE_MPI
#endif

#include "library/io/output-dump.h"
#include "elff/c/models/ibm/IBRuntime.h"

// ============================================================================
// Globals
// ============================================================================

static ib_runtime_t elff_runtime = NULL;

// ============================================================================
// Function declarations
// ============================================================================

static inline ib_runtime_t elff_runtime_get ();
static inline void elff_runtime_init ();
static inline int elff_runtime_register (ib_model_t model, int pid);
static inline int elff_dump (const char* fname);
static inline int elff_checkpoint_dump (const char* path, void* ctx);
static inline int elff_restore (const char* fname);
static inline int elff_checkpoint_restore (const char* path, void* ctx);
static inline void elff_runtime_free ();

// ============================================================================
// Function definitions
// ============================================================================

/**
 * @brief
 */
static inline ib_runtime_t elff_runtime_get () {
  if (!elff_runtime) {
    elff_runtime = ib_runtime_new ();
#if _MPI && defined(ELFF_USE_MPI)
    ib_runtime_set_communicator (elff_runtime, MPI_COMM_WORLD);
#endif
    static bool registered = false;
    if (!registered) {
      checkpointer_register (
        (Checkpointer) {.filename = ".elff",
                        .dump_phase = CKPT_PHASE_POST_DUMP,
                        .dump = elff_checkpoint_dump,
                        .restore_phase = CKPT_PHASE_POST_RESTORE,
                        .restore = elff_checkpoint_restore,
                        .ctx = NULL});
      registered = true;
    }
  }

  assert (elff_runtime);
  return elff_runtime;
}

/**
 * @brief
 */
static inline void elff_runtime_init () {
  (void) elff_runtime_get ();
}

/**
 * @brief
 */
static inline int elff_runtime_register (ib_model_t model, int pid = 0) {
  if (!model)
    return -1;

  ib_runtime_t runtime = elff_runtime_get ();
  int id = ib_runtime_register (runtime, model);
  if (id < 0)
    return id;

  if (ib_runtime_set_pid (runtime, model, pid) != 0)
    return -1;

  return id;
}

/**
 * @brief
 */
static inline int elff_dump (const char* fname) {
  if (!fname)
    return -1;
  if (!elff_runtime)
    return -1;

  return ib_runtime_checkpoint (elff_runtime_get (), fname);
}

/**
 * @brief
 */
static inline int elff_checkpoint_dump (const char* path, void* ctx) {
  (void) ctx;
  return elff_dump (path);
}

/**
 * @brief
 */
static inline int elff_restore (const char* fname) {
  if (!fname)
    return -1;
  if (!elff_runtime)
    return -1;

  return ib_runtime_restore (elff_runtime_get (), fname);
}

/**
 * @brief
 */
static inline int elff_checkpoint_restore (const char* path, void* ctx) {
  (void) ctx;
  return elff_restore (path);
}

/**
 * @brief
 */
static inline void elff_runtime_free () {
  if (!elff_runtime)
    return;
  ib_runtime_delete (elff_runtime);
  elff_runtime = NULL;
}
