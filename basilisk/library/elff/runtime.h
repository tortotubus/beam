#pragma once

#if _MPI
#ifndef ELFF_USE_MPI
#define ELFF_USE_MPI 1
#endif
#endif

#include "library/io/output-dump.h"
#include "library/io/vtk/vtkTimeSeries.h"
#include "elff/c/models/ibm/IBRuntime.h"
#include "elff/c/general/output.h"

// ============================================================================
// Globals
// ============================================================================

static ib_runtime_t elff_runtime = NULL;
static ib_model_view_array_t elff_model_views = {0};

// ============================================================================
// Function declarations
// ============================================================================

static inline ib_runtime_t elff_runtime_get ();
static inline void elff_runtime_init ();
static inline int elff_init (const char* case_path, ib_runtime_error* error);
static inline ib_model_view_array_t elff_models ();
static inline void elff_deinit ();
static inline int elff_dump (const char* fname);
static inline int elff_write_polydata (const char* fname, int overwrite);
static inline int output_hdf_elff_pd_series (const char* basename,
                                             int iter,
                                             double time,
                                             bool overwrite);
static inline int output_hdf_elff_pd (const char* basename,
                                      int iter,
                                      double time,
                                      bool overwrite);
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
  // elff_set_out_prefix("[info] ");
  return elff_runtime;
}

/** Create and seal the global runtime from the ordered <ELFF><Models> list. */
static inline int elff_init (const char* case_path, ib_runtime_error* error) {
  if (!case_path || elff_runtime)
    return -1;
  if (ib_runtime_create_from_case_file (case_path, &elff_runtime, error) != 0)
    return -1;
#if _MPI && defined(ELFF_USE_MPI)
  if (ib_runtime_set_communicator (elff_runtime, MPI_COMM_WORLD) != 0) {
    ib_runtime_delete (elff_runtime);
    elff_runtime = NULL;
    return -1;
  }
#endif
  if (ib_runtime_get_model_views (elff_runtime, &elff_model_views) != 0) {
    ib_runtime_delete (elff_runtime);
    elff_runtime = NULL;
    elff_model_views = (ib_model_view_array_t) {0};
    return -1;
  }

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
  return 0;
}

/** Return the stable borrowed view array for the global sealed runtime. */
static inline ib_model_view_array_t elff_models () {
  ib_model_view_array_t views = {0};
  if (elff_runtime &&
      ib_runtime_get_model_views (elff_runtime, &views) == 0)
    return views;
  return views;
}

/** Destroy all runtime-owned models after Basilisk adapters have been freed. */
static inline void elff_deinit () {
  if (elff_runtime)
    ib_runtime_delete (elff_runtime);
  elff_runtime = NULL;
  elff_model_views = (ib_model_view_array_t) {0};
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
static inline int elff_write_polydata (const char* fname, int overwrite = 1) {
  if (!fname)
    return -1;
  if (!elff_runtime)
    return -1;

  return ib_runtime_write_polydata (elff_runtime_get (), fname, overwrite);
}

/**
 * @brief Write runtime-exported ELFF PolyData and update the series file.
 */
static inline int output_hdf_elff_pd_series (const char* basename,
                                             int iter = i,
                                             double time = t,
                                             bool overwrite = true) {
  if (!basename)
    return -1;

  char dirname[4096];
  snprintf (dirname, sizeof (dirname), "%s/vtkhdf-elff-pd-series", basename);

  char series_entry[4096];
  snprintf (series_entry,
            sizeof (series_entry),
            "vtkhdf-elff-pd-series/polydata_%d.vtkhdf",
            iter);

  char series_filename[4096];
  snprintf (series_filename,
            sizeof (series_filename),
            "%s/elff-pd.vtkhdf.series",
            basename);

  if (pid () == 0) {
    if (create_path (basename) != 0)
      return -1;
    if (create_path (dirname) != 0)
      return -1;
    output_vtk_series (series_entry, series_filename, iter, time);
  }

  char fname[8192];
  snprintf (fname, sizeof (fname), "%s/polydata_%d.vtkhdf", dirname, iter);
  return elff_write_polydata (fname, overwrite ? 1 : 0);
}

/**
 * @brief Write runtime-exported ELFF PolyData using the standard series layout.
 */
static inline int output_hdf_elff_pd (const char* basename = NULL,
                                      int iter = i,
                                      double time = t,
                                      bool overwrite = true) {
  char basename_buff[64];

  if (!basename) {
    snprintf (
      basename_buff, sizeof (basename_buff), "%s", get_unique_basename ());
  } else {
    snprintf (basename_buff, sizeof (basename_buff), "%s", basename);
  }

  return output_hdf_elff_pd_series (basename_buff, iter, time, overwrite);
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
  elff_deinit ();
}
