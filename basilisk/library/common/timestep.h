#include "io/output-dump.h"

// ============================================================================
// Type definitions
// ============================================================================

struct CheckpointSidecar {
  int iter;
  int inext;
  double time;
  double dt;
  double tnext;
};

// ============================================================================
// Globals
// ============================================================================

static struct CheckpointSidecar cpsc = {-1, -1, 0., 0., 0.};

// ============================================================================
// Function declarations
// ============================================================================

static int write_checkpoint_sidecar (const char* file,
                                     double time,
                                     double dt_,
                                     int iter_,
                                     int inext_,
                                     double tnext_);
static int read_checkpoint_sidecar (const char* file,
                                    struct CheckpointSidecar* s);

static int timestep_checkpoint_dump (const char* path, void* ctx);
static int timestep_checkpoint_restore (const char* path, void* ctx);
double timestep (const face vector u, double dtmax);

// ============================================================================
// Events
// ============================================================================

event defaults (i = 0) {
  checkpointer_register (
    (Checkpointer) {.filename = ".dt",
                    .dump_phase = CKPT_PHASE_POST_DUMP,
                    .dump = timestep_checkpoint_dump,
                    .restore_phase = CKPT_PHASE_POST_RESTORE,
                    .restore = timestep_checkpoint_restore,
                    .ctx = NULL});
}

// ============================================================================
// Function definitions
// ============================================================================

/**
 * @brief
 */
static int write_checkpoint_sidecar (const char* file,
                                     double time,
                                     double dt_,
                                     int iter_,
                                     int inext_,
                                     double tnext_) {
  FILE* fp = fopen (file, "wb");
  if (!fp)
    return -1;

  struct CheckpointSidecar s = {
    .iter = iter_, .inext = inext_, .time = time, .dt = dt_, .tnext = tnext_};

  int ok = fwrite (&s, sizeof (s), 1, fp) == 1 ? 0 : -1;
  fclose (fp);
  return ok;
}

static int read_checkpoint_sidecar (const char* file,
                                    struct CheckpointSidecar* s) {
  FILE* fp = fopen (file, "rb");
  if (!fp)
    return -1;

  int ok = fread (s, sizeof (*s), 1, fp) == 1 ? 0 : -1;
  fclose (fp);

  return ok;
}

/**
 * @brief
 */
static int timestep_checkpoint_dump (const char* path, void* ctx) {
  (void) ctx;
  return write_checkpoint_sidecar (path, t, dt, 0, inext, tnext);
}

/**
 * @brief
 */
static int timestep_checkpoint_restore (const char* path, void* ctx) {
  (void) ctx;
  return read_checkpoint_sidecar (path, &cpsc);
}

/**
 * @brief
 */
double timestep (const face vector u, double dtmax) {
  static double previous = 0.;

  if (cpsc.iter >= 0) {
    previous = cpsc.dt;
    dt = cpsc.dt;
    cpsc.iter = -1;
    return dt;
  }

#if TIMESTEP_USE_DT_MAX
  previous = dtmax;
  return dtmax;
#endif

  if (t == 0.)
    previous = 0.;
  dtmax /= CFL;

  foreach_face (reduction (min : dtmax)) if (u.x[] != 0.) {
    double dt = Delta / fabs (u.x[]);
    assert (fm.x[]);
    dt *= fm.x[];
    if (dt < dtmax)
      dtmax = dt;
  }

  dtmax *= CFL;
  if (dtmax > previous)
    dtmax = (previous + 0.1 * dtmax) / 1.1;
  previous = dtmax;
  return dtmax;
}
