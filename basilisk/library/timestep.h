#include "library/io/output-dump.h"

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
// CFL growth history before dtnext() clips dt to scheduled output times.
static double timestep_previous = 0.;
static double checkpoint_previous = 0.;

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
  if (ok == 0 && fwrite (&timestep_previous, sizeof (double), 1, fp) != 1)
    ok = -1;
  if (fclose (fp) != 0)
    ok = -1;
  return ok;
}

static int read_checkpoint_sidecar (const char* file,
                                    struct CheckpointSidecar* s) {
  FILE* fp = fopen (file, "rb");
  if (!fp)
    return -1;

  int ok = fread (s, sizeof (*s), 1, fp) == 1 ? 0 : -1;
  if (ok == 0) {
    // Legacy sidecars end after the clock record. New files append the CFL
    // history, which is generally larger than the output-aligned saved dt.
    checkpoint_previous = s->dt;
    double previous;
    size_t n = fread (&previous, 1, sizeof (previous), fp);
    if (n == sizeof (previous))
      checkpoint_previous = previous;
    else if (n != 0 || ferror (fp))
      ok = -1;
  }
  fclose (fp);

  return ok;
}

/**
 * @brief
 */
static int timestep_checkpoint_dump (const char* path, void* ctx) {
  (void) ctx;
  // All ranks share the same clock; only rank zero writes the shared sidecar.
  int rc = pid () == 0 ?
    write_checkpoint_sidecar (path, t, dt, iter, inext, tnext) : 0;
#if _MPI
  MPI_Bcast (&rc, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  return rc;
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

  if (cpsc.iter >= 0) {
    timestep_previous = checkpoint_previous;
    dt = cpsc.dt;
    cpsc.iter = -1;
    return dt;
  }
  if (t == 0.)
    timestep_previous = 0.;
  dtmax /= CFL;

  foreach_face (reduction (min : dtmax)) if (u.x[] != 0.) {
    double dt = Delta / fabs (u.x[]);
    assert (fm.x[]);
    dt *= fm.x[];
    if (dt < dtmax)
      dtmax = dt;
  }

  dtmax *= CFL;
  if (dtmax > timestep_previous)
    dtmax = (timestep_previous + 0.1 * dtmax) / 1.1;
  timestep_previous = dtmax;
  return dtmax;
}
