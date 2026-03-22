/**
 * README: To use this framework, simply use @c restore_handler("my_folder") inside 
 * and your @c init event as well as @c checkpoint_hanlder("my_folder") in another
 * event. By default, the built-in basilisk @c dump() function is always called.
 * 
 * If you have written something which also requires dumping and restoring of state
 * you may also have it managed by these handlers: Simply write your dump and restore
 * functions using the signatures of @c checkpoint_function and register it using
 * @c checkpointer_register to have it managed for you.
 */

#include "library/io/output-common.h"
#include "output.h"

// ============================================================================
// Type definitions
// ============================================================================

typedef int (*checkpoint_function) (const char* path, void* ctx);

enum {
  CKPT_PHASE_PRE_DUMP = 0,
  CKPT_PHASE_POST_DUMP = 1,
  CKPT_PHASE_PRE_RESTORE = 2,
  CKPT_PHASE_POST_RESTORE = 3,
};

/**
 * @brief Generic checkpoint function + context
 */
typedef struct {
  const char* filename;
  int dump_phase;
  checkpoint_function dump;
  int restore_phase;
  checkpoint_function restore;
  void* ctx;
} Checkpointer;

/**
 * @brief Configuration struct for under what conditions a checkpoint is
 * triggered
 */
typedef struct {
  bool checkpoint_on_wall_time;
  double checkpoint_on_wall_time_seconds;
  bool checkpoint_on_sim_time;
  double checkpoint_on_sim_time_seconds;
  bool checkpoint_on_sim_iter;
  int checkpoint_on_sim_iter_iterations;
  bool checkpoint_on_signal;
  int checkpoint_on_signals[6];
} CheckpointConfiguraiton;

// ============================================================================
// Globals
// ============================================================================

static struct {
  int n;
  Checkpointer items[100];
} checkpointers = {0};

static CheckpointConfiguraiton checkpoint_configuration = {
  .checkpoint_on_signal = true,
  .checkpoint_on_signals = {SIGINT, SIGTERM},
  .checkpoint_on_wall_time = true,
  .checkpoint_on_wall_time_seconds = 300,
  .checkpoint_on_sim_time = false,
  .checkpoint_on_sim_iter = false};

// ============================================================================
// Macro definitions
// ============================================================================
 
// ============================================================================
// Function declarations
// ============================================================================

static inline int checkpointer_register (Checkpointer c);

int checkpoint_handler (double sim_time_current,
                        int sim_iter_current,
                        const char* basepath,
                        scalar* slist);

int restore_handler (const char* basepath, scalar* slist);

// ============================================================================
// Function definitions
// ============================================================================

/**
 * @brief
 * @relates Checkpointer
 */
static inline int checkpointer_register (Checkpointer c) {
  assert (checkpointers.n < 100);
  checkpointers.items[checkpointers.n] = c;
  return checkpointers.n++;
}

#include "library/io/input-signals.h"

/**
 * @brief Handler function can be called and will check for the conditions in
 * @c checkpoint_configuration and call registered dump functions when met. If
 * any registered signals (e.g. SIGTERM) are called it will abort after dumping
 * 
 * @relates CheckpointConfiguraiton
 */
int checkpoint_handler (double sim_time_current = t,
                        int sim_iter_current = i,
                        const char* basepath = NULL,
                        scalar* slist = all) {
  char fname[1024];
  char pname[1024];
  char tname[1024];

  if (!basepath) {
    snprintf (pname, 1024, "%s", get_unique_basename ());
    snprintf (fname, 1024, "%s/dump", get_unique_basename ());
    snprintf (tname, 1024, "%s/dump.dt", get_unique_basename ());
  } else {
    snprintf (pname, 1024, "%s", basepath);
    snprintf (fname, 1024, "%s/dump", basepath);
    snprintf (tname, 1024, "%s/dump.dt", basepath);
  }

  assert (!create_path (pname));

  static double wall_time_last = -1.;
  static double sim_time_last = 0.;
  static int sim_iter_last = 0;
  bool should_dump = false;
  int shutdown_sig = 0;

#if _MPI
  double wall_time_current = MPI_Wtime ();
  if (wall_time_last < 0.)
    wall_time_last = wall_time_current;

  shutdown_sig =
    checkpoint_configuration.checkpoint_on_signal ? shutdown_signal () : 0;
  int local_should_dump = shutdown_sig != 0;

  if (!local_should_dump && checkpoint_configuration.checkpoint_on_wall_time &&
      wall_time_current - wall_time_last >
        checkpoint_configuration.checkpoint_on_wall_time_seconds) {
    wall_time_last = wall_time_current;
    local_should_dump = true;
  } else if (!local_should_dump &&
             checkpoint_configuration.checkpoint_on_sim_time &&
             (double) (sim_time_current - sim_time_last) >
               checkpoint_configuration.checkpoint_on_sim_time_seconds) {
    sim_time_last = sim_time_current;
    local_should_dump = true;
  } else if (!local_should_dump &&
             checkpoint_configuration.checkpoint_on_sim_iter &&
             (double) (sim_iter_current - sim_iter_last) >
               checkpoint_configuration.checkpoint_on_sim_iter_iterations) {
    sim_iter_last = sim_iter_current;
    local_should_dump = true;
  }

  int global_should_dump = 0;
  MPI_Allreduce (
    &local_should_dump, &global_should_dump, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  should_dump = global_should_dump != 0;

  MPI_Allreduce (&shutdown_sig, &shutdown_sig, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#else
  double wall_time_current = (double) time (NULL);
  if (wall_time_last < 0.)
    wall_time_last = wall_time_current;

  // ======================
  // Checkpoint on signal
  // ======================
  shutdown_sig =
    checkpoint_configuration.checkpoint_on_signal ? shutdown_signal () : 0;
  if (shutdown_sig) {
    should_dump = true;
  }

  // ======================
  // Checkpoint on walltime
  // ======================
  else if (checkpoint_configuration.checkpoint_on_wall_time &&
           (long) (wall_time_current - wall_time_last) >
             checkpoint_configuration.checkpoint_on_wall_time_seconds) {
    wall_time_last = wall_time_current;
    should_dump = true;
  }

  // =============================
  // Checkpoint on simulation time
  // =============================
  else if (checkpoint_configuration.checkpoint_on_sim_time &&
           (double) (sim_time_current - sim_time_last) >
             checkpoint_configuration.checkpoint_on_sim_time_seconds) {
    sim_time_last = sim_time_current;
    should_dump = true;
  }

  // =============================
  // Checkpoint on simulation time
  // =============================
  else if (checkpoint_configuration.checkpoint_on_sim_iter &&
           (double) (sim_iter_current - sim_iter_last) >
             checkpoint_configuration.checkpoint_on_sim_iter_iterations) {
    sim_iter_last = sim_iter_current;
    should_dump = true;
  }
#endif

  if (should_dump) {
    printf ("\nWriting dumpfile %s \n", fname);
    for (int ci = 0; ci < checkpointers.n; ci++) {
      Checkpointer* c = &checkpointers.items[ci];
      if (!c->dump || c->dump_phase != CKPT_PHASE_PRE_DUMP)
        continue;

      char path[1024];
      snprintf (path,
                1024,
                "%s%s",
                fname,
                c->filename ? c->filename : "");

      if (c->dump (path, c->ctx) != 0) {
        fprintf (
          ferr, "checkpoint_handler(): pre-dump checkpointer failed: %s\n", path);
        exit (1);
      }
    }

    dump (fname, slist, NULL, false, true);

    for (int ci = 0; ci < checkpointers.n; ci++) {
      Checkpointer* c = &checkpointers.items[ci];
      if (!c->dump || c->dump_phase != CKPT_PHASE_POST_DUMP)
        continue;

      char path[1024];
      snprintf (path,
                1024,
                "%s%s",
                fname,
                c->filename ? c->filename : "");

      if (c->dump (path, c->ctx) != 0) {
        fprintf (ferr,
                 "checkpoint_handler(): post-dump checkpointer failed: %s\n",
                 path);
        exit (1);
      }
    }
  }

  if (shutdown_sig) {
    if (pid () == 0)
      fprintf (stderr, "warning: Shutdown requested (signal %d)\n", shutdown_sig);
    return 1;
  }

  return 0;
}

/**
 * @brief Restoration handler that handles restoration analagously to basilisk's
 * built-in @c restore() routine.
 * 
 * @relates CheckpointConfiguraiton
 */
int restore_handler (const char* basepath = NULL, scalar* slist = all) {
  install_shutdown_handlers ();

  char fname[1024];
  char pname[1024];
  char tname[1024];

  if (!basepath) {
    snprintf (pname, 1024, "%s", get_unique_basename ());
    snprintf (fname, 1024, "%s/dump", get_unique_basename ());
    snprintf (tname, 1024, "%s/dump.dt", get_unique_basename ());
  } else {
    snprintf (pname, 1024, "%s", basepath);
    snprintf (fname, 1024, "%s/dump", basepath);
    snprintf (tname, 1024, "%s/dump.dt", basepath);
  }

  if (file_exists (fname)) {
    printf ("Resuming simulation from checkpoint file %s\n", fname);

    for (int ci = 0; ci < checkpointers.n; ci++) {
      Checkpointer* c = &checkpointers.items[ci];
      if (!c->restore || c->restore_phase != CKPT_PHASE_PRE_RESTORE)
        continue;

      char path[1024];
      snprintf (path,
                1024,
                "%s%s",
                fname,
                c->filename ? c->filename : "");

      if (c->restore (path, c->ctx) != 0) {
        fprintf (ferr,
                 "restore_handler(): pre-restore checkpointer failed: %s\n",
                 path);
        exit (1);
      }
    }

    restore (fname, slist, NULL);

    for (int ci = 0; ci < checkpointers.n; ci++) {
      Checkpointer* c = &checkpointers.items[ci];
      if (!c->restore || c->restore_phase != CKPT_PHASE_POST_RESTORE)
        continue;

      char path[1024];
      snprintf (path,
                1024,
                "%s%s",
                fname,
                c->filename ? c->filename : "");

      if (c->restore (path, c->ctx) != 0) {
        fprintf (ferr,
                 "restore_handler(): post-restore checkpointer failed: %s\n",
                 path);
        exit (1);
      }
    }
    return 1;
  } else {
    return 0;
  }
}
