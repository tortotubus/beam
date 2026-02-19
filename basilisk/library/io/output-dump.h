#include "library/io/output-common.h"
#include "library/io/hdf/dump-hdf.h"

#ifndef OUTPUT_DUMP_H
#define OUTPUT_DUMP_H
#endif

#if _MPI
#define dump_function(fname, slist) dump ((fname), (slist), NULL, false, true)
#define restore_function(fname, slist) restore ((fname), (slist), NULL)
#else
#define dump_function(fname, slist) dump ((fname), (slist), NULL, false, true)
#define restore_function(fname, slist) restore ((fname), (slist), NULL)
// #define dump_function(fname, slist) dump_hdf ((fname), (slist))
// #define restore_function(fname, slist) restore_hdf ((fname), (slist))
#endif

static struct {
  bool restored;
  int restored_i;
  int restored_inext;
  double restored_t;
  double restored_tnext;
} checkpoint_state = {
  .restored = false,
  .restored_i = 0,
  .restored_inext = 0,
  .restored_t = 0,
  .restored_tnext = 0,
};

static struct {
  bool checkpoint_on_wall_time;
  double checkpoint_on_wall_time_seconds;
  bool checkpoint_on_sim_time;
  double checkpoint_on_sim_time_seconds;
  bool checkpoint_on_sim_iter;
  int checkpoint_on_sim_iter_iterations;
  bool checkpoint_on_signal;
  int checkpoint_on_signals[6];
} checkpoint_configuration = {.checkpoint_on_signal = true,
                              .checkpoint_on_signals = {SIGINT, SIGTERM},
                              .checkpoint_on_wall_time = true,
                              .checkpoint_on_wall_time_seconds = 300,
                              .checkpoint_on_sim_time = false,
                              .checkpoint_on_sim_iter = false};

#include "library/io/input-signals.h"

int checkpoint_handler (double sim_time_current = t,
                        int sim_iter_current = i,
                        const char* basename = NULL,
                        scalar* slist = all) {
  char fname[64];
  if (!basename) {
    snprintf (fname, 64, "%s.dump", get_executable_name ());
  } else {
    snprintf (fname, 64, "%s.dump", basename);
  }

#if _MPI
  static double t0 = 0, t1 = 0;
  if (t0 == 0) {
    t0 = MPI_Wtime ();
  }
  t1 = MPI_Wtime ();
#else

  static time_t wall_time_last = 0, wall_time_current = 0;
  static double sim_time_last = 0;
  static int sim_iter_last = 0;

  wall_time_current = time (NULL);
  if (wall_time_last == 0) {
    wall_time_last = time (NULL);
  }

  // ======================
  // Checkpoint on signal
  // ======================
  if (checkpoint_configuration.checkpoint_on_signal && shutdown_signal ()) {
    printf ("Writing dumpfile %s \n", fname);
    dump_function (fname, slist);
  }
  // ======================
  // Checkpoint on walltime
  // ======================
  else if (checkpoint_configuration.checkpoint_on_wall_time &&
           (long) (wall_time_current - wall_time_last) >
             checkpoint_configuration.checkpoint_on_wall_time_seconds) {
    wall_time_last = wall_time_current;
    printf ("Writing dumpfile %s \n", fname);
    dump_function (fname, slist);
  }
  // =============================
  // Checkpoint on simulation time
  // =============================
  else if (checkpoint_configuration.checkpoint_on_sim_time &&
           (double) (sim_time_current - sim_time_last) >
             checkpoint_configuration.checkpoint_on_sim_time_seconds) {
    sim_time_last = sim_time_current;
    printf ("Writing dumpfile %s \n", fname);
    dump_function (fname, slist);
  }
  // =============================
  // Checkpoint on simulation time
  // =============================
  else if (checkpoint_configuration.checkpoint_on_sim_iter &&
           (double) (sim_iter_current - sim_iter_last) >
             checkpoint_configuration.checkpoint_on_sim_iter_iterations) {
    sim_iter_last = sim_iter_current;
    printf ("Writing dumpfile %s \n", fname);
    dump_function (fname, slist);
  }

  if (shutdown_requested ()) {
    fprintf (stderr, "Shutdown requested (signal %d)", shutdown_signal ());
    return 1;
  } else {
    return 0;
  }
#endif
}

int restore_handler (const char* basename = NULL, scalar* slist = all) {
  install_shutdown_handlers ();

  char fname[64];
  if (!basename) {
    snprintf (fname, 64, "%s.dump", get_executable_name ());
  } else {
    snprintf (fname, 64, "%s.dump", basename);
  }

  if (file_exists (fname)) {
    printf ("Resuming simulation from checkpoint file %s\n", fname);
    restore_function (fname, slist);
    
    checkpoint_state.restored = true;
    checkpoint_state.restored_i = iter;
    checkpoint_state.restored_inext = inext;
    checkpoint_state.restored_t = t;
    checkpoint_state.restored_tnext = tnext;

    return 1;
  } else {
    return 0;
  }
}

// event stability (i++, first) {
//   if (checkpoint_state.restored) {
//     if (checkpoint_state.restored_i == i) {
//       double rt = checkpoint_state.restored_t;
//       double rtnext = checkpoint_state.restored_tnext;
//       dtnext(rtnext - rt);
//     }
//   }
// }

event checkpoint_event (i++, last) {
  return checkpoint_handler ();
}

