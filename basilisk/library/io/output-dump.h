#include "library/io/output-common.h"
#include "output.h"

#define OUTPUT_DUMP_H 1

@ifndef dump_function 
@define dump_function (fname, slist) dump ((fname), (slist), NULL, false, true) 
@endif

  
@ifndef restore_function 
@define restore_function (fname, slist) restore ((fname), (slist), NULL) 
@endif
 

static struct {
  bool checkpoint_on_wall_time;
  double checkpoint_on_wall_time_seconds;
  bool checkpoint_on_sim_time;
  double checkpoint_on_sim_time_seconds;
  bool checkpoint_on_sim_iter;
  int checkpoint_on_sim_iter_iterations;
  bool checkpoint_on_signal;
  int checkpoint_on_signals[6];
} checkpoint_configuration = {
  .checkpoint_on_signal = true,
  .checkpoint_on_signals = {SIGINT, SIGTERM},
  .checkpoint_on_wall_time = true,
  .checkpoint_on_wall_time_seconds = 300,
  .checkpoint_on_sim_time = false,
  .checkpoint_on_sim_iter = false
};

#include "library/io/input-signals.h"

int checkpoint_handler (double sim_time_current = t,
                        int sim_iter_current = i,
                        const char* basename = NULL,
                        scalar* slist = all) {
  char fname[1024];
  char pname[1024];
  char tname[1024];

  if (!basename) {
    snprintf (pname, 1024, "%s", get_unique_basename ());
    snprintf (fname, 1024, "%s/dump", get_unique_basename ());
    snprintf (tname, 1024, "%s/dump.dt", get_unique_basename ());
  } else {
    snprintf (pname, 1024, "%s", basename);
    snprintf (fname, 1024, "%s/dump", basename);
    snprintf (tname, 1024, "%s/dump.dt", basename);
  }

  assert (!create_path (pname));

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
    printf ("\nWriting dumpfile %s \n", fname);
    dump_function (fname, slist);
#ifdef CKPT_TIMESTEP
    write_checkpoint_sidecar(tname, sim_time_current, dt, sim_iter_current, inext, tnext);
#endif
  }

  // ======================
  // Checkpoint on walltime
  // ======================
  else if (checkpoint_configuration.checkpoint_on_wall_time &&
           (long) (wall_time_current - wall_time_last) >
             checkpoint_configuration.checkpoint_on_wall_time_seconds) {
    wall_time_last = wall_time_current;
    printf ("\nWriting dumpfile %s \n", fname);
    dump_function (fname, slist);
#ifdef CKPT_TIMESTEP
    write_checkpoint_sidecar(tname, sim_time_current, dt, sim_iter_current, inext, tnext);
#endif
  }

  // =============================
  // Checkpoint on simulation time
  // =============================
  else if (checkpoint_configuration.checkpoint_on_sim_time &&
           (double) (sim_time_current - sim_time_last) >
             checkpoint_configuration.checkpoint_on_sim_time_seconds) {
    sim_time_last = sim_time_current;
    printf ("\nWriting dumpfile %s \n", fname);
    dump_function (fname, slist);
#ifdef CKPT_TIMESTEP
    write_checkpoint_sidecar(tname, sim_time_current, dt, sim_iter_current, inext, tnext);
#endif
  }

  // =============================
  // Checkpoint on simulation time
  // =============================
  else if (checkpoint_configuration.checkpoint_on_sim_iter &&
           (double) (sim_iter_current - sim_iter_last) >
             checkpoint_configuration.checkpoint_on_sim_iter_iterations) {
    sim_iter_last = sim_iter_current;
    printf ("\nWriting dumpfile %s \n", fname);
    dump_function (fname, slist);
#ifdef CKPT_TIMESTEP
    write_checkpoint_sidecar(tname, sim_time_current, dt, sim_iter_current, inext, tnext);
#endif
  }

  if (shutdown_requested ()) {
    fprintf (
      stderr, "warning: Shutdown requested (signal %d)\n", shutdown_signal ());
    return 1;
  } else {
    return 0;
  }
#endif
}

int restore_handler (const char* basename = NULL, scalar* slist = all) {
  install_shutdown_handlers ();

  char fname[1024];
  char pname[1024];
  char tname[1024];

  if (!basename) {
    snprintf (pname, 1024, "%s", get_unique_basename ());
    snprintf (fname, 1024, "%s/dump", get_unique_basename ());
    snprintf (tname, 1024, "%s/dump.dt", get_unique_basename ());
  } else {
    snprintf (pname, 1024, "%s", basename);
    snprintf (fname, 1024, "%s/dump", basename);
    snprintf (tname, 1024, "%s/dump.dt", basename);
  }

  if (file_exists (fname)) {
    printf ("Resuming simulation from checkpoint file %s\n", fname);
    restore_function (fname, slist);
#ifdef CKPT_TIMESTEP
    read_checkpoint_sidecar(tname, &cpsc);
#endif
    return 1;
  } else {
    return 0;
  }
}
