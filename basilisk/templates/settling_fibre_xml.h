#pragma once
#include "fibre_case_runtime.h"

static bool settling_xml = false;
static FibreCaseConfiguration settling_configuration;
face vector settling_alpha[];
scalar settling_rho[];

static bool settling_is_xml (int argc, char **argv) {
  if (argc != 2) return false;
  size_t n = strlen (argv[1]);
  return n >= 4 && !strcmp (argv[1] + n - 4, ".xml");
}

static int settling_load_xml (const char *path, bool single) {
  if (fibre_case_init (path, &settling_configuration) != 0)
    return 1;
  ib_model_view_array_t models = elff_models ();
  if (!models.count || (single && models.count != 1)) {
    fprintf (stderr, "Expected %s in <ELFF><Models>.\n",
             single ? "exactly one fibre" : "at least one fibre");
    elff_deinit ();
    return 1;
  }
  settling_xml = true;
  L_fluid = settling_configuration.domain_size;
  dt_fluid = settling_configuration.time_step;
  minlevel = settling_configuration.min_level;
  maxlevel = settling_configuration.max_level;
  ibmlevel = settling_configuration.immersed_boundary_level;
  base_path = settling_configuration.output_directory;
  experiment_t_end = settling_configuration.end_time;
  experiment_output_interval = settling_configuration.output_interval;
  alpha = settling_alpha;
  rho = settling_rho;
  return 0;
}

static void settling_register_xml_models () {
  ib_model_view_array_t models = elff_models ();
  for (size_t m = 0; m < models.count; ++m) {
    int mesh_id = ibmeshmanager_add_mesh ();
    if (models.data[m].handle)
      ib_euler_beam_huang_set_implicit_bending (
        (ib_euler_beam_huang_t) models.data[m].handle,
        experiment_huang_implicit_bending);
    if (elff_register_model (mesh_id, (int)m, models.data[m], ibmlevel) != 0) {
      fprintf (stderr, "Cannot register XML fibre %zu.\n", m);
      exit (1);
    }
  }
}

// Sparse nodes can appear after grid adaptation or model advancement.
event advance_interface (i++) {
  if (settling_xml)
    foreach_ibnode () {
      ibval (gravity.x) = settling_configuration.uniform_force_density[0];
      ibval (gravity.y) = settling_configuration.uniform_force_density[1];
#if dimension > 2
      ibval (gravity.z) = settling_configuration.uniform_force_density[2];
#endif
    }
}

static void settling_xml_stats (int iteration, double time) {
  // Position and velocity sums are weighted by nodal measure, not replicas.
  size_t count = (size_t) ibmm.nm * 7;
  double *values = calloc (count, sizeof (double));
  foreach_ibnode_per_ibmesh () {
#if _MPI
    if (node->pid != pid ()) continue;
#endif
    double *v = values + 7*mesh_id;
    double w = ibval (nweight);
    v[0] += w;
    v[1] += w*ibval (npos.x); v[2] += w*ibval (npos.y);
    v[4] += w*ibval (nvel.x); v[5] += w*ibval (nvel.y);
#if dimension > 2
    v[3] += w*ibval (npos.z); v[6] += w*ibval (nvel.z);
#endif
  }
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, values, (int)count, MPI_DOUBLE, MPI_SUM,
                 MPI_COMM_WORLD);
#endif
  if (pid () == 0) {
    create_path (base_path);
    char filename[4096];
    snprintf (filename, sizeof(filename), "%s/fibre-centres.csv", base_path);
    FILE *fp = fopen (filename, iteration == 0 ? "w" : "a");
    if (!fp) { perror (filename); exit (1); }
    if (iteration == 0)
      fprintf (fp, "i,t,fibre_id,center_x,center_y,center_z,velocity_x,velocity_y,velocity_z,measure\n");
    for (int m = 0; m < ibmm.nm; ++m) {
      double *v = values + 7*m;
      double inv = v[0] > 0. ? 1./v[0] : 0.;
      fprintf (fp, "%d,%.17g,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
               iteration, time, m, v[1]*inv, v[2]*inv, v[3]*inv,
               v[4]*inv, v[5]*inv, v[6]*inv, v[0]);
    }
    fclose (fp);
    fprintf (stderr, "settling i=%d t=%g fibres=%d\n", iteration, time, ibmm.nm);
  }
  free (values);
}
