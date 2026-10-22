/** Shared XML-driven settling experiment for one or many fibres.
 * Banaei parameters and dimensional coefficients are evaluated in the case XML.
 */
#ifndef IBMESHMANAGER_POOL_SIZE_BYTES
#define IBMESHMANAGER_POOL_SIZE_BYTES (1 << 24)
#endif
#include "library/ibm/navier-stokes/centered-mdf.h"
#include "library/io/xml/xml.h"
#include "lambda2.h"
#if dimension != 3
#error "settling_fibre requires grid/octree.h"
#endif
#include "settling_fibre_parser.h"

static SettlingFibreCase settling_configuration;
static int stats_interval = 1;
static char *base_path;
scalar l2[], omega_z[];

/** @brief Refresh derived fields immediately before fluid output. */
static void settling_fibre_prepare_output () {
  lambda2 (u, l2);
  vorticity (u, omega_z);
}

static void settling_fibre_stats (int iteration, double time) {
  // Position and velocity sums are weighted by nodal measure, not replicas.
  size_t count = (size_t) ibmm.nm * 7;
  double *values = calloc (count, sizeof (double));
  if (!values) { fprintf (stderr, "Unable to allocate fibre statistics.\n"); exit (1); }
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

int main (int argc, char **argv) {
  if (argc != 2) {
    if (pid() == 0) fprintf (stderr, "Usage: settling_fibre case.xml\n");
    return 1;
  }
  if (settling_fibre_case_init (argv[1], &settling_configuration)) return 1;
  BasiliskXMLCase *configuration = &settling_configuration.basilisk;
  stats_interval = settling_configuration.stats_interval;
  base_path = configuration->output_directory;
  configuration->prepare_fluid_output = settling_fibre_prepare_output;
  ib_model_view_array_t models = elff_models ();
  if (pid() == 0)
    fprintf (stderr, "settling fibres=%zu periodic rho=%g mu=%g gravity_y=%g\n",
             models.count, configuration->fluid_density, configuration->dynamic_viscosity,
             configuration->force_density[1]);
  ib_force_relaxation = configuration->force_relaxation;
  ib_richardson_iters = configuration->force_iterations;
  int status = generic_ibm_case_run (configuration);
  elff_deinit ();
  basilisk_xml_case_destroy (configuration);
  return status;
}

event logfile (i += stats_interval; t <= generic_case.time_end) {
  settling_fibre_stats (i, t);
}
