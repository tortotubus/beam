#include "library/elff/elff.h"
#include "library/ibm/IBMeshManager.h"
#include "library/ibm/IBOutput.h"
#include "library/ibm/navier-stokes/centered-mdf.h"
#include "library/io/output-vtk.h"

#include "flapping_fibre_parser.h"
#include "fibre_case_runtime.h"
 
/*
  The default simulation parameters correspond to Figure 13 of
  https://doi.org/10.1016/j.jcp.2007.07.002
*/

static FlappingFibreConfiguration case_configuration = {0};
static char *base_path = NULL;
static double U0 = 0.;
static double fluid_density = 1.;
static double fluid_viscosity = 0.;
static double output_interval = 0.05;
static double end_time = 0.;
static int maxlevel = 0;
static int minlevel = 0;
static int ibmlevel = 0;

/* Additional fields */
face vector muv[], alphav[];
scalar rhov[];

/* Boundary conditions */

u.n[left] = dirichlet(U0);
u.t[left] = dirichlet(0.);
p[left] = neumann(0.);
pf[left] = neumann(0.);

u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);

int main(int argc, char **argv) {
  const char *case_path = argc > 1 ? argv[1] : "flapping_fibre.xml";
  if (fibre_case_init(case_path, &case_configuration) != 0)
    return 1;

  /* Setting relevant parameters for basilisk */
  L0 = case_configuration.domain_size;
  origin(case_configuration.origin[0], case_configuration.origin[1]);
  U0 = case_configuration.inlet_velocity[0];
  fluid_density = case_configuration.fluid_density;
  fluid_viscosity = case_configuration.dynamic_viscosity;
  minlevel = case_configuration.min_level;
  maxlevel = case_configuration.max_level;
  ibmlevel = case_configuration.immersed_boundary_level;
  base_path = case_configuration.output_directory;
  output_interval = case_configuration.output_interval;
  end_time = case_configuration.end_time;
#if TREE
  N = 1 << minlevel;
#else
  N = 1 << ibmlevel;
#endif
  DT = case_configuration.time_step;
  mu = muv;
  alpha = alphav;
  rho = rhov;

  run();
}

/* Reynolds control */
event properties(i++) {
  foreach_face() {
    muv.x[] = fm.x[] * fluid_viscosity;
    alphav.x[] = fm.x[] / fluid_density;
  }
  foreach ()
    rhov[] = cm[] * fluid_density;
}

static bool recover_missing_ib_force = false;
static int flapping_recover_missing_ib () {
  if (!elff_checkpoint_restored) return -1;
  recover_missing_ib_force = true;
  if (pid () == 0)
    fprintf (stderr, "Missing dump.ib: restored ELFF geometry; reconstructing coupling force from fluid state. This is an approximate restart.\n");
  return 0;
}

/* Beam setup */
event init(i = 0) {
  ib_model_view_array_t models = elff_models();
  assert(models.count > 0);
  for (size_t model_index = 0; model_index < models.count; ++model_index) {
    int mesh_index = ibmeshmanager_add_mesh();
    assert(mesh_index == (int)model_index);
    int rc = elff_register_model(mesh_index, (int)model_index,
                                 models.data[model_index], ibmlevel);
    assert(rc == 0);
  }

  foreach_ibnode_per_ibmesh() {
    node->depth = mesh->depth;
    if (case_configuration.has_uniform_force_density) {
      ibval(gravity.x) = case_configuration.uniform_force_density[0];
      ibval(gravity.y) = case_configuration.uniform_force_density[1];
#if dimension > 2
      ibval(gravity.z) = case_configuration.uniform_force_density[2];
#endif
    }
  }

  ib_checkpoint_recover_missing = flapping_recover_missing_ib;
  if (!restore_handler(base_path)) {
    foreach ()
      u.x[] = U0;
#if TREE
    adapt_wavelet_ibm(NULL, NULL, 0, 1, all, true);
#endif
  } else {
  }
}

/* Sparse nodes may materialize after adaptation; refresh their body force
   immediately before the solver consumes it. */
event advance_interface(i++) {
  if (recover_missing_ib_force) {
    // Recover the missing nodal force without applying a second fluid impulse.
    vector saved_u[], saved_ibmf[];
    foreach ()
      foreach_dimension () {
        saved_u.x[] = u.x[];
        saved_ibmf.x[] = ibmf.x[];
      }
    event ("interface_force_richardson");
    foreach ()
      foreach_dimension () {
        u.x[] = saved_u.x[];
        ibmf.x[] = saved_ibmf.x[];
      }
    boundary ((scalar*) {u, ibmf});
    recover_missing_ib_force = false;
  }
  if (case_configuration.has_uniform_force_density)
    foreach_ibnode () {
      ibval(gravity.x) = case_configuration.uniform_force_density[0];
      ibval(gravity.y) = case_configuration.uniform_force_density[1];
#if dimension > 2
      ibval(gravity.z) = case_configuration.uniform_force_density[2];
#endif
    }
}

event logfile(i++) {
  if (pid() == 0)
    fprintf(stderr, "%d %g\n", i, t);
}

event statsfile(i++) {
  double x_tip = 0.;
  double y_tip = 0.;
  {
    foreach_ibnode () {
      if (node->mesh_gid == 0 && node->node_lid == 0 && node->pid == pid ()) {
        x_tip = ibval(npos.x);
        y_tip = ibval(npos.y);
      }
    }
#if _MPI
    MPI_Allreduce(MPI_IN_PLACE, &x_tip, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &y_tip, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  }

  if (pid() == 0) {
    FILE *fp = NULL;
    if (!base_path) {
      fprintf(stderr, "warning: base_path is NULL; skipping tip output\n");
      return 0;
    }
    create_path(base_path);
    char fname[4096];
    snprintf(fname, sizeof(fname), "%s/tip.csv", base_path);
    fp = fopen(fname, "a");
    if (!fp) {
      fprintf(stderr, "warning: failed to open %s for append\n", fname);
      return 0;
    }
    if (i == 0) {
      fprintf(fp, "i,t,x_tip,y_tip\n");
    }
    fprintf(fp, "%d,%g,%g,%g\n", i, t, x_tip, y_tip);
    fclose(fp);
  }
}

event output(t += output_interval; t <= end_time) {
  scalar omega[];
  vorticity(u, omega);
#if TREE
  output_hdf_htg({omega, p}, {u, ibmf}, base_path);
#else
  output_hdf_imagedata({omega, p}, {u, ibmf}, base_path);
#endif
  output_hdf_pd(NULL, (IBvector[]){eulvel, nforce, nvel, {{-1}}}, base_path);
}

#if TREE
event adapt(i++) {
  adapt_wavelet_ibm({u}, (double[]){5e-3, 5e-3}, maxlevel, minlevel);
}
#endif

event checkpoint_event(i++, last) {
  return checkpoint_handler(t, i, base_path);
}

event elff_cleanup(i = end, last) {
  ibmeshmanager_free();
  elff_deinit();
}
