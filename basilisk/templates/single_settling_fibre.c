#include "library/ibm/IBMeshManager.h"
#include "library/ibm/navier-stokes/centered-mdf.h"
#include "library/elff/elff.h"
#include "library/ibm/IBOutput.h"
#include "library/io/output-vtk.h"
#include "library/io/params/params-cli.h"

#include "lambda2.h"

/* Default simulations parameters */

double dt_fluid = 0.001; // 5e-2,5e-3
double L_fluid = 16.;
double U0 = 0.;

int maxlevel = 11;
int minlevel = 4;
int ibmlevel = 11;

/*
 * Banaei et al. (2020) nondimensional groups:
 *   rp    = L/d
 *   gamma = EI/(rho_l r g L^3), with rho_l = Delta rho A_f
 *   r     = Delta rho/rho_0
 *   Ga    = sqrt(r g L^3)/nu
 */

double banaei_rp = 30.;
double banaei_gamma = 1.0;
double banaei_r = 0.1;
double banaei_Ga = 40.;

double b_length = 1.;
coord b_s0 = {-1. / 2., 0., 0.};
int b_nodes = 67;
double b_theta = 0.00;

/* Nonphysical experiment controls */
double experiment_ib_force_relaxation = 0.8;
int experiment_ib_richardson_iters = 2;
int experiment_huang_implicit_bending = 1;

int experiment_stats_interval = 1;
double experiment_output_interval = 1.0;
int experiment_output_iter = 1;

double experiment_t_end = 500.0;

char *base_path = "single_settling_fibre_output";

/* Derived parameters */
#define b_rho_0 (1.)
#define b_g (1.)

#define b_diameter (b_length / banaei_rp)
#define b_area (pi * b_diameter * b_diameter / 4.)
#define b_linear_density_difference (banaei_r * b_rho_0 * b_area)
#define b_mu ((1. + banaei_r) * b_rho_0 * b_area)
#define b_submerged_weight_per_length (b_linear_density_difference * b_g)
#define b_EI (banaei_gamma * banaei_r * b_submerged_weight_per_length * b_length * b_length * b_length)
#define b_gravity -(b_submerged_weight_per_length)

#define fluid_nu (sqrt(banaei_r * b_g * b_length * b_length * b_length) / banaei_Ga)
#define fluid_dynamic_viscosity (b_rho_0 * fluid_nu)
#define fluid_velocity_scale (sqrt(banaei_r * b_g * b_length))

/* Additional fields */

face vector muv[];

/* Quiescent tank boundary conditions */

u.n[left] = dirichlet(0.);
u.t[left] = dirichlet(0.);
p[left] = neumann(0.);
pf[left] = neumann(0.);

u.n[right] = dirichlet(0.);
u.t[right] = dirichlet(0.);
p[right] = neumann(0.);
pf[right] = neumann(0.);

u.n[top] = dirichlet(0.);
u.t[top] = dirichlet(0.);
p[top] = neumann(0.);
pf[top] = neumann(0.);

u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
p[bottom] = neumann(0.);
pf[bottom] = neumann(0.);

#include "single_settling_fibre_parser.h"

static bool settling_xml = false;
static SingleSettlingFibreCase settling_configuration;
face vector settling_alpha[];
scalar settling_rho[];

static bool single_settling_fibre_is_xml (int argc, char **argv) {
  if (argc != 2) return false;
  size_t n = strlen (argv[1]);
  return n >= 4 && !strcmp (argv[1] + n - 4, ".xml");
}

static int single_settling_fibre_load_xml (const char *path) {
  if (single_settling_fibre_case_init (path, &settling_configuration) != 0)
    return 1;
  settling_xml = true;
  L_fluid = settling_configuration.basilisk.domain_size;
  dt_fluid = settling_configuration.basilisk.time_step;
  minlevel = settling_configuration.basilisk.min_level;
  maxlevel = settling_configuration.basilisk.max_level;
  ibmlevel = settling_configuration.basilisk.immersed_boundary_level;
  base_path = settling_configuration.basilisk.output_directory;
  experiment_t_end = settling_configuration.basilisk.end_time;
  experiment_output_interval = settling_configuration.basilisk.output_interval;
  alpha = settling_alpha;
  rho = settling_rho;
  return 0;
}

static void single_settling_fibre_register_xml_models () {
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
      ibval (gravity.x) = settling_configuration.basilisk.uniform_force_density[0];
      ibval (gravity.y) = settling_configuration.basilisk.uniform_force_density[1];
#if dimension > 2
      ibval (gravity.z) = settling_configuration.basilisk.uniform_force_density[2];
#endif
    }
}

static void single_settling_fibre_xml_stats (int iteration, double time) {
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

int main(int argc, char **argv) {
  if (single_settling_fibre_is_xml (argc, argv)) {
    if (single_settling_fibre_load_xml (argv[1])) return 1;
  } else {
  /* Here we register runtime options for the simulation */
  input_file_register_option("basilisk.fluid", L_fluid, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.fluid", dt_fluid, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.fluid", minlevel, PARAM_VALUE_INT);
  input_file_register_option("basilisk.fluid", maxlevel, PARAM_VALUE_INT);
  input_file_register_option("basilisk.fluid", ibmlevel, PARAM_VALUE_INT);
  input_file_register_option("basilisk.output", base_path, PARAM_VALUE_STRING);
  input_file_register_option_named("banaei", "rp", banaei_rp, PARAM_VALUE_DOUBLE);
  input_file_register_option_named("banaei", "gamma", banaei_gamma, PARAM_VALUE_DOUBLE);
  input_file_register_option_named("banaei", "r", banaei_r, PARAM_VALUE_DOUBLE);
  input_file_register_option_named("banaei", "Ga", banaei_Ga,  PARAM_VALUE_DOUBLE);
  input_file_register_option("beam", b_length, PARAM_VALUE_DOUBLE);
  input_file_register_option("beam", b_nodes, PARAM_VALUE_INT);
  input_file_register_option("beam", b_theta, PARAM_VALUE_DOUBLE);
  input_file_register_option("experiment", experiment_ib_force_relaxation, PARAM_VALUE_DOUBLE);
  input_file_register_option("experiment", experiment_ib_richardson_iters, PARAM_VALUE_INT);
  input_file_register_option("experiment", experiment_huang_implicit_bending, PARAM_VALUE_INT);
  input_file_register_option("experiment", experiment_t_end, PARAM_VALUE_DOUBLE);
  input_file_register_option("experiment", experiment_stats_interval, PARAM_VALUE_INT);
  input_file_register_option("experiment", experiment_output_interval, PARAM_VALUE_DOUBLE);

  /* Here we parse the options given through the command-line */
  int input_file_parse_result = input_file_parse_cli(argc, argv);
  if (input_file_parse_result != 0)
    return input_file_parse_result > 0 ? 0 : input_file_parse_result;
  else {
    input_file_print_options();
    input_file_apply_options();
  }

  }

  /* Setting relevant parameters for basilisk */
  L0 = L_fluid;
  if (settling_xml)
    origin(settling_configuration.basilisk.origin[0], settling_configuration.basilisk.origin[1],
           settling_configuration.basilisk.origin[2]);
  else
    origin(-L_fluid / 2., -L_fluid / 2., -L_fluid / 2.);
#if TREE
  N = 1 << minlevel;
#else
  N = 1 << ibmlevel;
#endif
  DT = dt_fluid;
  mu = muv;
  ib_force_relaxation = experiment_ib_force_relaxation;
  ib_richardson_iters = experiment_ib_richardson_iters;
  display_control(banaei_Ga, 1, 1000);

  run();
}

event logfile(i++) {
  if (pid() == 0)
    fprintf(stderr, "%d %g\n", i, t);
}

/* Reynolds control */
event properties(i++) {
  foreach_face() {
    muv.x[] = fm.x[] * (settling_xml ? settling_configuration.basilisk.dynamic_viscosity : fluid_dynamic_viscosity);
    if (settling_xml) settling_alpha.x[] = fm.x[] / settling_configuration.basilisk.fluid_density;
  }
  if (settling_xml)
    foreach () settling_rho[] = cm[] * settling_configuration.basilisk.fluid_density;
}

/* Beam setup */
event init(i = 0) {
  if (settling_xml) {
    single_settling_fibre_register_xml_models ();
  } else {

  ib_euler_beam_bcs_t b_bcs = elff_euler_beam_bcs_types(IB_EULER_BEAM_BC_FREE, IB_EULER_BEAM_BC_FREE);
  int m_id = ibmeshmanager_add_mesh();
  IBMeshModel beam_model = elff_euler_beam_huang_new_theta(
    b_length, b_EI, b_mu, b_nodes, b_theta, b_bcs, b_s0);
  ib_euler_beam_huang_set_implicit_bending((ib_euler_beam_huang_t) elff_force_coupled_handle (beam_model),
                                           experiment_huang_implicit_bending);
  ibmeshmanager_set_model(m_id, beam_model);

  foreach_ibnode_per_ibmesh() {
    mesh->depth = ibmlevel;
    node->depth = ibmlevel;
    ibval(gravity.x) = 0.;
    ibval(gravity.y) = b_gravity;
  }

  }

  if (!restore_handler(base_path)) {
    foreach ()
      foreach_dimension() u.x[] = 0.;
#if TREE
    adapt_wavelet_ibm(NULL, NULL, 0, 1, all, true);
#endif
  } else {
    // if restored
  }
}

event marchetti_csv(i += experiment_stats_interval; t <= experiment_t_end) {
  if (settling_xml) { single_settling_fibre_xml_stats (i, t); return 0; }

  coord pos_first;
  coord pos_last;
  coord pos_middle;
  coord pos_min_vert = {0., Y0 + L0, 0.};

  double delta_first;
  double delta_last;

  foreach_ibnode() {
    if (node_id == 0) {
      foreach_dimension() { pos_first.x = ibval(npos.x); }
    } else if (node_id == ibmm.pool.active.size - 1) {
      foreach_dimension() { pos_last.x = ibval(npos.x); }
    } else if (node_id == b_nodes / 2) {
      foreach_dimension() { pos_middle.x = ibval(npos.x); }
    }
    if (pos_min_vert.y > ibval(npos.y)) {
      foreach_dimension() { pos_min_vert.x = ibval(npos.x); }
    }
  }

  delta_first = (pos_first.y - pos_min_vert.y) / (0.5 * b_length);
  delta_last = (pos_last.y - pos_min_vert.y) / (0.5 * b_length);

  if (pid() == 0) {
    fprintf(stderr, "%d %g %g %g %g\n", i, t, pos_min_vert.y, delta_first,
            delta_last);

    FILE *fp = NULL;

    if (!base_path) {
      fprintf(stderr, "warning: base_path is NULL; skipping tip output\n");
      return 0;
    }

    create_path(base_path);
    char fname[4096];
    snprintf(fname, sizeof(fname), "%s/banaei-marchetti-validation.csv",
             base_path);

    fp = fopen(fname, i == 0 ? "w" : "a");

    if (!fp) {
      fprintf(stderr, "warning: failed to open %s for append\n", fname);
      return 0;
    }

    if (i == 0) {
      fprintf(fp, "i,t,banaei_Ga,banaei_gamma,banaei_r,banaei_rp,B,x0,xm,xmy,xf,"
                  "y0,ym,ymy,yf,z0,zm,zmy,zf,delta_0,delta_f\n");
    }

    fprintf(fp, "%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
            i, t, banaei_Ga, banaei_gamma, banaei_r, banaei_rp, (1./(banaei_gamma * banaei_r)), pos_first.x,
            pos_middle.x, pos_min_vert.x, pos_last.x, pos_first.y, pos_middle.y,
            pos_min_vert.y, pos_last.y, pos_first.z, pos_middle.z,
            pos_min_vert.z, pos_last.z, delta_first, delta_last);
    fclose(fp);
  }
}

event output(t += experiment_output_interval; t <= experiment_t_end) {
  scalar l2[], omega_z[];
  lambda2(u, l2);
  vorticity(u, omega_z);

#if TREE
  output_hdf_htg({l2, omega_z, p}, {u, ibmf}, base_path);
#else
  output_hdf_imagedata({l2, omega_z, p}, {u, ibmf}, base_path);
#endif
  output_hdf_pd(NULL, (IBvector[]){eulvel, nforce, nvel, {{-1}}}, base_path);
}

#if TREE
event adapt(i++) {
  double adapt_rel_u_tol = 1e-4;
  double val = adapt_rel_u_tol * fluid_velocity_scale;
  adapt_wavelet_ibm({u}, (double[]){val, val, val}, maxlevel, minlevel);
}
#endif

event checkpoint_event(i++, last) {
  return checkpoint_handler(t, i, base_path);
}

event elff_cleanup(i = end, last) {
  ibmeshmanager_free();
  elff_deinit();
}

event stop(t = experiment_t_end) {
  if (pid () == 0) fprintf (stderr, "finished single-fibre settling at t=%g i=%d\n", t, i);
  return 1;
}
