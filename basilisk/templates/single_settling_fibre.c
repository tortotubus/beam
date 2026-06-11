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

int maxlevel = 10;
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
double banaei_gamma = 0.1;
double banaei_r = 0.1;
double banaei_Ga = 40.;

double b_length = 1.;
coord b_s0 = {-1. / 2., 0., 0.};
int b_nodes = 65;
double b_penalty = 1e3;
double b_theta = 0.00;
int b_pid = 0;

/* Nonphysical experiment controls */
double experiment_ib_force_relaxation = 0.4;
int experiment_ib_richardson_iters = 3;
int experiment_stats_interval = 1;
double experiment_output_interval = 0.05;
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
#define b_EI                                                                   \
  (banaei_gamma * banaei_r * b_submerged_weight_per_length * b_length *        \
   b_length * b_length)
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

int main(int argc, char **argv) {
  /* Here we register runtime options for the simulation */
  input_file_register_option("basilisk.fluid", L_fluid, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.fluid", dt_fluid, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.fluid", minlevel, PARAM_VALUE_INT);
  input_file_register_option("basilisk.fluid", maxlevel, PARAM_VALUE_INT);
  input_file_register_option("basilisk.fluid", ibmlevel, PARAM_VALUE_INT);
  input_file_register_option("basilisk.output", base_path, PARAM_VALUE_STRING);
  input_file_register_option_named("banaei", "rp", banaei_rp,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option_named("banaei", "gamma", banaei_gamma,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option_named("banaei", "r", banaei_r, PARAM_VALUE_DOUBLE);
  input_file_register_option_named("banaei", "Ga", banaei_Ga,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option("beam", b_length, PARAM_VALUE_DOUBLE);
  input_file_register_option("beam", b_nodes, PARAM_VALUE_INT);
  input_file_register_option_named("beam", "r_penalty", b_penalty,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option("beam", b_theta, PARAM_VALUE_DOUBLE);
  input_file_register_option("experiment", experiment_ib_force_relaxation,
                             PARAM_VALUE_DOUBLE);
  input_file_register_option("experiment", experiment_ib_richardson_iters,
                             PARAM_VALUE_INT);
  input_file_register_option("experiment", experiment_t_end,
                             PARAM_VALUE_DOUBLE);
  input_file_register_option("experiment", experiment_stats_interval,
                             PARAM_VALUE_INT);
  input_file_register_option("experiment", experiment_output_interval,
                             PARAM_VALUE_DOUBLE);

  /* Here we parse the options given through the command-line */
  int input_file_parse_result = input_file_parse_cli(argc, argv);
  if (input_file_parse_result != 0)
    return input_file_parse_result > 0 ? 0 : input_file_parse_result;
  else {
    input_file_print_options();
    input_file_apply_options();
  }

  /* Setting relevant parameters for basilisk */
  L0 = L_fluid;
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
  foreach_face() muv.x[] = fm.x[] * fluid_dynamic_viscosity;
}

/* Beam setup */
event init(i = 0) {
  ib_euler_beam_bcs_t b_bcs =
      elff_euler_beam_bcs_types(IB_EULER_BEAM_BC_FREE, IB_EULER_BEAM_BC_FREE);
  int m_id = ibmeshmanager_add_mesh();

  IBMeshModel beam_model = elff_euler_beam_addm_new_theta(
      b_length, b_EI, b_mu, b_nodes, b_penalty, b_theta, b_bcs, b_s0, b_pid);

  ibmeshmanager_set_model(m_id, beam_model);

  foreach_ibnode_per_ibmesh() {
    mesh->depth = ibmlevel;
    node->depth = ibmlevel;
    ibval(gravity.x) = 0.;
    ibval(gravity.y) = b_gravity;
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

    fprintf(fp, "%d,%g,%g,%g,%g,%g%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
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
