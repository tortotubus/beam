#include "library/ibm/IBMeshManager.h"
#include "library/ibm/navier-stokes/centered-mdf.h"
#include "library/elff/elff.h"
#include "library/ibm/IBOutput.h"
#include "library/io/output-vtk.h"
#include "library/io/params/params-cli.h"

#include "lambda2.h"

/* Default simulation parameters */

char *base_path = "static_fibre_test";

int maxlevel = 11;
int minlevel = 5;
double Reynolds = 0.05;
double U0 = 1.;
double L_fluid = 16.;
double dt_fluid = 0.0005;
double dt_output = 0.5;
double t_end = 60.;

int ibmlevel = 11;
double L_fibre = 1.;
double D_fibre = 1. / 30.;
double alpha_fibre = 1.5;
double theta_fibre = .5 * pi;

double experiment_ib_force_relaxation = 0.9;
int experiment_ib_richardson_iters = 3;

coord c_fibre = {0.};

/* Derived parameters */

#define h_fluid (L_fluid / (1 << ibmlevel))
#define ds_fibre (alpha_fibre * h_fluid)
#define N_fibre ((int)ceil(L_fibre / ds_fibre) + 1)
#define A_fibre (L_fibre * D_fibre)
#define aspect_ratio_fibre (L_fibre / D_fibre)
#define delta_perp_fibre                                                      \
  (0.839 + 0.185 / aspect_ratio_fibre +                                      \
   0.233 / sq(aspect_ratio_fibre))
#define Cd_sbt_fibre                                                          \
  (8. * pi / (Reynolds * (log(aspect_ratio_fibre) + delta_perp_fibre)))

/* Additional fields */

face vector muv[];

/* Boundary conditions */

u.n[left] = dirichlet(U0);
p[left] = neumann(0.);
pf[left] = neumann(0.);

u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);

int main(int argc, char **argv) {
  /* Here we register runtime options for the simulation */
  input_file_register_option("basilisk.fluid", L_fluid, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.fluid", dt_fluid, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.fluid", t_end, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.fluid", minlevel, PARAM_VALUE_INT);
  input_file_register_option("basilisk.fluid", maxlevel, PARAM_VALUE_INT);
  input_file_register_option("basilisk.fluid", Reynolds, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.output", base_path, PARAM_VALUE_STRING);
  input_file_register_option("basilisk.output", dt_output, PARAM_VALUE_DOUBLE);
  input_file_register_option("fibre", L_fibre, PARAM_VALUE_DOUBLE);
  input_file_register_option("fibre", D_fibre, PARAM_VALUE_DOUBLE);
  input_file_register_option("fibre", c_fibre.x, PARAM_VALUE_DOUBLE);
  input_file_register_option("fibre", c_fibre.y, PARAM_VALUE_DOUBLE);
  input_file_register_option("fibre", c_fibre.z, PARAM_VALUE_DOUBLE);
  input_file_register_option("fibre", theta_fibre, PARAM_VALUE_DOUBLE);
  input_file_register_option("fibre", ibmlevel, PARAM_VALUE_INT);
  input_file_register_option("fibre", alpha_fibre, PARAM_VALUE_DOUBLE);
  input_file_register_option("experiment", experiment_ib_force_relaxation,
                             PARAM_VALUE_DOUBLE);
  input_file_register_option("experiment", experiment_ib_richardson_iters,
                             PARAM_VALUE_INT);

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
  origin(-5 * L_fibre, -L0 / 2., -L0 / 2);
  N = 1 << minlevel;
  DT = dt_fluid;
  mu = muv;
  ib_force_relaxation = experiment_ib_force_relaxation;
  ib_richardson_iters = experiment_ib_richardson_iters;
  display_control(Reynolds, 10, 1000);

  run();
}

event properties(i++) {
  foreach_face() muv.x[] = fm.x[] * D_fibre * U0 / Reynolds;
}

event init(i = 0) {
  int m_id = ibmeshmanager_add_mesh();

  const double q_w = cos(theta_fibre / 2.);
  const coord q_xyz = {0., 0., sin(theta_fibre / 2.)};
  IBMeshModel fibre_model = elff_pinned_rigid_body_fibre_new(
      L_fibre, D_fibre, N_fibre, 1., c_fibre, q_w, q_xyz, 0);
  ibmeshmanager_set_model(m_id, fibre_model);

  foreach_ibnode_per_ibmesh() {
    node->depth = ibmlevel;
    mesh->depth = ibmlevel;
  }

  if (!restore_handler(base_path)) {
    foreach ()
      u.x[] = U0;

#if TREE
    adapt_wavelet_ibm(NULL, NULL, 0, 1, all, true);
#endif
  } else {
  }
}

event logfile(i++) {
  fprintf(stderr, "%d %g %d %d %d\n", i, t, mgp.i, mgu_a.i, mgu_b.i);
}

event csvfile(i++) {
  static int have_prev_Cd = 0;
  static double prev_Cd = 0.;

  double Cd = 0., Cl = 0.;
  double fx = 0., fy = 0., fz = 0.;
  foreach (reduction(+ : fx) reduction(+ : fy) reduction(+ : fz)) {
    fx += -ibmf.x[] * dv();
    fy += -ibmf.y[] * dv();
    fz += -ibmf.z[] * dv();
  }

  double qA = 0.5 * sq(U0) * A_fibre;
  Cd = fx / qA;
  Cl = sqrt(sq(fy) + sq(fz)) / qA;
  if (pid() == 0) {
    double dCd_rel = have_prev_Cd && fabs(prev_Cd) > 0.
                         ? fabs(Cd - prev_Cd) / fabs(prev_Cd)
                         : 0.;
    double Cl_over_Cd = fabs(Cd) > 0. ? Cl / fabs(Cd) : 0.;

    fprintf(stderr,
            "stats i=%d t=%g fx=%g Cd=%g Cd_target=%g Cd/Cd_target=%g "
            "dCd_rel=%g Cl=%g Cl/Cd=%g\n",
            i,
            t,
            fx,
            Cd,
            Cd_sbt_fibre,
            Cd / Cd_sbt_fibre,
            dCd_rel,
            Cl,
            Cl_over_Cd);
    fflush(stderr);
    prev_Cd = Cd;
    have_prev_Cd = 1;

    if (!base_path) {
      fprintf(stderr, "warning: base_path is NULL; skipping fibre stats output\n");
      return 0;
    }
    if (create_path(base_path) != 0) {
      fprintf(stderr, "warning: failed to create output path %s; skipping fibre stats output\n", base_path);
      return 0;
    }

    char fname[4096];
    snprintf(fname, sizeof(fname), "%s/fibre_stats_Re_%f.csv", base_path, Reynolds);

    FILE *fp = NULL;
    if (i == 0)
      fp = fopen(fname, "w");
    else
      fp = fopen(fname, "a");
    if (!fp) {
      fprintf(stderr, "warning: failed to open %s for write\n", fname);
      return 0;
    }
    if (i == 0)
      fprintf(fp,
              "t,fx,fy,fz,Cd,Cl,Cd_sbt,Cd_over_Cd_sbt,h,ds,"
              "support_over_D,ds_over_D,aspect_ratio\n");
    fprintf(fp,
            "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
            t,
            fx,
            fy,
            fz,
            Cd,
            Cl,
            Cd_sbt_fibre,
            Cd / Cd_sbt_fibre,
            h_fluid,
            ds_fibre,
            2. * h_fluid / D_fibre,
            ds_fibre / D_fibre,
            aspect_ratio_fibre);
    fclose(fp);
  }
}

event output(t += dt_output; t <= t_end) {
  scalar l2[];
  lambda2(u, l2);

  scalar *slist = {l2, p};
  vector *vlist = {u, ibmf};

#if TREE
  output_hdf_htg(slist, vlist, base_path);
#else
  output_hdf_imagedata(slist, vlist, base_path);
#endif
  output_hdf_pd(NULL, NULL, base_path);
}

#if TREE
event adapt(i++) {
  adapt_wavelet_ibm({u}, (double[]){0.02, 0.02, 0.02}, maxlevel, minlevel);
}
#endif

event checkpoint_event(i++, last) {
  return checkpoint_handler(t, i, base_path);
}
