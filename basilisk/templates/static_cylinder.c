/* Immersed-boundary framework */
#include "library/ibm/IBMeshManager.h"

/* Multi-direct forcing solver */
#include "library/ibm/navier-stokes/centered-mdf.h"

/* Structural models */
#include "library/elff/elff.h"

/* Input/output headers */
#include "library/ibm/IBOutput.h"
#include "library/io/output-vtk.h"
#include "library/io/params/params-cli.h"

#include "tracer.h"
#include <stdio.h>

/* Default simulation parameters */

char *base_path = "static_cylinder_test";

int maxlevel = 10;
int minlevel = 6;
double Reynolds = 100.;
double U0 = 1.;
double L_fluid = 8;
double dt_fluid = 0.0005;
double t_end = 60.;

int ibmlevel = 10;
double R_cylinder = 0.15;
double alpha_cylinder = 1.0;
coord c_cylinder = {0.};

/* Derived parameters */

#define h_fluid (L_fluid / (1 << ibmlevel))
#define N_cylinder ((int)ceil(1.45 * pi * R_cylinder / (alpha_cylinder * h_fluid)))
#define ds_cylinder ((2. * pi * R_cylinder) / (double)N_cylinder)

/* Additional fields */

scalar f[];
scalar *tracers = {f};
face vector muv[];

/* Boundary conditions */

u.n[left] = dirichlet(U0);
p[left] = neumann(0.);
pf[left] = neumann(0.);
f[left] = dirichlet(y < 0);

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
  input_file_register_option("cylinder", R_cylinder, PARAM_VALUE_DOUBLE);
  input_file_register_option("cylinder", c_cylinder.x, PARAM_VALUE_DOUBLE);
  input_file_register_option("cylinder", c_cylinder.y, PARAM_VALUE_DOUBLE);
  input_file_register_option("cylinder", ibmlevel, PARAM_VALUE_INT);
  input_file_register_option("cylinder", alpha_cylinder, PARAM_VALUE_DOUBLE);

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
  origin(-1.85, -L0 / 2.);
#if TREE
  N = 1 << minlevel;
#else 
  N = 1 << ibmlevel;
#endif 
  DT = dt_fluid;
  mu = muv;

  ib_force_relaxation = 0.3;
  ib_richardson_iters = 4;

  display_control(Reynolds, 10, 1000);

  run();
}

event properties(i++) {
  foreach_face() muv.x[] = fm.x[] * 2 * R_cylinder * U0 / Reynolds;
}

event init(i = 0) {
  int m_id = ibmeshmanager_add_mesh();

  IBMeshModel cylinder_model = elff_pinned_rigid_body_circle_new(
      R_cylinder, N_cylinder, 1., c_cylinder, 0., 0);
  ibmeshmanager_set_model(m_id, cylinder_model);

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

int write_gp_file(const char * base_path) {
  if (pid() == 0) {
    if (!base_path) {
      fprintf(stderr, "warning: base_path is NULL; skipping gp output\n");
      return 0;
    }

    if (create_path(base_path) != 0) {
      fprintf(stderr, "warning: failed to create output path %s; skipping gp output\n", base_path);
      return 0;
    }

    char fname[4096];
    snprintf(fname, sizeof(fname), "%s/static-cylinder-drag.gp", base_path);

    FILE *fp = NULL;

    fp = fopen(fname, "w");

    const char *gp_str = 
      "if (!exists(\"csv_file\")) csv_file = \"cylinderstats.csv\"\n"
      "if (!exists(\"plot_output\")) plot_output = \"static-cylinder-drag.tex\"\n"
      "set datafile separator comma\n"
      "set terminal cairolatex pdf size 4.8in,3.2in color colortex font \",10\"\n"
      "set output plot_output\n"
      "set xlabel \"$t$\"\n"
      "set ylabel \"\"\n"
      "set grid\n"
      "set key off\n"
      "plot \\\n"
      "\tcsv_file every ::1 using 1:5 with linespoints title \"$C_l$\"\\\n";

    fclose(fp);
  }

  return 0;
}

event statsfile(i++) {
  double Cd = 0., Cl = 0.;
  double fx = 0., fy = 0.;
  foreach (reduction(+ : fx) reduction(+ : fy)) {
    fx += -ibmf.x[] * dv();
    fy += -ibmf.y[] * dv();
  }

  Cd = fx / (0.5 * sq(U0) * 2 * R_cylinder);
  Cl = fy / (0.5 * sq(U0) * 2 * R_cylinder);

  if (pid() == 0) {
    if (!base_path) {
      fprintf(stderr, "warning: base_path is NULL; skipping cylinder stats output\n");
      return 0;
    }
    if (create_path(base_path) != 0) {
      fprintf(stderr, "warning: failed to create output path %s; skipping cylinder stats output\n", base_path);
      return 0;
    }

    char fname[4096];
    snprintf(fname, sizeof(fname), "%s/cylinderstats.csv", base_path);

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
      fprintf(fp, "t,fx,fy,Cd,Cl\n");
    fprintf(fp, "%f,%f,%f,%f,%f\n", t, fx, fy, Cd, Cl);
    fclose(fp);
  }
}

event output(i += 200; t <= t_end)
// output(i += 1; i <= 5.)
{
  scalar omega[];
  vorticity(u, omega);
#if TREE
  output_hdf_htg(NULL, NULL, base_path);
#else
  output_hdf_imagedata(NULL, NULL, base_path);
#endif
  output_hdf_pd(NULL, NULL, base_path);
}

#if TREE
event adapt(i++) {
  adapt_wavelet_ibm({u, f}, (double[]){3e-2, 3e-2, 3e-2}, maxlevel, minlevel);
}
#endif

event checkpoint_event(i++, last) {
  return checkpoint_handler(t, i, base_path);
}
