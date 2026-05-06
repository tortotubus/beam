#include "library/ibm/IBMeshManager.h"
#include "library/ibm/navier-stokes/centered-mdf.h"
#include "library/elff/elff.h"
#include "library/ibm/IBOutput.h"
#include "library/io/output-vtk.h"
#include "library/io/params/params-cli.h"

#include "lambda2.h"

/* Default simulation parameters */

char *base_path = "static_sphere_test";

int maxlevel = 10;
int minlevel = 5;
double Reynolds = 100.;
double U0 = 1.;
double L_fluid = 40;
double dt_fluid = 0.0001;
double t_end = 60.;

int ibmlevel = 10;
double D_sphere = 1.;
double alpha_sphere = 2.1;
// double alpha_sphere = 1.8;

coord c_sphere = {0.};

/* Derived parameters */

#define h_fluid (L_fluid / (1 << ibmlevel))
#define R_sphere (D_sphere / 2.)
#define ds_sphere (alpha_sphere * h_fluid)
#define dA_sphere ((sqrt(3.) / 2.) * sq(ds_sphere))
#define N_sphere ((int)ceil(pi * sq(D_sphere) / dA_sphere))
// #define N_sphere 10

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
  input_file_register_option("sphere", D_sphere, PARAM_VALUE_DOUBLE);
  input_file_register_option("sphere", c_sphere.x, PARAM_VALUE_DOUBLE);
  input_file_register_option("sphere", c_sphere.y, PARAM_VALUE_DOUBLE);
  input_file_register_option("sphere", c_sphere.z, PARAM_VALUE_DOUBLE);
  input_file_register_option("sphere", ibmlevel, PARAM_VALUE_INT);
  input_file_register_option("sphere", alpha_sphere, PARAM_VALUE_DOUBLE);

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
  origin(-5 * D_sphere, -L0 / 2., -L0 / 2);
  N = 1 << minlevel;
  // DT = dt_fluid;
  mu = muv;
  display_control(Reynolds, 10, 1000);

  run();
}

event properties(i++) {
  foreach_face() muv.x[] = fm.x[] * D_sphere * U0 / Reynolds;
}

event init(i = 0) {
  int m_id = ibmeshmanager_add_mesh();

  IBMeshModel sphere_model =
      elff_pinned_rigid_body_sphere_new(R_sphere, N_sphere, 1., c_sphere, 1., {0.}, 0);
  ibmeshmanager_set_model(m_id, sphere_model);

  foreach_ibnode_per_ibmesh() {
    node->depth = ibmlevel;
    mesh->depth = ibmlevel;
  }

  if (!restore_handler(base_path)) {
    foreach ()
      u.x[] = U0;
  } else {
  }
}

event logfile(i++) {
  fprintf(stderr, "%d %g %d %d %d\n", i, t, mgp.i, mgu_a.i, mgu_b.i);
}

event output(i += 1; t<=t_end) 
{ 
  scalar l2[];
  lambda2(u,l2);

  scalar * slist  = {l2, p};
  vector * vlist = {u,ibmf};

#if TREE
  output_hdf_htg(slist, vlist, base_path);
#else
  output_hdf_imagedata(slist, vlist, base_path);
#endif
  output_hdf_pd(NULL, NULL, base_path);
}

#if TREE
event adapt(i++) {
  adapt_wavelet_ibm({u}, (double[]){0.05, 0.05, 0.05}, maxlevel, minlevel);
}
#endif

event checkpoint_event(i++, last) {
  return checkpoint_handler(t, i, base_path);
}
