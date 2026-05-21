#include "library/ibm/IBMeshManager.h"
#include "library/ibm/navier-stokes/centered-mdf.h"
#include "library/elff/elff.h"
#include "library/ibm/IBOutput.h"
#include "library/io/output-vtk.h"
#include "library/io/params/params-cli.h"

/* Default simulations parameters */

double dt_fluid = 0.0003;
double L_fluid = 16.;
double U0 = 1.;

int maxlevel = 10;
int minlevel = 6;
int ibmlevel = 11;

double huang_Re = 200.;
double huang_gamma = 0.001;
double huang_Fr = 0.5;
double huang_rho = 1.5;

double b_length = 1.;
int b_nodes = 65;
double b_r = 1e4;
double b_theta = (0.1 * pi);

char * base_path = "huang_output";

/* Derived parameters */

#define b_gravity (b_mu * huang_Fr)
#define b_EI huang_gamma
#define b_mu huang_rho

/* Additional fields */
face vector muv[]; 

/* Boundary conditions */

u.n[left] = dirichlet(U0);
u.t[left] = dirichlet(0.);
p[left] = neumann(0.);
pf[left] = neumann(0.);

u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);


int
main(int argc, char ** argv)
{
  /* Here we register runtime options for the simulation */
  input_file_register_option("basilisk.fluid", L_fluid, PARAM_VALUE_DOUBLE); 
  input_file_register_option("basilisk.fluid", dt_fluid, PARAM_VALUE_DOUBLE); 
  input_file_register_option("basilisk.fluid", minlevel, PARAM_VALUE_INT);
  input_file_register_option("basilisk.fluid", maxlevel, PARAM_VALUE_INT);
  input_file_register_option("basilisk.fluid", ibmlevel, PARAM_VALUE_INT);
  input_file_register_option("basilisk.output", base_path, PARAM_VALUE_STRING);
  input_file_register_option("huang", huang_Re, PARAM_VALUE_DOUBLE);
  input_file_register_option("huang", huang_gamma, PARAM_VALUE_DOUBLE);
  input_file_register_option("huang", huang_Fr, PARAM_VALUE_DOUBLE);
  input_file_register_option("huang", huang_rho, PARAM_VALUE_DOUBLE);
  input_file_register_option("beam", b_length, PARAM_VALUE_DOUBLE);
  input_file_register_option("beam", b_nodes, PARAM_VALUE_INT);
  input_file_register_option("beam", b_r, PARAM_VALUE_DOUBLE);

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
  origin(-2., -L0 / 2.);
#if TREE
  N = 1 << minlevel;
#else 
  N = 1 << ibmlevel;
#endif   
  DT = dt_fluid;
  mu = muv;
  display_control(huang_Re, 10, 1000);

  run();
}

/* Reynolds control */
event
properties(i++)
{
  foreach_face() muv.x[] = fm.x[] * b_length * U0 / huang_Re;
}

/* Beam setup */
event
init(i = 0)
{
  int m_id = ibmeshmanager_add_mesh();
  ib_euler_beam_bcs_t b_bcs = elff_euler_beam_bcs_theta_pin((coord){0});
  // IBMeshModel beam_model = elff_euler_beam_huang_new_theta(b_length, b_EI, b_mu, b_nodes, b_theta, b_bcs);
  IBMeshModel beam_model = elff_euler_beam_addm_new_theta(b_length, b_EI, b_mu, b_nodes, b_r, b_theta, b_bcs);
  ibmeshmanager_set_model(m_id, beam_model);

  foreach_ibnode_per_ibmesh()
  {
    mesh->depth = ibmlevel;
    node->depth = ibmlevel;
    ibval(gravity.x) = b_gravity;
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

event
logfile(i++) {
  if (pid() == 0)
    fprintf(stderr, "%d %g\n", i, t);
}

event
statsfile(i++)
{
  double x_tip = 0.;
  double y_tip = 0.;
  {
    IBNode* node = ibmm.pool.active.ptrs[0];
    x_tip = ibval(npos.y);
    y_tip = ibval(npos.y);
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
    if (i==0) {
      fprintf(fp, "i,t,x_tip,y_tip\n");
    }
    fprintf(fp, "%d,%g,%g,%g\n", i, t, x_tip, y_tip);
    fclose(fp);
  }
}


event
output(t += 0.05; t <= 50)
{
  scalar omega[];
  vorticity(u, omega);
#if TREE
  output_hdf_htg({omega,p}, {u,ibmf}, base_path);
#else
  output_hdf_imagedata({omega,p}, {u,ibmf}, base_path);
#endif
  output_hdf_pd(NULL, (IBvector[]){eulvel, nforce, nvel, {{-1}}}, base_path);
}

#if TREE
event
adapt(i++)
{
  adapt_wavelet_ibm({ u }, (double[]){ 5e-3, 5e-3 }, maxlevel, minlevel);
}
#endif

event
checkpoint_event(i++, last)
{
  return checkpoint_handler(t, i, base_path);
}
