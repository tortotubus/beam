#include "library/ibm/IBMeshManager.h"
#include "library/ibm/navier-stokes/centered-ibm.h"
#include "library/ibm/IBViscosityRatio.h"
#include "library/elff/elff.h"
#include "library/ibm/IBOutput.h"
#include "library/io/output-vtk.h"
#include "library/io/params/params-cli.h"

#ifndef RBC_L0
#define RBC_L0 4.
#endif
#ifndef RBC_RADIUS
#define RBC_RADIUS 1.
#endif
#ifndef RBC_SHEAR_RATE
#define RBC_SHEAR_RATE 1.
#endif
#ifndef RBC_RE
#define RBC_RE 0.01
#endif
#ifndef RBC_CA
#define RBC_CA 0.1
#endif
#ifndef RBC_ND_EB
#define RBC_ND_EB 0.01
#endif
#ifndef RBC_AREA_DILATATION_MODULUS
#define RBC_AREA_DILATATION_MODULUS 50.
#endif
#ifndef RBC_C0
#define RBC_C0 (-2.09/RBC_RADIUS)
#endif

char *base_path = "rbc_shear_test";

int maxlevel = 6;
int minlevel = 4;
int ibmlevel = 6;
int capsule_refinements = 4;
double t_end = 30. / RBC_SHEAR_RATE;
double dt_fluid = 1.e-4;
double dt_output = 0.5;
double output_interval = 10.;
double Reynolds = RBC_RE;
double shear_rate = RBC_SHEAR_RATE;
double radius = RBC_RADIUS;
double capillary = RBC_CA;
double area_dilatation_modulus = RBC_AREA_DILATATION_MODULUS;
double nondim_bending_modulus = RBC_ND_EB;
double reference_curvature = RBC_C0;
double outer_viscosity = 1.;
double inner_viscosity = 5.;
coord capsule_center = {0., 0., 0.};

#define elastic_modulus (radius * shear_rate / capillary)
#define bending_modulus (nondim_bending_modulus * elastic_modulus * sq(radius))
#define fluid_density (Reynolds / (shear_rate * sq(radius)))

face vector muv[];
scalar capsule_viscosity_indicator[];
const scalar myrho[] = fluid_density;
const face vector myalpha[] = {1. / fluid_density,
                               1. / fluid_density,
                               1. / fluid_density};

ib_capsule_t capsule = NULL;
int capsule_mesh_id = -1;

u.n[left] = dirichlet(0.);
u.n[right] = dirichlet(0.);
u.r[left] = dirichlet(0.);
u.r[right] = dirichlet(0.);
u.t[left] = dirichlet(-shear_rate*x);
u.t[right] = dirichlet(-shear_rate*x);
uf.n[left] = dirichlet(0.);
uf.n[right] = dirichlet(0.);

int main(int argc, char **argv) {
  input_file_register_option("basilisk.fluid", t_end, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.fluid", dt_fluid, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.fluid", minlevel, PARAM_VALUE_INT);
  input_file_register_option("basilisk.fluid", maxlevel, PARAM_VALUE_INT);
  input_file_register_option("basilisk.fluid", Reynolds, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.output", base_path, PARAM_VALUE_STRING);
  input_file_register_option("basilisk.output", dt_output, PARAM_VALUE_DOUBLE);
  input_file_register_option("basilisk.output", output_interval,
                             PARAM_VALUE_DOUBLE);
  input_file_register_option("capsule", ibmlevel, PARAM_VALUE_INT);
  input_file_register_option("capsule", capsule_refinements, PARAM_VALUE_INT);
  input_file_register_option("capsule", radius, PARAM_VALUE_DOUBLE);
  input_file_register_option("capsule", capillary, PARAM_VALUE_DOUBLE);
  input_file_register_option("capsule", area_dilatation_modulus,
                             PARAM_VALUE_DOUBLE);
  input_file_register_option("capsule", nondim_bending_modulus,
                             PARAM_VALUE_DOUBLE);
  input_file_register_option("capsule", reference_curvature,
                             PARAM_VALUE_DOUBLE);
  input_file_register_option("capsule", outer_viscosity, PARAM_VALUE_DOUBLE);
  input_file_register_option("capsule", inner_viscosity, PARAM_VALUE_DOUBLE);
  input_file_register_option("capsule", capsule_center.x, PARAM_VALUE_DOUBLE);
  input_file_register_option("capsule", capsule_center.y, PARAM_VALUE_DOUBLE);
  input_file_register_option("capsule", capsule_center.z, PARAM_VALUE_DOUBLE);

  int input_file_parse_result = input_file_parse_cli(argc, argv);
  if (input_file_parse_result != 0)
    return input_file_parse_result > 0 ? 0 : input_file_parse_result;

  input_file_print_options();
  input_file_apply_options();

  L0 = RBC_L0;
  origin(-0.5 * L0, -0.5 * L0, -0.5 * L0);
  periodic(top);
  periodic(front);

  N = 1 << minlevel;
  DT = dt_fluid;
  rho = myrho;
  alpha = myalpha;
  mu = muv;
  stokes = true;
  ib_capsule_viscosity_set_wall_boundary_conditions(
    capsule_viscosity_indicator);

  run();
}

event properties(i++) {
  if (capsule && capsule_mesh_id >= 0) {
    ib_capsule_construct_viscosity_indicator(capsule_viscosity_indicator,
                                             capsule, capsule_mesh_id);
    ib_capsule_set_viscosity_from_indicator(muv, capsule_viscosity_indicator,
                                            outer_viscosity, inner_viscosity);
  } else {
    foreach_face()
      muv.x[] = outer_viscosity * fm.x[];
  }
}

event init(i = 0) {
  capsule_mesh_id = ibmeshmanager_add_mesh();

  capsule = ib_capsule_biconcave_new(
    radius, elff_vertex_from_coord(capsule_center), capsule_refinements,
    elastic_modulus);
  ib_capsule_set_skalak_law(capsule, elastic_modulus,
                            area_dilatation_modulus);
  ib_capsule_set_helfrich_bending_law(capsule, bending_modulus);
  ib_capsule_set_constant_reference_curvature(capsule, reference_curvature);

  elff_runtime_register((ib_model_t)capsule, 0);
  IBMeshModel capsule_model =
    elff_velocity_coupled_model(capsule, elff_capsule_destroy);
  ibmeshmanager_set_model(capsule_mesh_id, capsule_model);

  foreach_ibnode_per_ibmesh() {
    node->depth = ibmlevel;
    mesh->depth = ibmlevel;
  }

  if (!restore_handler(base_path)) {
    foreach() {
      u.x[] = 0.;
      u.y[] = -shear_rate*x;
      u.z[] = 0.;
    }

#if TREE
    adapt_wavelet_ibm(NULL, NULL, 0, 1, all, true);
#endif
  }
}

event logfile(i++) {
  fprintf(stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);
}

event csvfile(i += 10) {
  double fx = 0., fy = 0., fz = 0.;
  foreach(reduction(+ : fx) reduction(+ : fy) reduction(+ : fz)) {
    fx += ibmf.x[] * dv();
    fy += ibmf.y[] * dv();
    fz += ibmf.z[] * dv();
  }

  if (pid() == 0) {
    if (!base_path || create_path(base_path) != 0)
      return 0;

    char fname[4096];
    snprintf(fname, sizeof(fname), "%s/rbc_shear_stats.csv", base_path);
    FILE *fp = i == 0 ? fopen(fname, "w") : fopen(fname, "a");
    if (!fp)
      return 0;
    if (i == 0)
      fprintf(fp, "t,fx,fy,fz\n");
    fprintf(fp, "%g,%g,%g,%g\n", t, fx, fy, fz);
    fclose(fp);
  }
}

event output(i += output_interval; t <= t_end) {
  scalar *slist = {p, capsule_viscosity_indicator};
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
  adapt_wavelet_ibm({u}, (double[]){0.01, 0.01, 0.01},
                    maxlevel, minlevel);
}
#endif

event checkpoint_event(i++, last) {
  return checkpoint_handler(t, i, base_path);
}

event end(t = t_end) {
  return 0.;
}
