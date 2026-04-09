
#include "grid/quadtree.h" 
#include "library/ibm/IBMeshManager.h"
#include "library/ibm/navier-stokes/centered-mdf.h"
#include "library/elff/elff.h"
#include "library/ibm/IBOutput.h"
#include "library/io/output-vtk.h"

#define basenamestr "huang_two_fibres_Re300_Mu0001_T0003_d0.3"

#define L_fluid 16.
#define maxlevel 10
#define minlevel 6
#define ibmlevel 11
#define h_fluid (L_fluid / (1 << maxlevel))

#define huang_Re 300.
#define huang_gamma 0.001
#define huang_Fr 0.5
#define huang_rho 1.5

#define b_d 0.15
#define b_length 1.
#define b_EI huang_gamma
#define b_mu huang_rho
#define b_nodes ((int)65)               
#define b_ds (b_length / (b_nodes - 1)) 
#define b_r 1e3
#define b_theta (0.1 * pi)
#define b_gravity (b_mu * huang_Fr)


face vector muv[];
int Reynolds = huang_Re;
double U0 = 1.;

int
main()
{
  L0 = L_fluid;
  origin(-2., -L0 / 2.);
  N = 1 << ibmlevel;
  mu = muv;
  display_control(Reynolds, 10, 1000);
  DT = 0.0003;

  install_shutdown_handlers();
  run();
}

event
properties(i++)
{
  foreach_face() muv.x[] = fm.x[] * b_length * U0 / Reynolds;
}

u.n[left] = dirichlet(U0);
u.t[left] = dirichlet(0.);
p[left] = neumann(0.);
pf[left] = neumann(0.);

u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);

event
init(i = 0)
{
  int m_id1 = ibmeshmanager_add_mesh();
  IBMeshModel beam_model1 = elff_euler_beam_new_theta(b_length, b_EI, b_mu, b_nodes, b_r, b_theta, {0.,b_d,0.});
  ibmeshmanager_set_model(m_id1, beam_model1);

  int m_id2 = ibmeshmanager_add_mesh();
  IBMeshModel beam_model2 = elff_euler_beam_new_theta(b_length, b_EI, b_mu, b_nodes, b_r, b_theta, {0.,-b_d,0.});
  ibmeshmanager_set_model(m_id2, beam_model2);

  foreach_ibnode_per_ibmesh()
  {
    mesh->depth = ibmlevel;
    node->depth = ibmlevel;
    ibval(gravity.x) = b_gravity;
    ibval(nweight) = b_ds;
    if (node_id == 0 || node_id == mesh->nodes.size - 1)
      ibval(nweight) *= 0.5;
  }

  if (!restore_handler(basenamestr)) {
    foreach ()
      u.x[] = U0;
  } else {
  }
}

event
logfile(i++) {
  if (pid() == 0)
    fprintf(stderr, "[basilisk] %d %g\n", i, t);
}

event
statsfile(i++)
{
  double x_tip[2];
  double y_tip[2];

  foreach_ibmesh()
  {
    IBNode* node = mesh->nodes.ptrs[0];
    x_tip[mesh_id] = ibval(npos.x);
    y_tip[mesh_id] = ibval(npos.y);
  }

  if (pid() == 0) {
    FILE *fp = NULL;
    create_path(basenamestr);
    char fname[4096];
    snprintf(fname, sizeof(fname), "%s/tips.txt", basenamestr);
    fp = fopen(fname, "a");
    fprintf(fp, "%d %g %g %g %g %g\n", i, t, x_tip[0], y_tip[0], x_tip[1], y_tip[1]);
    fprintf(stderr, "%d %g %g %g %g %g\n", i, t, x_tip[0], y_tip[0], x_tip[1], y_tip[1]);
    fclose(fp);
  }
}

event
output(t += 0.05; t <= 50)
{
  scalar omega[];
  vorticity(u, omega);
#if TREE
  output_hdf_htg({omega,p}, {u,ibmf}, basenamestr);
#else
  output_hdf_imagedata({omega,p}, {u,ibmf}, basenamestr);
#endif
  output_hdf_pd(NULL, (IBvector[]){eulvel, nforce, nvel, {{-1}}}, basenamestr);
}


#if TREE
event
adapt(i++)
{
  adapt_wavelet_ibm({ u }, (double[]){ 3e-3, 3e-3 }, maxlevel, minlevel);
}
#endif

event
checkpoint_event(i++, last)
{
  return checkpoint_handler(t, i, basenamestr);
}
