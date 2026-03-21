#include "grid/quadtree.h"
// #include "grid/multigrid.h"

#include "library/ibm/IBMeshManager.h"
// #include "library/ibm/navier-stokes/centered-dlmfd.h"
#include "library/ibm/navier-stokes/unserious/centered-split-rich.h"
#include "library/elff/elff.h"
#include "library/io/output-vtk.h"

face vector muv[];

#define L_fluid 16.
#define maxlevel 10
#define minlevel 6
#define ibmlevel 11
#define h_fluid (L_fluid / (1 << maxlevel))

#define huang_Re 200.
#define huang_gamma 0.001
#define huang_Fr 0.5
#define huang_rho 1.5

#define b_length 1.
#define b_EI huang_gamma
#define b_mu huang_rho
#define b_nodes ((int)65)               //((b_length / b_ds) + 1.))
#define b_ds (b_length / (b_nodes - 1)) // (2. * h_fluid)
#define b_r 1e3
// #define b_r 1e4
#define b_theta (0.1 * pi)
#define b_gravity (1.)

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

// u.t[top] = dirichlet(U0);
// u.n[top] = dirichlet(0.);
// u.t[bottom] = dirichlet(U0);
// u.n[bottom] = dirichlet(0.);

u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);

event
init_ib(i = 0)
{
  // ib_force_relaxation = 0.1;
  int m_id = ibmeshmanager_add_mesh();
  IBMeshModel beam_model =
    elff_beam_new_theta(b_length, b_EI, b_mu, b_nodes, b_r, b_theta);
  ibmeshmanager_set_model(m_id, beam_model);
  foreach_ibnode_per_ibmesh()
  {
    mesh->depth = ibmlevel;
    node->depth = ibmlevel;
    ibval(gravity.x) = b_gravity;
    ibval(nweight) = b_ds;
    if (node_id == 0 || node_id == mesh->nodes.size - 1)
      ibval(nweight) *= 0.5;
  }
}

event init(t = 0) foreach () u.x[] = U0;

event
logfile(i++)
{
  double y_tip = 0.;
  {
    IBNode *node = ibmm.pool.active.ptrs[0];
    y_tip = ibval(npos.y);
  }
  // fprintf(stderr, "%d %g %d\n", i, t, cgiter);
  if (pid() == 0)
  fprintf(stderr, "%d %g %g\n", i, t, y_tip);
}

scalar omega[];

event
output(t+=0.05; t <= 50)
{

  vorticity(u, omega);
#if TREE
  output_hdf_htg();
#else
  output_hdf_imagedata();
#endif
  output_hdf_pd();
}

#if TREE
event
adapt(i++)
{
  adapt_wavelet_ibm({ u }, (double[]){ 3e-3, 3e-3 }, maxlevel, minlevel);
}
#endif
