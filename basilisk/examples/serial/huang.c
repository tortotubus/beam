#include "grid/quadtree.h"
// #include "grid/multigrid.h"

#include "library/ibm/IBMeshManager.h"

// #include "library/ibm/navier-stokes/centered-dlmfd.h"
#include "library/ibm/navier-stokes/unserious/centered-split-rich.h"

#include "library/elff/elff.h"

#include "library/io/output-vtk.h"

face vector muv[];


#define L_fluid 8.
#define maxlevel 10
#define minlevel 6
#define ibmlevel 10
#define h_fluid (L_fluid / (1 << maxlevel))

#define huang_Re 200.
#define huang_gamma 0.001
#define huang_Fr 0.5
#define huang_rho 1.5

#define b_length 1.
#define b_EI huang_gamma
#define b_mu huang_rho
#define b_nodes ((int) 64 ) //((b_length / b_ds) + 1.))
#define b_ds (b_length / (b_nodes-1)) // (2. * h_fluid)
#define b_r 1e2
// #define b_r 1e4
#define b_theta (0.1 * pi)
#define b_gravity (b_mu * huang_Fr)

int Reynolds = huang_Re;
double U0 = 1.;

int
main()
{

  printf("b_nodes: %d\n", b_nodes);
 

  L0 = L_fluid;
  origin(-L0 / 4., -L0 / 2.);

#if TREE
  N = 1 << minlevel;
#else
  N = 1 << ibmlevel;
#endif

  mu = muv;
  display_control(Reynolds, 10, 1000);
  DT = 0.0003;

  // cgiter_max = 10000;
  // cgtol = 1e-7;

  run();
}

event
properties(i++)
{
  foreach_face() muv.x[] = fm.x[]*b_length*U0 / Reynolds;
}


u.n[left] = dirichlet(U0);
u.t[left] = dirichlet(0.);
p[left] = neumann(0.);
pf[left] = neumann(0.);

u.t[top] = dirichlet(U0);
u.n[top] = dirichlet(0.);
u.t[bottom] = dirichlet(U0);
u.n[bottom] = dirichlet(0.);

u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);


event
init_ib(i = 0)
{
  // ib_force_relaxation = 0.1;
  int m_id = ibmeshmanager_add_mesh();
  IBMeshModel beam_model = elff_beam_new_theta(b_length, b_EI, b_mu, b_nodes, b_r, b_theta);
  ibmeshmanager_set_model(m_id, beam_model);
  foreach_ibnode_per_ibmesh() {
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
  // fprintf(stderr, "%d %g %d\n", i, t, cgiter);
  fprintf(stderr, "%d %g\n", i, t);
}

event
output(t+=0.05; t <= 50)
{ 

  scalar omega[];
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
  adapt_wavelet_ibm(
    { u}, (double[]){ 3e-2, 3e-2}, maxlevel, minlevel);
}
#endif
