#include "grid/quadtree.h"
#include "library/ibm/IBAdapt.h"
// #include "library/ibm/IBAdapt2.h"
#include "library/ibm/IBMeshManager.h"
#include "library/io/output-vtk.h"

#include "run.h" 

scalar dummy[];

coord
circle(int n, int N, coord centre, double radius)
{
  double rad = 2. * pi * ((double)n / (double)N);
  coord c = centre;
  c.x += radius * cos(rad);
  c.y += radius * sin(rad);
  return c;
}

const int maxlevel = 6;
const int minlevel = 6;

const double L_fluid = 8;
const int lvl_circ = 10;
const double R_circ = 0.15;
const coord c_circ = { 0 };
const double h_fluid = L_fluid / (1 << lvl_circ);
const int N_circ = (int)(1.5 * pi * R_circ / h_fluid);
const double ds_circ = (double)(2 * pi * R_circ) / (N_circ);
double dV_circ = 0;


int
main()
{
  periodic(right);
  L0 = L_fluid;
  origin(-L0 / 2, -L0 / 2.);

#if TREE
  N = 1 << lvl_circ;
#else
  N = 1 << lvl_circ;
#endif

  display_control(Reynolds, 10, 1000);
 

  run();
}

event
init_ib(i = 0)
{
  ibmeshmanager_init(0);

  int new_id = ibmeshmanager_add_mesh();
  ibmeshmanager_add_nodes(new_id, N_circ);
  foreach_ibnode_per_ibmesh()
  {
    mesh->depth = lvl_circ;
    node->depth = lvl_circ;
  }
  foreach_ibnode_per_ibmesh()
  {
    foreach_dimension()
      ibval(npos.x)  = circle(node_id, N_circ, c_circ, R_circ).x;
  }

  ibmeshmanager_update_pid();

  adapt_wavelet_ibm(NULL, NULL, maxlevel);
}



event
output(i += 1; i <= 2000)
{ 
  if (pid() == 0)
    printf("%d\n", i);
#if TREE
  output_hdf_htg();
#else
  output_hdf_imagedata();
#endif
  output_hdf_pd();
}

event
update_ib(i++)
{
  foreach_ibnode()
  {
    ibval(npos.x) += 0.005;
  }
  ibmeshmanager_update_pid();
}

#if TREE
event
adapt(i++)
{
  foreach() {
    dummy[] = 0.0;
  }

  boundary({dummy});

  adapt_wavelet_ibm(NULL, NULL, maxlevel);
}
#endif
