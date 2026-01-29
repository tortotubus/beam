

#define DEBUG 0

#define LEVEL 9
#define DT_SAVE 0.1
#define T_END 0.1

scalar f[];
scalar* tracers = { f };
double Reynolds = 200.;
int maxlevel = 9;
int Nf = 1 << 9;
face vector muv[];

double U0 = 1.;

#include "grid/quadtree.h"

#include "library/ibm/IBcentered.h"
#include "library/ibm/IBoutput.h"
#include "library/output_htg.h"

#include "embed.h"
#include "tracer.h"

int
main()
{
  L0 = 8;
  origin(-2, -4.);
  mu = muv;
  init_grid(Nf);
  DT = 0.0005;

  int n_meshes = 1;
  ibnode_list_init(100);
  ib_mesh_manager_init(n_meshes);

  // double f_dx = L0 / Nf;

  double b_length = 1;
  // // double b_EI = 0.0001;
  // // double b_mu = 1.5;
  int b_nodes = 65;
  double b_ds = b_length / b_nodes;
  // // double b_r = 1e4;
  double b_theta = 0.1 * pi;

  ib_mesh_init_nn(ibmesh(0), b_nodes);

  foreach_ibnode_mesh(ibmesh(0))
  {
    double radius = (node_index * b_ds);
    node->lagpos.x = radius * cos(b_theta);
    node->lagpos.y = radius * sin(b_theta);
    node->weight = b_ds;
  }

  run();
}

u.n[left] = dirichlet(U0);
p[left] = neumann(0.);
pf[left] = neumann(0.);
f[left] = dirichlet(y < 0);

u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);

event
init(i = 0)
{

  foreach ()
    u.x[] = cs[] ? U0 : 0.;
}

event
properties(i++)
{
  foreach_face() muv.x[] = 1. / Reynolds;
}

event
progress_output(i += 20)
{
  if (pid() == 0) {
    foreach_ibmesh()
    {
      coord centroid = ib_mesh_get_centroid(mesh);
      printf("[Mesh #%d, Time: %f]: (%f,%f,%f)\n",
             mesh_index,
             t,
             centroid.x,
             centroid.y,
             centroid.z);
    }
  }
}

event
vtk(t += DT_SAVE; t <= T_END)
{
  scalar omega[];
  vorticity(u, omega);

  scalar proc[];
  foreach () {
    proc[] = pid();
  }

  output_ibnodes("linetest_serial_nodes", i, t);
  output_hdf_htg({ omega, p, stencils, f, proc },
                 { u, forcing },
                 "linetest_serial_fluid");
}

event
adapt(i++)
{
  // ib_mesh_manager_tag_stencil();
  adapt_wavelet({ cs, u, f, stencils },
                (double[]){ 1e-2, 3e-2, 3e-2, 3e-2, 3e-2 },
                maxlevel,
                4);
}

event
end(t = T_END)
{
  return 0;
}