// #include "grid/multigrid.h"

#include "library/ibm/IBMeshManager.h"

#include "library/ibm/navier-stokes/centered-split.h"
#include "tracer.h"

#include "library/ibm/IBOutput.h"
#include "library/io/output-vtk.h"

coord
circle(int n, int N, coord centre, double radius)
{
  double rad = 2. * pi * ((double)n / (double)N);
  coord c = centre;
  c.x += radius * cos(rad);
  c.y += radius * sin(rad);
  return c;
}

scalar f[];
scalar* tracers = { f };
face vector muv[];

double Reynolds = 100.;
double U0 = 1.;

const int maxlevel = 9;
const int minlevel = 4;

const int N_circ = 90;
const double R_circ = 0.15;
const int L_circ = 11;
const coord c_circ = { 0 };

int
main()
{
  L0 = 16;
  origin(-1.85, -L0/2.);
  N = 1 << 5;
  mu = muv;

  display_control(Reynolds, 10, 1000);

  // DT = 0.003;

  run();
}

event
properties(i++)
{
  foreach_face() muv.x[] = fm.x[] * 2 * R_circ * U0 / Reynolds;
}

u.n[left] = dirichlet(U0);
p[left] = neumann(0.);
pf[left] = neumann(0.);
f[left] = dirichlet(y < 0);

u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);

event
init(t = 0)
{
  int new_id = ibmeshmanager_add_mesh();
  ibmeshmanager_add_nodes(new_id, N_circ);

  double ds = (2 * pi * R_circ) / N_circ;

  foreach_ibmesh()
  {
    mesh->refinement_level = L_circ;
  }

  foreach_ibnode_per_ibmesh()
  {
    node->lagpos = circle(node_id, N_circ, c_circ, R_circ);
    // node->dV = L0 / (1 << L_circ) * ds;
  }

  foreach ()
    u.x[] = U0;
}

event logfile(i++) { 
  fprintf(stderr, "%d %g %d %d %d\n", i, t, mgp.i, mgu_a.i, mgu_b.i); 
}

event
output(t += 0.25; t <= 30.)
// output(i += 1; t <= 30.)
{
  scalar omega[];
  vorticity(u, omega);
  output_hdf_htg({ p, omega, f }, { u, ibmf });
  output_ibnodes("uhllman_cylinder", i, t);
}

event
adapt(i++)
{
  adapt_wavelet_ibm({ u, f }, (double[]){ 3e-2, 3e-2, 3e-2 }, maxlevel, minlevel);
}
