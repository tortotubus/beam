#include "grid/quadtree.h"
// #include "grid/multigrid.h"

#include "library/ibm/IBMeshManager.h"
#include "library/ibm/navier-stokes/unserious/centered-split-rich.h"
#include "tracer.h"
#include "library/io/output-vtk.h"
#include "library/io/output-dump.h"

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

const int maxlevel = 10;
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
  L0 = L_fluid;
  origin(-1.85, -L0 / 2.);

#if TREE
  // N = 1 << minlevel;
  N = 1 << lvl_circ;
#else
  N = 1 << lvl_circ;
#endif

  dV_circ = L_fluid / (1 << lvl_circ) * ds_circ;
  mu = muv;

  display_control(Reynolds, 10, 1000);

  DT = 0.0005;

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
init_ib(i = 0)
{
  int new_id = ibmeshmanager_add_mesh();
  ibmeshmanager_add_nodes(new_id, N_circ);
  foreach_ibnode_per_ibmesh() {
    node->depth = lvl_circ; 
    mesh->depth = lvl_circ;
    ibval(nweight) = ds_circ;
  }
  foreach_ibnode_per_ibmesh() {
    coord cpos = circle(node_id, N_circ, c_circ, R_circ);
    foreach_dimension() {
      ibval(npos.x) = cpos.x;
    }
  }
}

event init(t = 0) foreach () u.x[] = U0;

event
logfile(i++)
{
  fprintf(stderr, "%d %g %d %d %d\n", i, t, mgp.i, mgu_a.i, mgu_b.i);
}

event
statsfile(i++)
{

  coord f = { 0 };
  double Cd = 0., Cl = 0.;
  double fx = 0., fy = 0.;
  foreach (reduction(+:fx) reduction(+:fy)) {
    fx += -ibmf.x[] * dv();
    fy += -ibmf.y[] * dv();
  }

  Cd = fx / (0.5 * sq(U0) * 2 * R_circ);
  Cl = fy / (0.5 * sq(U0) * 2 * R_circ);

  if (pid() == 0) {
    if (i == 0) {
      FILE* fp = fopen("cylinderstats.txt", "w");
      fprintf(fp, "");
      fclose(fp);
    }

    FILE* fp = fopen("cylinderstats.txt", "a");

    fprintf(fp, "%f %f %f %f %f\n", t, fx, fy, Cd, Cl);
    fclose(fp);
  }
}

event
output(i += 200; t <= 60.)
// output(i += 1; i <= 5.)
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
    { u, f }, (double[]){ 3e-2, 3e-2, 3e-2 }, maxlevel, minlevel);
}
#endif
