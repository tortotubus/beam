#define DT_SAVE 0.1
#define T_END 10.0

const double huang_Re = 200;
const double huang_gamma = 0.001;
const double huang_Fr = 0.5;
const double huang_rho = 1.5;

const double U0 = 1.;
const double Reynolds = 200.;

const double b_length = 1;
const double b_EI = huang_gamma;
const double b_mu = huang_rho;
const int b_nodes = 50;
const double b_ds = b_length / b_nodes;
const double b_r = 1e4;
const double b_theta = 0.1 * pi;
const double b_gravity = b_mu * huang_Fr;



const int maxlevel = 9;
const int Nf = 1 << maxlevel;

face vector muv[];


#include "grid/quadtree.h"

#include "interface/ibm/beam.h"
#include "library/ibm/uzawa2.h"
#include "library/output_htg.h"

#include "embed.h"

int
main()
{
  L0 = 8;
  origin(-L0 / 4., -L0 / 2.);
  init_grid(Nf);

  mu = muv;
  DT = 0.006;

  run();
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
init(i = 0)
{
  foreach ()
    u.x[] = cs[] ? U0 : 0.;
}

event
defaults(i = 0)
{
  int n_meshes = 1;
  ib_mesh_manager_init(n_meshes);

  // printf("rat: %f\n", b_ds / f_dx);

  ib_beam_t beam1 =
    ib_beam_new_theta(b_length, b_EI, b_mu, b_nodes, b_r, b_theta);
  ib_mesh_model(ibmesh(0), beam1);

  foreach_ibnode(ibmesh(0))
  {
    node->weight = b_ds;
    node->gravity.x = b_gravity;
  }
}

event
properties(i++)
{
  foreach_face() muv.x[] = fm.x[]*b_length*U0 / Reynolds;
}

// event
// progress_output(i+=10)
// {
//   for (int mi = 0; mi < ib_mesh_manager.nm; mi++) {
//     vertex_t c = ib_mesh_get_centroid(ibmesh(mi));
// #if dimension == 2
//     printf("[Mesh #%d, Time: %f]: (%f,%f)\n", mi, t, c.x, c.y);
// #else
//     printf("[Mesh #%d, Time: %f]: (%f,%f,%f)\n", mi, t, c.x, c.y, c.z);
// #endif
//   }
// }

event
tip_logging(i++)
{
  vertex_t tip = ibmesh(0)->nodes[0].lagpos;
  fprintf(stdout, "%g %g %g\n", t, tip.x, tip.y);
}

#include "library/ibm/IBoutput.h"

scalar omega[];

event
vtk(t += 0.20; t <= T_END)
{
  vorticity(u, omega);
  output_hdf_htg({ omega, p, stencils },
                 { u, forcing, tmp_vel, tmp_force },
                 "huang_serial_fluid");
  output_ibnodes("huang_serial", i, t);
}

event
adapt(i++)
{
  ib_mesh_manager_tag_stencil();
  adapt_wavelet({ cs, u, omega, stencils },
                (double[]){
                  1e-3,
                  3e-2,
                  3e-2,
                  3e-2,
                  3e-2,
                },
                maxlevel,
                4);
}

event
end(t = T_END)
{
  return 0;
}