#define DEBUG 0

#define LEVEL 9
#define DT_SAVE 0.1
#define T_END 11.6

scalar f[];
scalar* tracers = { f };
double Reynolds = 200.;
int maxlevel = 9;
int Nf = 1 << 9;
face vector muv[];

double U0 = 1.;

#include "grid/quadtree.h"

#include "interface/ibm/beam.h"
#include "library/ibm/uzawa2.h"
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
defaults(i = 0)
{
  int n_meshes = 1;
  ib_mesh_manager_init(n_meshes);

  double f_dx = L0 / Nf;
  // double b_dx = f_dx * 1.683;

  double b_length = 1;
  double b_EI = 0.0001;
  double b_mu = 1.5;
  int b_nodes = 65;
  double b_ds = b_length / b_nodes;
  double b_r = 1e4;
  double b_theta = 0.1 * pi;

  printf("rat: %f\n", b_ds / f_dx);

  ib_beam_t beam1 =
    ib_beam_new_theta(b_length, b_EI, b_mu, b_nodes, b_r, b_theta);
  ib_mesh_model(ibmesh(0), beam1);

  foreach_ibnode(ibmesh(0))
  {
    node->weight = b_mu * b_ds;
    node->gravity.x = 1;
  }
}

event
properties(i++)
{
  foreach_face() muv.x[] = 1. / Reynolds;
}

event
progress_output(i += 20)
{
  for (int mi = 0; mi < ib_mesh_manager.nm; mi++) {
    vertex_t c = ib_mesh_get_centroid(ibmesh(mi));
#if dimension == 2
    printf("[Mesh #%d, Time: %f]: (%f,%f)\n", mi, t, c.x, c.y);
#else
    printf("[Mesh #%d, Time: %f]: (%f,%f,%f)\n", mi, t, c.x, c.y, c.z);
#endif
  }
}

event
vtk(t += DT_SAVE; t <= T_END)
{
  scalar omega[];
  vorticity(u, omega);
  output_hdf_htg({ omega, p, stencils, f },
                 { u, forcing, tmp_vel, tmp_force },
                 "huang_serial_fluid");
}

// Write all IB node positions to a legacy VTK PolyData file
event
output_ib_points(t += DT_SAVE;
                 t <= T_END) // change 10 to whatever stride you want
{
  if (pid() == 0) {
    // Count total active points
    int total_points = 0;
    foreach_ibmesh()
    {
      if (mesh->isactive) {
        total_points += mesh->nn;
      }
    }

    // if (total_points == 0)
    // return;

    char fname[128];
    sprintf(fname, "huang_serial_nodes_%d.vtk", i);
    FILE* fp = fopen(fname, "w");

    // if (!fp) {
    //   perror("fopen ib_points vtk");
    //   return;
    // }

    // Legacy VTK header
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "IB points at step %d\n", i);
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET POLYDATA\n");
    fprintf(fp, "POINTS %d double\n", total_points);

    // Write point coordinates
    for (int mi = 0; mi < ib_mesh_manager.nm; mi++) {
      IBMesh* mesh = ibmesh(mi);
      if (!mesh->isactive)
        continue;

      for (int ni = 0; ni < mesh->nn; ni++) {
        IBNode* node = &mesh->nodes[ni];

#if dimension == 1
        double x = node->lagpos.x;
        fprintf(fp, "%.16g %.16g %.16g\n", x, 0.0, 0.0);
#elif dimension == 2
        double x = node->lagpos.x;
        double y = node->lagpos.y;
        fprintf(fp, "%.16g %.16g %.16g\n", x, y, 0.0);
#else // dimension == 3
        double x = node->lagpos.x;
        double y = node->lagpos.y;
        double z = node->lagpos.z;
        fprintf(fp, "%.16g %.16g %.16g\n", x, y, z);
#endif
      }
    }

    // Optional: make them actual VERTICES so ParaView shows them naturally
    fprintf(fp, "VERTICES %d %d\n", total_points, 2 * total_points);
    int idx = 0;
    for (int mi = 0; mi < ib_mesh_manager.nm; mi++) {
      IBMesh* mesh = ibmesh(mi);
      if (!mesh->isactive)
        continue;
      for (int ni = 0; ni < mesh->nn; ni++) {
        fprintf(fp, "1 %d\n", idx++);
      }
    }

    fclose(fp);
  }
}

event
adapt(i++)
{
  ib_mesh_manager_tag_stencil();
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