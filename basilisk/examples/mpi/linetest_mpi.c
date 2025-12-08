

#include "grid/quadtree.h"

#include "interface/ibm/beam.h"
#include "library/ibm/uzawa3.h"

event
defaults(i = 0)
{
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face() a.x[] = 0.;
  }
}

event
acceleration(i++)
{
  face vector ae = a;

  // 0. Clear the forcing term on the Eulerian grid
  foreach () {
    if (cm[] > 1e-20) {
      foreach_dimension()
      {
        forcing.x[] = 0.;
      }
    }
  }

  // 1. Advance lagrangian mesh
  // ib_mesh_manager_advance_lagrangian_mesh ();

  // 2. Generate the stencil cache
  ib_mesh_manager_generate_stencil_cache();

  // 3. Interpolate \f(u*\f) with \f(J u*\f)
  ib_mesh_manager_interpolate_eulerian_velocities();

  // 4. Build the right-hand side vector \f(b = v_{\rm lag} - J u^* \f)
  // ib_mesh_manager_compute_constraint_rhs();

  // 5. Solve A \lambda = b via CG/Uzawa on nodes
  // ib_mesh_manager_solve_lambda_CG(dt);

  foreach_ibmesh() {
    foreach_ibnode(mesh) {
      node->force.x = 1.;
      node->force.y = 1.;
    }
  }

  // 6. Spread lambda to grid forcing \f(J^T \lambda\f)
  ib_mesh_manager_spread_eulerian_forcing(forcing);

  // 7. Add forcing to the acceleration field
  foreach_face()
  {
    if (fm.x[] > 1e-20) {
      ae.x[] += .5 * alpha.x[] * (forcing.x[] + forcing.x[-1]);
    }
  }
}

event
acceleration(i++)
{
  // ib_mesh_manager_sort_triplets ();
}

/** At the end of the simulation, we free the allocated memory.*/
event
cleanup(t = end)
{
  ib_mesh_manager_free();
}
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
    // node->lagpos.x = 0.;
    node->lagpos.y = -.5 + node_index * b_ds;
    node->weight = b_ds;
    node->gravity.x = 1;
  }
}
event
defaults(i = 0)
{
  int n_meshes = 1;
  ib_mesh_manager_init(n_meshes);

  ib_mesh_init_nn(ibmesh(0), 65);

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
    for (int mi = 0; mi < ib_mesh_manager.nm; mi++) {
      vertex_t c = ib_mesh_get_centroid(ibmesh(mi));
#if dimension == 2
      printf("[Mesh #%d, Time: %f]: (%f,%f)\n", mi, t, c.x, c.y);
#else
      printf("[Mesh #%d, Time: %f]: (%f,%f,%f)\n", mi, t, c.x, c.y, c.z);
#endif
    }
  }
}

event
vtk(t += DT_SAVE; t <= T_END)
{
  scalar omega[];
  vorticity(u, omega);
  output_hdf_htg({ omega, p, stencils, f },
                 { u, forcing, tmp_vel, tmp_force },
                 "linetest_mpi_fluid");
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
    sprintf(fname, "linetest_mpi_nodes_%d.vtk", i);
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