#define DEBUG_PRINT_CALL_OR_EVENT 0

#define LEVEL 8
#define LENGTH 10.
#define DT_SAVE 0.05
#define T_END 0.5

#include "grid/quadtree.h"

#include "library/ibm/immersedboundary.h"
#include "interface/ibm/spring.h"
#include "library/output_htg.h"

int
main()
{
  L0 = LENGTH;
  origin(-.5 * L0, -.5 * L0);
  init_grid(1 << LEVEL);
  periodic(left);
  //stokes = true;
  DT = 5.e-3;
  run();
}

event
defaults(i = 0)
{
  int n_meshes = 2;
  ib_mesh_manager_init(n_meshes);

  int n_points;
  double K;
  double radius_x, radius_y;
  double center_x, center_y;

  n_points = 100;
  K = .1;
  radius_x = 0.5, radius_y = 0.8;
  center_x = -.1, center_y = 0.1;
  ib_spring_circle_t spr_c0 = ib_spring_circle_create_ellipse(
    n_points, K, radius_x, radius_y, center_x, center_y);
  ib_mesh_set_velocity_model(&IBMESH(0), spr_c0);

  n_points = 100;
  K = .1;
  radius_x = 0.6, radius_y = 0.4;
  center_x = 2., center_y = 1.;
  ib_spring_circle_t spr_c1 = ib_spring_circle_create_ellipse(
    n_points, K, radius_x, radius_y, center_x, center_y);
  ib_mesh_set_velocity_model(&IBMESH(1), spr_c1);
}

event
progress_output(i++)
{
  for (int mi = 0; mi < ib_mesh_manager.nm; mi++) {
    vertex_t c = ib_mesh_get_centroid(&IBMESH(mi));
#if dimension == 2
    printf("[Mesh #%d, Time: %f]: (%f,%f)\n", mi, t, c.x, c.y);
#else
    printf("[Mesh #%d, Time: %f]: (%f,%f,%f)\n", mi, t, c.x, c.y, c.z);
#endif
  }
}

event
adapt(i++)
{
  tag_stencils();
  adapt_wavelet({ stencils }, (double[]){ 1.e-2 }, maxlevel = LEVEL);
  generate_stencil_caches();
}

// event
// movies(i += 5; t <= T_END)
//{
//   scalar omega[];
//   vorticity(u, omega);
//   output_ppm(omega, file = "vort.mp4", linear = true);
// }

event
vtk(t += DT_SAVE; t <= T_END)
{
  scalar omega[];
  vorticity(u, omega);
  output_hdf_htg({ omega, p, stencils }, { u }, "ex1");
  // output_hdf_polydata("ex1");
}


// Write all IB node positions to a legacy VTK PolyData file
event output_ib_points (t += DT_SAVE; t <= T_END) // change 10 to whatever stride you want
{
  // Count total active points
  int total_points = 0;
  for (int mi = 0; mi < ib_mesh_manager.nm; mi++) {
    IBMesh *mesh = &IBMESH(mi);
    if (!mesh->isactive)
      continue;
    total_points += mesh->nn;
  }

  // if (total_points == 0)
    // return;

  char fname[128];
  sprintf(fname, "ib_points-%05d.vtk", i);
  FILE *fp = fopen(fname, "w");

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
    IBMesh *mesh = &IBMESH(mi);
    if (!mesh->isactive)
      continue;

    for (int ni = 0; ni < mesh->nn; ni++) {
      IBNode *node = &mesh->nodes[ni];

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
  fprintf(fp, "VERTICES %d %d\n", total_points, 2*total_points);
  int idx = 0;
  for (int mi = 0; mi < ib_mesh_manager.nm; mi++) {
    IBMesh *mesh = &IBMESH(mi);
    if (!mesh->isactive)
      continue;
    for (int ni = 0; ni < mesh->nn; ni++) {
      fprintf(fp, "1 %d\n", idx++);
    }
  }

  fclose(fp);
}

event
end(t = T_END)
{
  return 0;
}