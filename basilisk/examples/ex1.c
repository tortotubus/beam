#define DEBUG_PRINT_CALL_OR_EVENT 1

#define LEVEL 5
#define LENGTH 10.
#define T_END 10

#include "grid/quadtree.h"

#include "library/ibm.h"
#include "interface/spring.h"

#include "library/output.h"

int
main()
{
  L0 = LENGTH;
  origin(-.5 * L0, -.5 * L0);
  init_grid(1 << LEVEL);
  periodic(left);
  DT = 5.e-3;

  run();
}

event defaults(i = 0) {
  int n_meshes = 1;
  ib_mesh_manager_init(n_meshes);
  
  int n_points = 30;
  double K = 0.0;
  double radius_x = 0.5, radius_y = 0.6;
  double center_x = -.1, center_y = 0.1;
  ib_spring_circle_t spr_c = ib_spring_circle_create_ellipse(n_points,K,radius_x, radius_y, center_x, center_y);
  ib_mesh_set_model(&IBMESH(0), spr_c);

  printf("Dimension: %d\n", dimension);

  for (int mi = 0; mi < ib_mesh_manager->nm; mi++) {
    printf("[Mesh #%d]: ", mi);
    ib_mesh_print_pos(&IBMESH(mi));
  }
}


event
impose_u(i++)
{
  foreach () {
    u.x[] = 1.;
    u.y[] = 0.;
  }
}

event
progress_output(i++)
{
  for (int mi = 0; mi < ib_mesh_manager->nm; mi++) {
    vertex_t c = ib_mesh_get_centroid(&IBMESH(mi));
    #if dimension == 2
    printf("[Mesh #%d, Time: %f]: (%f,%f)\n", mi, t, c.x, c.y);
    #else 
    printf("[Mesh #%d, Time: %f]: (%f,%f,%f)\n", mi, t, c.x, c.y, c.z);
    #endif
  }
}

event adapt (i++) {
  tag_stencils();
  adapt_wavelet({stencils}, (double []){1.e-2}, maxlevel = LEVEL);
  generate_stencils();
}

//event
//movies(i += 5; t <= T_END)
//{
//  scalar omega[];
//  vorticity(u, omega);
//  output_ppm(omega, file = "vort.mp4", linear = true);
//}

event
vtk(i += 5; t <= T_END)
{
  scalar omega[];
  vorticity(u, omega);
  output_hdf_htg({omega,p},{forcing},"ex1");
  output_hdf_polydata("ex1");
}

event
end(t = T_END)
{
  return 0;
}