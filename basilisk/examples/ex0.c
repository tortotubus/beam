#include <stdlib.h>
#include <stdio.h>

#include "interface/dummy.h"
#include "interface/immersedboundary.h"
#include "interface/spring.h"

int
main()
{
  /*
    Test dummy.h
  */
  int base = 1;
  int sq = 2;
  int result = 0;

  beam_handle_t* dummy = beam_dummy_create(base);
  result = beam_dummy_compute(dummy, sq);
  beam_dummy_destroy(dummy);

  printf("Dummy: %d\n", result);

  /*
    Test immersedboundary.h and spring.h
  */

  int n_nodes = 10;
  int n_dims = 3;

  double K = 1.;
  double radius_x = 1., radius_y = .8;
  double center_x = 0., center_y = 0.;

  ib_spring_circle_t spring_circle =  ib_spring_circle_create_ellipse(n_nodes, K, radius_x, radius_y, center_x, center_y);
  vertex_t *velocity = calloc(n_nodes, sizeof(vertex_t));

  //ib_structure_mesh_t mesh = ib_structure_model_get_current(spring_circle);
  ib_structure_mesh_t mesh = ib_structure_model_get_midpoint(spring_circle, velocity, n_nodes, 0.);

  printf("Spring Circle Points: {");
  for (int i = 0; i < n_nodes; i++) {
    printf("(%f,%f,%f), ", mesh.points[i].x, mesh.points[i].y, mesh.points[i].z);
  }
  printf("}\n");

  printf("Spring Circle Forces: {");
  for (int i = 0; i < n_nodes; i++) {
    printf("(%f,%f,%f), ", mesh.forces[i].x, mesh.forces[i].y, mesh.forces[i].z);
  }
  printf("}\n");
  
  //ib_structure_mesh_free(&mesh);
  
  return 0;
}