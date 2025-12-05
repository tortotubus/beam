#pragma once

#include "immersedboundary.h"

#ifdef __cplusplus
extern "C"
{
#endif

  typedef void* ib_spring_circle_t;

  // Constructor: returns a handle to a new IBSpringCircle
  ib_spring_circle_t ib_spring_circle_create_circle(int npoints,
                                                    double K,
                                                    double radius,
                                                    double center_x,
                                                    double center_y);

  ib_spring_circle_t ib_spring_circle_create_ellipse(int npoints,
                                                     double K,
                                                     double radius_x,
                                                     double radius_y,
                                                     double center_x,
                                                     double center_y);
  // ib_structure_mesh_t ib_spring_circle_get_current(ib_spring_circle_t
  // handle); ib_structure_mesh_t
  // ib_structure_model_get_midpoint(ib_spring_circle_t handle, vertex_t*
  // velocity, int n, double dt); ib_structure_mesh_t
  // ib_structure_model_get_next(ib_spring_circle_t handle, vertex_t* velocity,
  // int n, double dt);

  // Destructor: destroys the IBSpringCircle object
  void ib_spring_circle_destroy(ib_spring_circle_t handle);

#ifdef __cplusplus
}
#endif