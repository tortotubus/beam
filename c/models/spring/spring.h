#pragma once

#include "immersedboundary.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /**
   * @name ib_spring_circle_t
   * @brief C handle and related functions for spring-circle models.
   * @{
   */
  /**
   * @brief Opaque handle to an immersed-boundary spring-circle model.
   */
  typedef void* ib_spring_circle_t;

  /**
   * @brief Creates a circular spring network.
   *
   * @param npoints Number of Lagrangian points
   * @param K Spring stiffness
   * @param radius Circle radius
   * @param center_x Circle center x-coordinate
   * @param center_y Circle center y-coordinate
   * @return Opaque spring-circle handle
   */
  ib_spring_circle_t ib_spring_circle_create_circle(int npoints,
                                                    double K,
                                                    double radius,
                                                    double center_x,
                                                    double center_y);

  /**
   * @brief Creates an elliptical spring network.
   *
   * @param npoints Number of Lagrangian points
   * @param K Spring stiffness
   * @param radius_x Ellipse radius along x
   * @param radius_y Ellipse radius along y
   * @param center_x Ellipse center x-coordinate
   * @param center_y Ellipse center y-coordinate
   * @return Opaque spring-circle handle
   */
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

  /**
   * @brief Destroys a spring-circle model.
   *
   * @param handle Opaque spring-circle handle
   */
  void ib_spring_circle_destroy(ib_spring_circle_t handle);
  /** @} */

#ifdef __cplusplus
}
#endif
