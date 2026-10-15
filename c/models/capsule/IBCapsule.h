#pragma once

#include "elff/c/models/ibm/IBMesh.h"
#include "elff/c/models/ibm/IBVelocityCoupled.h"

#ifdef __cplusplus
extern "C"
{
#endif

  typedef void* ib_capsule_t;

  ib_capsule_t ib_capsule_sphere_new(double radius,
                                     vertex_t center,
                                     int refinements,
                                     double elastic_modulus);

  ib_capsule_t ib_capsule_ellipsoid_new(vertex_t radii,
                                        vertex_t center,
                                        int refinements,
                                        double elastic_modulus);

  ib_capsule_t ib_capsule_biconcave_new(double radius,
                                        vertex_t center,
                                        int refinements,
                                        double elastic_modulus);

  void ib_capsule_set_neo_hookean_law(ib_capsule_t handle,
                                      double elastic_modulus);

  void ib_capsule_set_skalak_law(ib_capsule_t handle,
                                 double elastic_modulus,
                                 double area_dilatation_modulus);

  void ib_capsule_set_linear_bending_law(ib_capsule_t handle,
                                         double bending_modulus);

  void ib_capsule_set_helfrich_bending_law(ib_capsule_t handle,
                                           double bending_modulus);

  void ib_capsule_set_constant_reference_curvature(
    ib_capsule_t handle,
    double reference_curvature);

  void ib_capsule_clear_bending_law(ib_capsule_t handle);

  int ib_capsule_get_triangle_count(ib_capsule_t handle);

  void ib_capsule_get_triangle_node_ids(ib_capsule_t handle,
                                        int triangle_index,
                                        int node_ids[3]);

  ib_velocity_coupled_t ib_capsule_as_velocity_coupled(ib_capsule_t handle);

  void ib_capsule_destroy(ib_capsule_t handle);

#ifdef __cplusplus
}
#endif
