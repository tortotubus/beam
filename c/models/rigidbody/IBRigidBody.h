#pragma once

#include "elff/c/models/ibm/IBMesh.h"

#ifdef __cplusplus
extern "C"
{
#endif

  typedef struct
  {
    double w, x, y, z;
  } quaternion_t;

  typedef void* ib_rigid_body_circle_t;
  typedef void* ib_rigid_body_fibre_t;
  typedef void* ib_rigid_body_sphere_t;
  typedef void* ib_pinned_rigid_body_circle_t;
  typedef void* ib_pinned_rigid_body_fibre_t;
  typedef void* ib_pinned_rigid_body_sphere_t;

  ib_rigid_body_circle_t ib_rigid_body_circle_new(double radius,
                                                  int point_count,
                                                  double density,
                                                  vertex_t center,
                                                  vertex_t velocity,
                                                  double angle,
                                                  double angular_velocity_z);

  void ib_rigid_body_circle_destroy(ib_rigid_body_circle_t handle);

  ib_rigid_body_fibre_t ib_rigid_body_fibre_new(
    double length,
    double diameter,
    int nodes,
    double linear_density,
    vertex_t center,
    vertex_t velocity,
    quaternion_t q_bw,
    vertex_t angular_velocity_body);

  void ib_rigid_body_fibre_destroy(ib_rigid_body_fibre_t handle);

  ib_pinned_rigid_body_circle_t ib_pinned_rigid_body_circle_new(
    double radius,
    int point_count,
    double density,
    vertex_t center,
    double angle);

  void ib_pinned_rigid_body_circle_destroy(
    ib_pinned_rigid_body_circle_t handle);

  ib_pinned_rigid_body_fibre_t ib_pinned_rigid_body_fibre_new(
    double length,
    double diameter,
    int nodes,
    double linear_density,
    vertex_t center,
    quaternion_t q_bw);

  void ib_pinned_rigid_body_fibre_destroy(
    ib_pinned_rigid_body_fibre_t handle);

  ib_rigid_body_sphere_t ib_rigid_body_sphere_new(
    double radius,
    int point_count,
    double density,
    vertex_t center,
    vertex_t velocity,
    quaternion_t q_bw,
    vertex_t angular_velocity_body);

  void ib_rigid_body_sphere_destroy(ib_rigid_body_sphere_t handle);

  ib_pinned_rigid_body_sphere_t ib_pinned_rigid_body_sphere_new(
    double radius,
    int point_count,
    double density,
    vertex_t center,
    quaternion_t q_bw);

  void ib_pinned_rigid_body_sphere_destroy(
    ib_pinned_rigid_body_sphere_t handle);

#ifdef __cplusplus
}
#endif
