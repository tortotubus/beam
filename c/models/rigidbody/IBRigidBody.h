#pragma once

#include "elff/c/models/ibm/IBMesh.h"
#include "elff/c/models/ibm/IBCoupling.h"

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

  ib_rigid_body_circle_t ib_rigid_body_circle_new(ib_coupling_t runtime,
                                                  double radius,
                                                  int point_count,
                                                  double density,
                                                  vertex_t center,
                                                  vertex_t velocity,
                                                  double angle,
                                                  double angular_velocity_z);

  ib_rigid_body_fibre_t ib_rigid_body_fibre_new(
    ib_coupling_t runtime,
    double length,
    double diameter,
    int nodes,
    double linear_density,
    vertex_t center,
    vertex_t velocity,
    quaternion_t q_bw,
    vertex_t angular_velocity_body);

  ib_pinned_rigid_body_circle_t ib_pinned_rigid_body_circle_new(
    ib_coupling_t runtime,
    double radius,
    int point_count,
    double density,
    vertex_t center,
    double angle);

  ib_pinned_rigid_body_fibre_t ib_pinned_rigid_body_fibre_new(
    ib_coupling_t runtime,
    double length,
    double diameter,
    int nodes,
    double linear_density,
    vertex_t center,
    quaternion_t q_bw);

  ib_rigid_body_sphere_t ib_rigid_body_sphere_new(
    ib_coupling_t runtime,
    double radius,
    int point_count,
    double density,
    vertex_t center,
    vertex_t velocity,
    quaternion_t q_bw,
    vertex_t angular_velocity_body);

  ib_pinned_rigid_body_sphere_t ib_pinned_rigid_body_sphere_new(
    ib_coupling_t runtime,
    double radius,
    int point_count,
    double density,
    vertex_t center,
    quaternion_t q_bw);

#ifdef __cplusplus
}
#endif
