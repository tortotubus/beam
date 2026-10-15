#pragma once

#include "elff/c/models/ibm/IBMesh.h"

#ifdef __cplusplus
extern "C"
{
#endif

  typedef enum {
    IB_EULER_BEAM_BC_FREE = 0,
    IB_EULER_BEAM_BC_SIMPLE = 1,
    IB_EULER_BEAM_BC_CLAMPED = 2,
    IB_EULER_BEAM_BC_POINT_FORCE = 3,
    IB_EULER_BEAM_BC_POINT_TORQUE = 4,
  } ib_euler_beam_bc_type_t;

  typedef enum {
    IB_EULER_BEAM_BC_LEFT = 0,
    IB_EULER_BEAM_BC_RIGHT = 1,
  } ib_euler_beam_bc_end_t;

  typedef struct {
    vertex_t position;
    vertex_t slope;
    vertex_t force;
    vertex_t torque;
  } ib_euler_beam_bc_vals_t;

  typedef struct {
    int end[2];
    int type[2];
    ib_euler_beam_bc_vals_t vals[2];
  } ib_euler_beam_bcs_t;

#ifdef __cplusplus
}
#endif
