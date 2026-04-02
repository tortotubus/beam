#pragma once

#include "elff/c/models/ibm/IBMesh.h"

#ifdef __cplusplus
extern "C"
{
#endif

  typedef void* ib_euler_beam_huang_t;

  ib_euler_beam_huang_t ib_euler_beam_huang_new(vertex_t s0,
                                                int bc_type_1,
                                                int bc_type_2,
                                                double length,
                                                double EI,
                                                double mu,
                                                int nodes);

  ib_euler_beam_huang_t ib_euler_beam_huang_new_theta(vertex_t s0,
                                                      int bc_type_1,
                                                      int bc_type_2,
                                                      double length,
                                                      double EI,
                                                      double mu,
                                                      int nodes,
                                                      double theta);

  void ib_euler_beam_huang_destroy(ib_euler_beam_huang_t handle);

#ifdef __cplusplus
}
#endif
