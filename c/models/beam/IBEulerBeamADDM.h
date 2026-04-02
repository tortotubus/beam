#pragma once

#include "elff/c/models/ibm/IBMesh.h"

#ifdef __cplusplus
extern "C"
{
#endif

  typedef void* ib_euler_beam_addm_t;

  ib_euler_beam_addm_t ib_euler_beam_addm_new(vertex_t s0,
                                              int bc_type_1,
                                              int bc_type_2,
                                              double length,
                                              double EI,
                                              double mu,
                                              int nodes,
                                              double r_penalty);

  ib_euler_beam_addm_t ib_euler_beam_addm_new_theta(vertex_t s0,
                                                    int bc_type_1,
                                                    int bc_type_2,
                                                    double length,
                                                    double EI,
                                                    double mu,
                                                    int nodes,
                                                    double r_penalty,
                                                    double theta);

  void ib_euler_beam_addm_destroy(ib_euler_beam_addm_t handle);

#ifdef __cplusplus
}
#endif
