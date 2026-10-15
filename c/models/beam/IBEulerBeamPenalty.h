#pragma once

#include "elff/c/models/beam/IBEulerBeamBCs.h"

#ifdef __cplusplus
extern "C"
{
#endif

  typedef void* ib_euler_beam_penalty_t;

  ib_euler_beam_penalty_t ib_euler_beam_penalty_new(vertex_t s0,
                                                    ib_euler_beam_bcs_t bcs,
                                                    double length,
                                                    double EI,
                                                    double mu,
                                                    int nodes,
                                                    double r_penalty);

  ib_euler_beam_penalty_t ib_euler_beam_penalty_new_theta(vertex_t s0,
                                                          ib_euler_beam_bcs_t bcs,
                                                          double length,
                                                          double EI,
                                                          double mu,
                                                          int nodes,
                                                          double r_penalty,
                                                          double theta);

  void ib_euler_beam_penalty_destroy(ib_euler_beam_penalty_t handle);

#ifdef __cplusplus
}
#endif
