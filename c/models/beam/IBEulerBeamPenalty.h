#pragma once

#include "elff/c/models/beam/IBEulerBeamBCs.h"
#include "elff/c/models/ibm/IBRuntime.h"

#ifdef __cplusplus
extern "C"
{
#endif

  typedef void* ib_euler_beam_penalty_t;

  ib_euler_beam_penalty_t ib_euler_beam_penalty_new(ib_runtime_t runtime,
                                                    vertex_t s0,
                                                    ib_euler_beam_bcs_t bcs,
                                                    double length,
                                                    double EI,
                                                    double mu,
                                                    int nodes,
                                                    double r_penalty);

  ib_euler_beam_penalty_t ib_euler_beam_penalty_new_theta(ib_runtime_t runtime,
                                                          vertex_t s0,
                                                          ib_euler_beam_bcs_t bcs,
                                                          double length,
                                                          double EI,
                                                          double mu,
                                                          int nodes,
                                                          double r_penalty,
                                                          double theta);

#ifdef __cplusplus
}
#endif
