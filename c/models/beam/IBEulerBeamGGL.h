#pragma once

#include "elff/c/models/beam/IBEulerBeamBCs.h"
#include "elff/c/models/ibm/IBCoupling.h"

#ifdef __cplusplus
extern "C"
{
#endif

  typedef void* ib_euler_beam_ggl_t;

  ib_euler_beam_ggl_t ib_euler_beam_ggl_new(ib_coupling_t runtime,
                                            vertex_t s0,
                                            ib_euler_beam_bcs_t bcs,
                                            double length,
                                            double EI,
                                            double mu,
                                            int nodes,
                                            double r_penalty);

  ib_euler_beam_ggl_t ib_euler_beam_ggl_new_theta(ib_coupling_t runtime,
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
