#pragma once

#include "elff/c/models/beam/IBEulerBeamBCs.h"

#ifdef __cplusplus
extern "C"
{
#endif

  typedef void* ib_euler_beam_huang_t;

  ib_euler_beam_huang_t ib_euler_beam_huang_new(vertex_t s0,
                                                ib_euler_beam_bcs_t bcs,
                                                double length,
                                                double EI,
                                                double mu,
                                                int nodes);

  ib_euler_beam_huang_t ib_euler_beam_huang_new_theta(vertex_t s0,
                                                      ib_euler_beam_bcs_t bcs,
                                                      double length,
                                                      double EI,
                                                      double mu,
                                                      int nodes,
                                                      double theta);

  void ib_euler_beam_huang_destroy(ib_euler_beam_huang_t handle);

#ifdef __cplusplus
}
#endif
