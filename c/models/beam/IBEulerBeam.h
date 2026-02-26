#pragma once

#ifdef __cplusplus
extern "C"
{
#endif

  typedef void* ib_beam_t;

  ib_beam_t ib_beam_new(double length,
                        double EI,
                        double mu,
                        int nodes,
                        double r_penalty);

  ib_beam_t ib_beam_new_theta(double length,
                              double EI,
                              double mu,
                              int nodes,
                              double r_penalty,
                              double theta);

  void ib_beam_destroy(ib_beam_t handle);

#ifdef __cplusplus
}
#endif