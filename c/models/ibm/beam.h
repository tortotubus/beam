#pragma once

#include "immersedboundary.h"

#ifdef __cplusplus
extern "C"
{
#endif

  typedef void* ib_beam_t;

  // Constructor: returns a handle to a new IBSpringCircle
  ib_beam_t ib_beam_new(double length,
                        double EI,
                        double mu,
                        int nodes,
                        double r_penalty);

  ib_beam_t ib_beam_new_theta(
                      double length,
                      double EI,
                      double mu,
                      int nodes,
                      double r_penalty,
                      double theta
                    );


  // Destructor: destroys the IBSpringCircle object
  void ib_beam_destroy(ib_beam_t handle);

#ifdef __cplusplus
}
#endif