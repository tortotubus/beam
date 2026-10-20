#pragma once

#include "elff/c/models/beam/IBEulerBeamBCs.h"
#include "elff/c/models/ibm/IBCoupling.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /**
   * @name ib_euler_beam_t
   * @brief C handle and related functions for immersed-boundary beam models.
   * @{
   */
  /**
   * @brief Opaque handle to a beam model exposed through the C API.
   */
  typedef void* ib_euler_beam_t;

  /**
   * @brief Creates a new beam model with the given initial slope.
   *
   * @param s0 Initial slope vector
   * @param bcs Boundary conditions at both beam ends
   * @param length Beam length
   * @param EI Flexural rigidity
   * @param mu Mass per unit length
   * @param nodes Number of discretization nodes
   * @param r_penalty Penalty parameter used in the inextensibility constraint
   * @return Opaque beam handle
   */
  ib_euler_beam_t ib_euler_beam_new(ib_coupling_t runtime,
                        vertex_t s0,
                        ib_euler_beam_bcs_t bcs,
                        double length,
                        double EI,
                        double mu,
                        int nodes,
                        double r_penalty);

  /**
   * @brief Creates a new beam model with an additional initial angle.
   *
   * @param s0 Initial slope vector
   * @param bcs Boundary conditions at both beam ends
   * @param length Beam length
   * @param EI Flexural rigidity
   * @param mu Mass per unit length
   * @param nodes Number of discretization nodes
   * @param r_penalty Penalty parameter used in the inextensibility constraint
   * @param theta Initial angle parameter
   * @return Opaque beam handle
   */
  ib_euler_beam_t ib_euler_beam_new_theta(ib_coupling_t runtime,
                              vertex_t s0,
                              ib_euler_beam_bcs_t bcs,
                              double length,
                              double EI,
                              double mu,
                              int nodes,
                              double r_penalty,
                              double theta);

  /** @} */

#ifdef __cplusplus
}
#endif
