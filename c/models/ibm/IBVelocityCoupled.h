#pragma once

#include "elff/c/models/ibm/IBMesh.h"


#ifdef __cplusplus
extern "C"
{
#endif

  /**
   * @name ib_velocity_coupled_t
   * @brief C handle and related functions for velocity-coupled IB models.
   * @{
   */
  /** 
   * @brief Opaque handle to an @ref ELFF::Models::IBVelocityCoupled instance.
   */
  typedef void* ib_velocity_coupled_t;

  /**
   * @brief Get the number of nodes of the model
   *
   * @param handle Opaque model handle
   * @return Number of nodes
   */
  int ib_velocity_coupled_get_number_of_nodes(
    ib_velocity_coupled_t handle);

  /**
   * @brief Get the "current" timestep mesh. 
   *
   * @param handle Opaque model handle
   * @return Current mesh state
   */
  ib_mesh_t ib_velocity_coupled_get_current(
    ib_velocity_coupled_t handle);

  /**
   * @brief Advances the model to the midpoint state.
   *
   * @param handle Opaque model handle
   * @param velocity Array of nodal velocities
   * @param n Number of velocity entries
   * @param dt Time-step size
   * @return Midpoint mesh state
   */
  ib_mesh_t ib_velocity_coupled_get_midpoint(
    ib_velocity_coupled_t handle,
    vertex_t* velocity,
    int n,
    double dt);

  /**
   * @brief Advances the model to the next state.
   *
   * @param handle Opaque model handle
   * @param velocity Array of nodal velocities
   * @param n Number of velocity entries
   * @param dt Time-step size
   * @return Next mesh state
   */
  ib_mesh_t ib_velocity_coupled_get_next(
    ib_velocity_coupled_t handle,
    vertex_t* velocity,
    int n,
    double dt);
  /** @} */

#ifdef __cplusplus
}
#endif
