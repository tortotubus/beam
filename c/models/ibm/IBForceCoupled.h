#pragma once

#include "elff/c/models/ibm/IBMesh.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /**
   * @name ib_force_coupled_t
   * @brief C handle and related functions for force-coupled IB models.
   * @{
   */
  /**
   * @brief Opaque handle to an @ref ELFF::Models::IBForceCoupled instance.
   */
  typedef void* ib_force_coupled_t;

  /**
   * @brief Returns the number of nodes in the coupled model.
   *
   * @param handle Opaque model handle
   * @return Number of nodes
   */
  int ib_force_coupled_get_number_of_nodes(
    ib_force_coupled_t handle);

  /**
   * @brief Returns the current immersed-boundary mesh.
   *
   * @param handle Opaque model handle
   * @return Current mesh state
   */
  ib_mesh_t ib_force_coupled_get_current(
    ib_force_coupled_t handle);

  /**
   * @brief Advances the model and returns the next mesh state.
   *
   * @param handle Opaque model handle
   * @param force Array of nodal forces
   * @param n Number of force entries
   * @param dt Time-step size
   * @return Next mesh state
   */
  ib_mesh_t ib_force_coupled_get_next(
    ib_force_coupled_t handle,
    vertex_t* force,
    int n,
    double dt);
  /** @} */

#ifdef __cplusplus
}
#endif
