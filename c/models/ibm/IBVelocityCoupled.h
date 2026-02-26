#pragma once

#include "c/models/ibm/IBMesh.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /**
   * @class ib_velocity_coupled_t
   *
   * @brief Pointer handle to any object of the @ref
   * ELFF::IBIBVelocityCoupledStructureModel.
   */
  typedef void* ib_velocity_coupled_t;

  /**
   * @brief Get the number of nodes of the model
   *
   * @memberof ib_velocity_coupled_t
   */
  int ib_velocity_coupled_get_number_of_nodes(
    ib_velocity_coupled_t handle);

  /**
   * @brief Get the "current" timestep mesh. @sa ELFF::IBVelocityStructureModel
   *
   * @memberof ib_velocity_coupled_t
   */
  ib_mesh_t ib_velocity_coupled_get_current(
    ib_velocity_coupled_t handle);

  /**
   * @memberof ib_velocity_coupled_t
   */
  ib_mesh_t ib_velocity_coupled_get_midpoint(
    ib_velocity_coupled_t handle,
    vertex_t* velocity,
    int n,
    double dt);

  /**
   * @memberof ib_velocity_coupled_t
   */
  ib_mesh_t ib_velocity_coupled_get_next(
    ib_velocity_coupled_t handle,
    vertex_t* velocity,
    int n,
    double dt);

#ifdef __cplusplus
}
#endif