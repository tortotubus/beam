#pragma once

#include "c/models/ibm/IBMesh.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /**
   * @class ib_force_coupled_t
   *
   * @brief Pointer handle to any object of the @ref
   * ELFF::Model::IBForceCoupled.
   */
  typedef void* ib_force_coupled_t;

  /**
   * @memberof ib_force_coupled_t
   */
  int ib_force_structure_model_get_number_of_nodes(
    ib_force_coupled_t handle);

  /**
   * @memberof ib_force_coupled_t
   */
  ib_mesh_t ib_force_structure_model_get_current(
    ib_force_coupled_t handle);

  /**
   * @memberof ib_force_coupled_t
   */
  ib_mesh_t ib_force_structure_model_get_next(
    ib_force_coupled_t handle,
    vertex_t* velocity,
    int n,
    double dt);

#ifdef __cplusplus
}
#endif