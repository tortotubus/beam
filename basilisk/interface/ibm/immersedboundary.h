#pragma once

#ifdef __cplusplus
extern "C"
{
#endif

  /**
   * @brief Helper class for storing coordinates
   */
  typedef struct
  {
    double x, y, z;
  } vertex_t;

  /**
   * @brief Struct for representing the @ref ELFF::IBStructureMesh
   */
  typedef struct
  {
    int n;
    vertex_t* position;
    vertex_t* velocity;
    vertex_t* forces;
  } ib_structure_mesh_t;

  /**
   * @memberof ib_structure_mesh_t
   */
  void ib_structure_mesh_free(ib_structure_mesh_t* mesh);

  /**
   * @class ib_velocity_structure_model_t
   *
   * @brief Pointer handle to any object of the @ref
   * ELFF::IBVelocityCoupledStructureModel.
   */
  typedef void* ib_velocity_structure_model_t;

  /**
   * @brief Get the number of nodes of the model
   *
   * @memberof ib_velocity_structure_model_t
   */
  int ib_velocity_structure_model_get_number_of_nodes(
    ib_velocity_structure_model_t handle);

  /**
   * @brief Get the "current" timestep mesh. @sa ELFF::IBVelocityStructureModel
   *
   * @memberof ib_velocity_structure_model_t
   */
  ib_structure_mesh_t ib_velocity_structure_model_get_current(
    ib_velocity_structure_model_t handle);

  /**
   * @memberof ib_velocity_structure_model_t
   */
  ib_structure_mesh_t ib_velocity_structure_model_get_midpoint(
    ib_velocity_structure_model_t handle,
    vertex_t* velocity,
    int n,
    double dt);

  /**
   * @memberof ib_velocity_structure_model_t
   */
  ib_structure_mesh_t ib_velocity_structure_model_get_next(
    ib_velocity_structure_model_t handle,
    vertex_t* velocity,
    int n,
    double dt);

  /**
   * @class ib_force_structure_model_t
   *
   * @brief Pointer handle to any object of the @ref
   * ELFF::IBForceCoupledStructureModel.
   */
  typedef void* ib_force_structure_model_t;

  /**
   * @memberof ib_force_structure_model_t
   */
  int ib_force_structure_model_get_number_of_nodes(
    ib_force_structure_model_t handle);

  /**
   * @memberof ib_force_structure_model_t
   */
  ib_structure_mesh_t ib_force_structure_model_get_current(
    ib_force_structure_model_t handle);

  /**
   * @memberof ib_force_structure_model_t
   */
  ib_structure_mesh_t ib_force_structure_model_get_next(
    ib_force_structure_model_t handle,
    vertex_t* velocity,
    int n,
    double dt);

#ifdef __cplusplus
}
#endif