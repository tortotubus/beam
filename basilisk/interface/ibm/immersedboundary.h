#pragma once

#ifdef __cplusplus
extern "C"
{
#endif

  /**
   * Wrappers for beam::IBStructureMesh
   */

  typedef struct {
    double x, y, z;
  } vertex_t;

  typedef struct {
    int n;
    vertex_t *position;
    vertex_t *velocity;
    vertex_t *forces;
  } ib_structure_mesh_t;

  void ib_structure_mesh_free(ib_structure_mesh_t *mesh);

  /*
   * Wrappers for beam::IBVelocityCoupledStructureModel
   */

  typedef void* ib_velocity_structure_model_t;

  // int ib_velocity_structure_model_get_dimension       (ib_velocity_structure_model_t handle);
  int ib_velocity_structure_model_get_number_of_nodes (ib_velocity_structure_model_t handle);
  
  ib_structure_mesh_t ib_velocity_structure_model_get_current  (ib_velocity_structure_model_t handle);
  ib_structure_mesh_t ib_velocity_structure_model_get_midpoint (ib_velocity_structure_model_t handle, vertex_t* velocity, int n, double dt);
  ib_structure_mesh_t ib_velocity_structure_model_get_next     (ib_velocity_structure_model_t handle, vertex_t* velocity, int n, double dt);

  /**
   * Wrappers for beam::IBForceCoupledStructureModel
   */
  typedef void* ib_force_structure_model_t;

  int ib_force_structure_model_get_number_of_nodes(ib_force_structure_model_t handle);
  ib_structure_mesh_t ib_force_structure_model_get_current(ib_force_structure_model_t handle);
  ib_structure_mesh_t ib_force_structure_model_get_next(ib_force_structure_model_t handle, vertex_t* velocity, int n, double dt);

#ifdef __cplusplus
}
#endif