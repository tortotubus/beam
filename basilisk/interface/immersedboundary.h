#pragma once

// #include "vertex.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /*
    Wrappers for beam::IBStructureMesh
  */

  typedef struct {
    double x, y, z;
  } vertex_t;

  typedef struct {
    int n;
    vertex_t *points;
    vertex_t *forces;
  } ib_structure_mesh_t;

  void ib_structure_mesh_free(ib_structure_mesh_t *mesh);

  /*
    Wrappers for beam::IBStructureModel
  */

  typedef void* ib_structure_model_t;

  int ib_structure_model_get_dimension       (ib_structure_model_t handle);
  int ib_structure_model_get_number_of_nodes (ib_structure_model_t handle);
  
  ib_structure_mesh_t ib_structure_model_get_current  (ib_structure_model_t handle);
  ib_structure_mesh_t ib_structure_model_get_midpoint (ib_structure_model_t handle, vertex_t* velocity, int n, double dt);
  ib_structure_mesh_t ib_structure_model_get_next     (ib_structure_model_t handle, vertex_t* velocity, int n, double dt);

#ifdef __cplusplus
}
#endif