#pragma once

#ifdef __cplusplus
extern "C"
{
#endif

  /*
    Wrappers for beam::fem::Vertex
  */
  typedef void* vertex_handle_t;

  // Constructor functions
  // vertex_handle_t vertex_create(double x, double y, double z);
  
  // Destructor
  // void vertex_destroy(vertex_handle_t handle);
  
  // Accessor functions
  // double vertex_get_component(vertex_t handle, int i);
  // void vertex_set_component(vertex_t handle, int i, double value);

#ifdef __cplusplus
}
#endif