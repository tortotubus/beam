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
   * @brief Struct for representing the @ref ELFF::Models::IBMesh
   */
  typedef struct
  {
    int n;
    vertex_t* position;
    vertex_t* velocity;
    vertex_t* forces;
  } ib_mesh_t;

  /**
   * @memberof ib_mesh_t
   */
  void ib_mesh_free(ib_mesh_t* mesh);

#ifdef __cplusplus
}
#endif