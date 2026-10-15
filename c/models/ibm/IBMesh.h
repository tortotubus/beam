#pragma once

#ifdef __cplusplus
extern "C"
{
#endif

  /**
   * @name vertex_t
   * @brief Coordinate container used by the C immersed-boundary API.
   * @{
   */
  /**
   * @brief Helper class for storing coordinates
   */
  typedef struct
  {
    double x, y, z;
  } vertex_t;
  /** @} */

  /**
   * @name ib_mesh_t
   * @brief Plain C representation of @ref ELFF::Models::IBMesh.
   * @{
   */
  /**
   * @brief Struct for representing the @ref ELFF::Models::IBMesh
   */
  typedef struct
  {
    int n;
    vertex_t* position;
    vertex_t* velocity;
    vertex_t* forces;
    double* measure;
  } ib_mesh_t;

  /**
   * @brief Releases heap-allocated arrays owned by an `ib_mesh_t`.
   *
   * @param mesh Mesh structure to free
   */
  void ib_mesh_free(ib_mesh_t* mesh);
  /** @} */

#ifdef __cplusplus
}
#endif
