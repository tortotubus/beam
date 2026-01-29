#include "immersedboundary.h"
#include "models/ibm/immersedboundary.hpp"
#include <stdlib.h>

// using ELFF::IBForceCoupledStructureModel;
// using ELFF::IBStructureMesh;
// using ELFF::IBVelocityCoupledStructureModel;

using namespace ELFF;

extern "C"
{

  /* ==================================================================================================================
   * Wrappers for ELFF::IBStructureMesh
   * ==================================================================================================================
   */

  /**
   * @brief
   */
  void ib_structure_mesh_free(ib_structure_mesh_t* mesh)
  {
    free(mesh->position);
    free(mesh->velocity);
    free(mesh->forces);
  }

  /* ==================================================================================================================
   * Wrappers for ELFF::IBVelocityCoupledStructureModel
   * ==================================================================================================================
   */

  /**
   * @brief
   */
  int ib_velocity_structure_model_get_number_of_nodes(
    ib_velocity_structure_model_t handle)
  {
    IBVelocityCoupledStructureModel* model =
      reinterpret_cast<IBVelocityCoupledStructureModel*>(handle);
    return static_cast<int>(model->GetNumberOfPoints());
  }

  /**
   * @brief
   */
  ib_structure_mesh_t ib_velocity_structure_model_get_current(
    ib_velocity_structure_model_t handle)
  {
    IBVelocityCoupledStructureModel* model =
      reinterpret_cast<IBVelocityCoupledStructureModel*>(handle);
    IBStructureMesh& mesh = model->GetCurrent();

    std::vector<IBStructureMesh::IBVertex>& position = mesh.GetPoints();
    std::vector<IBStructureMesh::IBVertex>& forces = mesh.GetForces();
    int nm = mesh.GetNumberOfPoints();

    ib_structure_mesh_t mesh_str = {
      .n = nm, .position = NULL, .velocity = NULL, .forces = NULL
    };

    mesh_str.position = (vertex_t*)calloc(nm, sizeof(vertex_t));
    mesh_str.velocity = (vertex_t*)calloc(nm, sizeof(vertex_t));
    mesh_str.forces = (vertex_t*)calloc(nm, sizeof(vertex_t));

    for (int i = 0; i < nm; ++i) {
      mesh_str.position[i].x = position[i].x;
      mesh_str.position[i].y = position[i].y;
      mesh_str.position[i].z = position[i].z;
      mesh_str.forces[i].x = forces[i].x;
      mesh_str.forces[i].y = forces[i].y;
      mesh_str.forces[i].z = forces[i].z;
    }

    return mesh_str;
  }

  /**
   * @brief
   */
  ib_structure_mesh_t ib_velocity_structure_model_get_midpoint(
    ib_velocity_structure_model_t handle,
    vertex_t* velocity,
    int nv,
    double dt)
  {
    // Get the pointer from the handle
    IBVelocityCoupledStructureModel* model =
      reinterpret_cast<IBVelocityCoupledStructureModel*>(handle);

    // Pack the raw C array into ELFF::Array<ELFF::fem:Vertex>
    // Array<Vertex> velocity_arr(nv);
    std::vector<IBStructureMesh::IBVertex> velocity_arr(nv);

    for (int i = 0; i < nv; ++i) {
      velocity_arr[i].x = velocity[i].x;
      velocity_arr[i].y = velocity[i].y;
      velocity_arr[i].z = velocity[i].z;
    }

    // Call GetMidpoint() and receive a reference to the proteced member
    IBStructureMesh& mesh = model->GetMidpoint(velocity_arr, dt);
    std::vector<IBStructureMesh::IBVertex>& position = mesh.GetPoints();
    std::vector<IBStructureMesh::IBVertex>& forces = mesh.GetForces();
    int nn = mesh.GetNumberOfPoints();

    // Pack the IBStructureMesh into our C struct
    ib_structure_mesh_t mesh_str = {
      .n = nn, .position = NULL, .velocity = NULL, .forces = NULL
    };

    mesh_str.position = (vertex_t*)calloc(nn, sizeof(vertex_t));
    mesh_str.velocity = (vertex_t*)calloc(nn, sizeof(vertex_t));
    mesh_str.forces = (vertex_t*)calloc(nn, sizeof(vertex_t));

    for (int i = 0; i < nn; ++i) {
      mesh_str.position[i].x = position[i].x;
      mesh_str.position[i].y = position[i].y;
      mesh_str.position[i].z = position[i].z;
      mesh_str.forces[i].x = forces[i].x;
      mesh_str.forces[i].y = forces[i].y;
      mesh_str.forces[i].z = forces[i].z;
    }

    return mesh_str;
  }

  /**
   * @brief
   */
  ib_structure_mesh_t ib_velocity_structure_model_get_next(
    ib_velocity_structure_model_t handle,
    vertex_t* velocity,
    int nv,
    double dt)
  {
    // Get the pointer from the handle
    IBVelocityCoupledStructureModel* model =
      reinterpret_cast<IBVelocityCoupledStructureModel*>(handle);

    // Pack the raw C array
    std::vector<IBStructureMesh::IBVertex> velocity_arr(nv);

    for (int i = 0; i < nv; ++i) {
      velocity_arr[i].x = velocity[i].x;
      velocity_arr[i].y = velocity[i].y;
      velocity_arr[i].z = velocity[i].z;
    }

    // Call GetMidpoint() and receive a reference to the proteced member
    IBStructureMesh& mesh = model->GetNext(velocity_arr, dt);
    std::vector<IBStructureMesh::IBVertex>& position = mesh.GetPoints();
    std::vector<IBStructureMesh::IBVertex>& forces = mesh.GetForces();
    int nm = mesh.GetNumberOfPoints();

    // Pack the IBStructureMesh into our C struct
    ib_structure_mesh_t mesh_str = {
      .n = nm, .position = NULL, .velocity = NULL, .forces = NULL
    };

    mesh_str.position = (vertex_t*)calloc(nm, sizeof(vertex_t));
    mesh_str.velocity = (vertex_t*)calloc(nm, sizeof(vertex_t));
    mesh_str.forces = (vertex_t*)calloc(nm, sizeof(vertex_t));

    for (int i = 0; i < nm; ++i) {
      mesh_str.position[i].x = position[i].x;
      mesh_str.position[i].y = position[i].y;
      mesh_str.position[i].z = position[i].z;
      mesh_str.forces[i].x = forces[i].x;
      mesh_str.forces[i].y = forces[i].y;
      mesh_str.forces[i].z = forces[i].z;
    }

    return mesh_str;
  }

  /* ================================================================================================================== 
   * Wrappers for ELFF::IBForceCoupledStructureModel
   * ==================================================================================================================
   */

  /**
   * @brief
   */
  int ib_force_structure_model_get_number_of_nodes(
    ib_force_structure_model_t handle)
  {
    IBForceCoupledStructureModel* model =
      reinterpret_cast<IBForceCoupledStructureModel*>(handle);

      int n = model->GetNumberOfPoints();

    return static_cast<int>(model->GetNumberOfPoints());
  }

  /**
   * @brief
   */
  ib_structure_mesh_t ib_force_structure_model_get_current(
    ib_force_structure_model_t handle)
  {
    IBForceCoupledStructureModel* model =
      reinterpret_cast<IBForceCoupledStructureModel*>(handle);
    IBStructureMesh& mesh = model->GetCurrent();

    std::vector<IBStructureMesh::IBVertex>& position = mesh.GetPoints();
    std::vector<IBStructureMesh::IBVertex>& velocity = mesh.GetVelocity();

    int nm = mesh.GetNumberOfPoints();

    ib_structure_mesh_t mesh_str = {
      .n = nm, .position = NULL, .velocity = NULL, .forces = NULL
    };

    mesh_str.position = (vertex_t*)calloc(nm, sizeof(vertex_t));
    mesh_str.velocity = (vertex_t*)calloc(nm, sizeof(vertex_t));
    mesh_str.forces = (vertex_t*)calloc(nm, sizeof(vertex_t));

    for (int i = 0; i < nm; ++i) {
      mesh_str.position[i].x = position[i].x;
      mesh_str.position[i].y = position[i].y;
      mesh_str.position[i].z = position[i].z;
      mesh_str.velocity[i].x = velocity[i].x;
      mesh_str.velocity[i].y = velocity[i].y;
      mesh_str.velocity[i].z = velocity[i].z;
    }

    return mesh_str;
  }

  /**
   * @brief
   */
  ib_structure_mesh_t ib_force_structure_model_get_next(
    ib_force_structure_model_t handle,
    vertex_t* force,
    int nv,
    double dt)
  {
    // Get the pointer from the handle
    IBForceCoupledStructureModel* model =
      reinterpret_cast<IBForceCoupledStructureModel*>(handle);

    // Pack the raw C array into ELFF::Array<ELFF::fem:Vertex>
    std::vector<IBStructureMesh::IBVertex> force_arr(nv);

    for (int i = 0; i < nv; ++i) {
      force_arr[i].x = force[i].x;
      force_arr[i].y = force[i].y;
      force_arr[i].z = force[i].z;
    }

    // Call GetNext() and receive a reference to the proteced member
    IBStructureMesh& mesh = model->GetNext(force_arr, dt);
    std::vector<IBStructureMesh::IBVertex>& position = mesh.GetPoints();
    std::vector<IBStructureMesh::IBVertex>& velocity = mesh.GetVelocity();
    int nm = mesh.GetNumberOfPoints();

    // Pack the IBStructureMesh into our C struct
    ib_structure_mesh_t mesh_str = {
      .n = nm, .position = NULL, .velocity = NULL, .forces = NULL
    };

    mesh_str.position = (vertex_t*)calloc(nm, sizeof(vertex_t));
    mesh_str.velocity = (vertex_t*)calloc(nm, sizeof(vertex_t));
    mesh_str.forces = (vertex_t*)calloc(nm, sizeof(vertex_t));

    for (int i = 0; i < nm; ++i) {
      mesh_str.position[i].x = position[i].x;
      mesh_str.position[i].y = position[i].y;
      mesh_str.position[i].z = position[i].z;
      mesh_str.velocity[i].x = velocity[i].x;
      mesh_str.velocity[i].y = velocity[i].y;
      mesh_str.velocity[i].z = velocity[i].z;
    }

    return mesh_str;
  }

} // extern "C"
