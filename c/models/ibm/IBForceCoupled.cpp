#include <stdlib.h>

#include "c/models/ibm/IBForceCoupled.h"
#include "models/ibm/IBForceCoupled.hpp"

using namespace ELFF::Models;

extern "C"
{

  int ib_force_coupled_get_number_of_nodes(
    ib_force_coupled_t handle)
  {
    IBForceCoupled* model = reinterpret_cast<IBForceCoupled*>(handle);
    return static_cast<int>(model->GetNumberOfPoints());
  }

  ib_mesh_t ib_force_coupled_get_current(
    ib_force_coupled_t handle)
  {
    IBForceCoupled* model = reinterpret_cast<IBForceCoupled*>(handle);
    IBMesh& mesh = model->GetCurrent();

    std::vector<IBMesh::IBVertex>& position = mesh.GetPoints();
    std::vector<IBMesh::IBVertex>& velocity = mesh.GetVelocity();
    int nm = mesh.GetNumberOfPoints();

    ib_mesh_t mesh_str = {
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

  ib_mesh_t ib_force_coupled_get_next(
    ib_force_coupled_t handle,
    vertex_t* force,
    int nv,
    double dt)
  {
    IBForceCoupled* model = reinterpret_cast<IBForceCoupled*>(handle);

    std::vector<IBMesh::IBVertex> force_arr(nv);

    for (int i = 0; i < nv; ++i) {
      force_arr[i].x = force[i].x;
      force_arr[i].y = force[i].y;
      force_arr[i].z = force[i].z;
    }

    IBMesh& mesh = model->GetNext(force_arr, dt);
    std::vector<IBMesh::IBVertex>& position = mesh.GetPoints();
    std::vector<IBMesh::IBVertex>& velocity = mesh.GetVelocity();
    int nm = mesh.GetNumberOfPoints();

    ib_mesh_t mesh_str = {
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
