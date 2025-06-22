#include "immersedboundary.h"
#include "../../models/immersedboundary.hpp"
#include <stdlib.h>

using beam::Array;
using beam::IBStructureMesh;
using beam::IBStructureModel;
using beam::fem::Vertex;

extern "C"
{

  /*
    Wrappers for beam::IBStructureMesh
  */

  void ib_structure_mesh_free(ib_structure_mesh_t* mesh)
  {
    free(mesh->points);
    free(mesh->forces);
  }

  /*
      Wrappers for beam::IBStructureModel
  */

  int ib_structure_model_get_dimension(ib_structure_model_t handle)
  {
    IBStructureModel* model = reinterpret_cast<IBStructureModel*>(handle);
    return model->GetDimension();
  }
  
  int ib_structure_model_get_number_of_nodes(ib_structure_model_t handle)
  {
    IBStructureModel* model = reinterpret_cast<IBStructureModel*>(handle);
    return model->GetNumberOfPoints();
  }

  ib_structure_mesh_t ib_structure_model_get_current(
    ib_structure_model_t handle)
  {
    IBStructureModel* model = reinterpret_cast<IBStructureModel*>(handle);
    IBStructureMesh& mesh = model->GetCurrent();

    Array<Vertex>& points = mesh.GetPoints();
    Array<Vertex>& forces = mesh.GetForces();
    int nm = mesh.GetNumberOfPoints();

    ib_structure_mesh_t mesh_str = { .n = nm, .points = NULL, .forces = NULL };

    mesh_str.points = (vertex_t*)calloc(nm, sizeof(vertex_t));
    mesh_str.forces = (vertex_t*)calloc(nm, sizeof(vertex_t));

    for (int i = 0; i < nm; ++i) {
      mesh_str.points[i].x = points[i](0);
      mesh_str.points[i].y = points[i](1);
      mesh_str.points[i].z = points[i](2);
      mesh_str.forces[i].x = forces[i](0);
      mesh_str.forces[i].y = forces[i](1);
      mesh_str.forces[i].z = forces[i](2);
    }

    return mesh_str;
  }

  ib_structure_mesh_t ib_structure_model_get_midpoint(
    ib_structure_model_t handle,
    vertex_t* velocity,
    int nv,
    double dt)
  {
    // Get the pointer from the handle
    IBStructureModel* model = reinterpret_cast<IBStructureModel*>(handle);

    // Pack the raw C array into beam::Array<beam::fem:Vertex>
    Array<Vertex> velocity_arr(nv);

    for (int i = 0; i < nv; ++i) {
      velocity_arr[i](0) = velocity[i].x;
      velocity_arr[i](1) = velocity[i].y;
      velocity_arr[i](2) = velocity[i].z;
    }

    // Call GetMidpoint() and receive a reference to the proteced member
    IBStructureMesh& mesh = model->GetMidpoint(velocity_arr, dt);
    Array<Vertex>& points = mesh.GetPoints();
    Array<Vertex>& forces = mesh.GetForces();
    int nm = mesh.GetNumberOfPoints();

    // Pack the IBStructureMesh into our C struct
    ib_structure_mesh_t mesh_str = { .n = nm, .points = NULL, .forces = NULL };

    mesh_str.points = (vertex_t*)calloc(nm, sizeof(vertex_t));
    mesh_str.forces = (vertex_t*)calloc(nm, sizeof(vertex_t));

    for (int i = 0; i < nm; ++i) {
      mesh_str.points[i].x = points[i](0);
      mesh_str.points[i].y = points[i](1);
      mesh_str.points[i].z = points[i](2);
      mesh_str.forces[i].x = forces[i](0);
      mesh_str.forces[i].y = forces[i](1);
      mesh_str.forces[i].z = forces[i](2);
    }

    return mesh_str;
  }

  ib_structure_mesh_t ib_structure_model_get_next(ib_structure_model_t handle,
                                                  vertex_t* velocity,
                                                  int nv,
                                                  double dt)
  {
    // Get the pointer from the handle
    IBStructureModel* model = reinterpret_cast<IBStructureModel*>(handle);

    // Pack the raw C array into beam::Array<beam::fem:Vertex>
    Array<Vertex> velocity_arr(nv);

    for (int i = 0; i < nv; ++i) {
      velocity_arr[i](0) = velocity[i].x;
      velocity_arr[i](1) = velocity[i].y;
      velocity_arr[i](2) = velocity[i].z;
    }

    // Call GetMidpoint() and receive a reference to the proteced member
    IBStructureMesh& mesh = model->GetMidpoint(velocity_arr, dt);
    Array<Vertex>& points = mesh.GetPoints();
    Array<Vertex>& forces = mesh.GetForces();
    int nm = mesh.GetNumberOfPoints();

    // Pack the IBStructureMesh into our C struct
    ib_structure_mesh_t mesh_str = { .n = nm, .points = NULL, .forces = NULL };

    mesh_str.points = (vertex_t*)calloc(nm, sizeof(vertex_t));
    mesh_str.forces = (vertex_t*)calloc(nm, sizeof(vertex_t));

    for (int i = 0; i < nm; ++i) {
      mesh_str.points[i].x = points[i](0);
      mesh_str.points[i].y = points[i](1);
      mesh_str.points[i].z = points[i](2);
      mesh_str.forces[i].x = forces[i](0);
      mesh_str.forces[i].y = forces[i](1);
      mesh_str.forces[i].z = forces[i](2);
    }

    return mesh_str;
  }

} // extern "C"
