
#include "models/beam/IBEulerBeam.h"
#include "models/ibm/IBForceCoupled.h"
#include "library/ibm/IBMeshModel.h"

/**
 * ELFF Force-Coupled Operations
 */

int elff_fc_node_count (void* ctx);
void elff_fc_sync (void* ctx, void* mesh);
void elff_fc_advance (void* ctx, void* mesh, double dt);

/**
 * @brief
 */
int elff_fc_node_count (void* ctx) {
  ib_force_coupled_t fc_ptr = (ib_force_coupled_t) ctx;
  return ib_force_coupled_get_number_of_nodes (fc_ptr);
}

/**
 * @brief 
 */
void elff_fc_sync (void* ctx, void* mesh) {
  // Basilisk pointers
  IBMesh* ib_mesh = (IBMesh*) mesh;
  int ib_nodes_count = ib_mesh->nodes.size;
  IBNode** ib_nodes = ib_mesh->nodes.ptrs;

  // ELFF Pointers
  ib_force_coupled_t fc_ptr = (ib_force_coupled_t) ctx;
  ib_mesh_t elff_mesh = ib_force_coupled_get_current (fc_ptr);

  for (int ni = 0; ni < ib_nodes_count; ni++) {
    foreach_dimension () {
      ib_nodes[ni]->pos.x = elff_mesh.position[ni].x;
      ib_nodes[ni]->vel.x = elff_mesh.velocity[ni].x;
    }
  }

  ib_mesh_free (&elff_mesh);
}

/**
 * @brief 
 */
void elff_fc_advance (void* ctx, void* mesh, double dt) {
  // Basilisk pointers
  IBMesh* ib_mesh = (IBMesh*) mesh;
  int ib_nodes_count = ib_mesh->nodes.size;
  IBNode** ib_nodes = ib_mesh->nodes.ptrs;

  // Pack nodal forces into vertex_t array
  vertex_t* forces = calloc (ib_nodes_count, sizeof (vertex_t));
  for (int ni = 0; ni < ib_nodes_count; ni++) {
    foreach_dimension () {
      forces[ni].x = ib_nodes[ni]->f.x;
    }
  }

  // ELFF Pointers
  ib_force_coupled_t fc_ptr = (ib_force_coupled_t) ctx;
  ib_mesh_t elff_mesh =
    ib_force_coupled_get_next (fc_ptr, forces, ib_nodes_count, dt);

  for (int ni = 0; ni < ib_nodes_count; ni++) {
    foreach_dimension () {
      ib_nodes[ni]->pos.x = elff_mesh.position[ni].x;
      ib_nodes[ni]->vel.x = elff_mesh.velocity[ni].x;
    }
  }

  free (forces);
  ib_mesh_free (&elff_mesh);
}

/*
* Euler-Bernouli beam
*/

IBMeshModel elff_beam_new (double length, double EI, double mu, int nodes, double r_penalty);
IBMeshModel elff_beam_new_theta (double length, double EI, double mu, int nodes, double r_penalty, double theta);
void elff_beam_destroy (void* ctx);


/**
 * @brief 
 */
IBMeshModel elff_beam_new (
  double length, double EI, double mu, int nodes, double r_penalty) {
  ib_beam_t beam_ptr = ib_beam_new (length, EI, mu, nodes, r_penalty);
  IBMeshModel ib_model = ibmeshmodel_force_coupled_init ();

  ib_model.ctx = beam_ptr;
  ib_model.force_ops->node_count = elff_fc_node_count;
  ib_model.force_ops->sync = elff_fc_sync;
  ib_model.force_ops->advance = elff_fc_advance;
  ib_model.force_ops->destroy = elff_beam_destroy;

  return ib_model;
}

IBMeshModel elff_beam_new_theta (
  double length, double EI, double mu, int nodes, double r_penalty, double theta) {
  ib_beam_t beam_ptr = ib_beam_new_theta (length, EI, mu, nodes, r_penalty, theta);
  IBMeshModel ib_model = ibmeshmodel_force_coupled_init ();

  ib_model.ctx = beam_ptr;
  ib_model.force_ops->node_count = elff_fc_node_count;
  ib_model.force_ops->sync = elff_fc_sync;
  ib_model.force_ops->advance = elff_fc_advance;
  ib_model.force_ops->destroy = elff_beam_destroy;

  return ib_model;
}

/**
 * @brief 
 */
void elff_beam_destroy (void* ctx) {
  ib_beam_t handle = (ib_beam_t) ctx;
  ib_beam_destroy (handle);
}
