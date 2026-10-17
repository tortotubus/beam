#if _MPI
#ifndef ELFF_USE_MPI
#define ELFF_USE_MPI 1
#endif
#endif

#include "elff/c/models/beam/IBEulerBeam.h"
#include "elff/c/models/beam/IBEulerBeamADDM.h"
#include "elff/c/models/beam/IBEulerBeamGGL.h"
#include "elff/c/models/beam/IBEulerBeamHuang.h"
#include "elff/c/models/beam/IBEulerBeamPenalty.h"
#include "elff/c/models/capsule/IBCapsule.h"
#include "elff/c/models/ibm/IBForceCoupled.h"
#include "elff/c/models/ibm/IBVelocityCoupled.h"
#include "elff/c/models/rigidbody/IBRigidBody.h"

#include "library/ibm/IBMeshModel.h"
#include "library/elff/runtime.h"

static inline ib_euler_beam_bcs_t elff_euler_beam_bcs_types (int left_type,
                                                             int right_type) {
  ib_euler_beam_bcs_t bcs = {0};
  bcs.end[0] = IB_EULER_BEAM_BC_LEFT;
  bcs.end[1] = IB_EULER_BEAM_BC_RIGHT;
  bcs.type[0] = left_type;
  bcs.type[1] = right_type;
  return bcs;
}

static inline ib_euler_beam_bcs_t elff_euler_beam_bcs_default () {
  ib_euler_beam_bcs_t bcs = {0};
  bcs.end[0] = IB_EULER_BEAM_BC_LEFT;
  bcs.end[1] = IB_EULER_BEAM_BC_RIGHT;
  bcs.type[0] = IB_EULER_BEAM_BC_FREE;
  bcs.type[1] = IB_EULER_BEAM_BC_SIMPLE;
  return bcs;
}

static inline ib_euler_beam_bcs_t elff_euler_beam_bcs_theta_pin (coord s0) {
  ib_euler_beam_bcs_t bcs = elff_euler_beam_bcs_default ();
  bcs.vals[1].position.x = s0.x;
  bcs.vals[1].position.y = s0.y;
  bcs.vals[1].position.z = s0.z;
  return bcs;
}

// ============================================================================
// ELFF Velocity-Coupled Operations
// ============================================================================

int elff_vc_node_count (void* ctx);
int elff_vc_sync (void* ctx, void* mesh);
int elff_vc_midpoint (void* ctx, void* mesh, double dt);
int elff_vc_advance (void* ctx, void* mesh, double dt);

static inline vertex_t elff_vertex_from_coord (coord value) {
  vertex_t out = {0};
  out.x = value.x;
  out.y = value.y;
  out.z = value.z;
  return out;
}

static inline IBMeshModel elff_velocity_coupled_model (void* ctx) {
  IBMeshModel ib_model = ibmeshmodel_velocity_coupled_init ();
  ib_model.ctx = ctx;
  ib_model.velocity_ops->node_count = elff_vc_node_count;
  ib_model.velocity_ops->sync = elff_vc_sync;
  ib_model.velocity_ops->midpoint = elff_vc_midpoint;
  ib_model.velocity_ops->advance = elff_vc_advance;
  ib_model.velocity_ops->destroy = NULL;
  return ib_model;
}

/**
 * @brief
 */
int elff_vc_node_count (void* ctx) {
  ib_velocity_coupled_t vc_ptr = (ib_velocity_coupled_t) ctx;
  return ib_velocity_coupled_get_number_of_nodes (vc_ptr);
}

/**
 * @brief
 */
trace int elff_vc_sync (void* ctx, void* mesh) {
  IBMesh* ib_mesh = (IBMesh*) mesh;
  int ib_nodes_count = ib_mesh->nodes.size;
  IBNode** ib_nodes = ib_mesh->nodes.ptrs;

  ib_velocity_coupled_t vc_ptr = (ib_velocity_coupled_t) ctx;
  ib_mesh_t elff_mesh = ib_velocity_coupled_get_current (vc_ptr);

  if (elff_mesh.n != ib_nodes_count) {
    ib_mesh_free (&elff_mesh);
    return -1;
  }

  for (int ni = 0; ni < ib_nodes_count; ni++) {
    IBNode* node = ib_nodes[ni];
    const double nodal_measure =
      elff_mesh.measure ? elff_mesh.measure[ni] : 1.0;
    foreach_dimension () {
      ibval (npos.x) = elff_mesh.position[ni].x;
      ibval (nvel.x) = elff_mesh.velocity[ni].x;
      ibval (nforce.x) = elff_mesh.forces ? elff_mesh.forces[ni].x : 0.;
    }
    ibval (nweight) = nodal_measure;
  }

  ib_mesh_free (&elff_mesh);
  return 0;
}

/**
 * @brief
 */
trace int elff_vc_midpoint (void* ctx, void* mesh, double dt) {
  IBMesh* ib_mesh = (IBMesh*) mesh;
  int ib_nodes_count = ib_mesh->nodes.size;
  IBNode** ib_nodes = ib_mesh->nodes.ptrs;

  vertex_t* velocity = calloc (ib_nodes_count, sizeof (vertex_t));
  for (int ni = 0; ni < ib_nodes_count; ni++) {
    IBNode* node = ib_nodes[ni];
    foreach_dimension () {
      velocity[ni].x = ibval (nvel.x);
    }
  }

  ib_velocity_coupled_t vc_ptr = (ib_velocity_coupled_t) ctx;
  ib_mesh_t elff_mesh =
    ib_velocity_coupled_get_midpoint (vc_ptr, velocity, ib_nodes_count, dt);

  free (velocity);

  if (elff_mesh.n != ib_nodes_count) {
    ib_mesh_free (&elff_mesh);
    return -1;
  }

  for (int ni = 0; ni < ib_nodes_count; ni++) {
    IBNode* node = ib_nodes[ni];
    const double nodal_measure =
      elff_mesh.measure ? elff_mesh.measure[ni] : 1.0;
    foreach_dimension () {
      ibval (npos.x) = elff_mesh.position[ni].x;
      ibval (nvel.x) = elff_mesh.velocity[ni].x;
      ibval (nforce.x) = elff_mesh.forces ? elff_mesh.forces[ni].x : 0.;
    }
    ibval (nweight) = nodal_measure;
  }

  ib_mesh_free (&elff_mesh);
  return 0;
}

/**
 * @brief
 */
trace int elff_vc_advance (void* ctx, void* mesh, double dt) {
  IBMesh* ib_mesh = (IBMesh*) mesh;
  int ib_nodes_count = ib_mesh->nodes.size;
  IBNode** ib_nodes = ib_mesh->nodes.ptrs;

  vertex_t* velocity = calloc (ib_nodes_count, sizeof (vertex_t));
  for (int ni = 0; ni < ib_nodes_count; ni++) {
    IBNode* node = ib_nodes[ni];
    foreach_dimension () {
      velocity[ni].x = ibval (nvel.x);
    }
  }

  ib_velocity_coupled_t vc_ptr = (ib_velocity_coupled_t) ctx;
  ib_mesh_t elff_mesh =
    ib_velocity_coupled_get_next (vc_ptr, velocity, ib_nodes_count, dt);

  free (velocity);

  if (elff_mesh.n != ib_nodes_count) {
    ib_mesh_free (&elff_mesh);
    return -1;
  }

  for (int ni = 0; ni < ib_nodes_count; ni++) {
    IBNode* node = ib_nodes[ni];
    const double nodal_measure =
      elff_mesh.measure ? elff_mesh.measure[ni] : 1.0;
    foreach_dimension () {
      ibval (npos.x) = elff_mesh.position[ni].x;
      ibval (nvel.x) = elff_mesh.velocity[ni].x;
      ibval (nforce.x) = elff_mesh.forces ? elff_mesh.forces[ni].x : 0.;
    }
    ibval (nweight) = nodal_measure;
  }

  ib_mesh_free (&elff_mesh);
  return 0;
}

// ============================================================================
// ELFF Force-Coupled Operations
// ============================================================================

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
trace void elff_fc_sync (void* ctx, void* mesh) {

  // Basilisk pointers
  IBMesh* ib_mesh = (IBMesh*) mesh;
  int ib_nodes_count = ib_mesh->nodes.size;
  IBNode** ib_nodes = ib_mesh->nodes.ptrs;

  // ELFF Pointers
  ib_force_coupled_t fc_ptr = (ib_force_coupled_t) ctx;
  ib_mesh_t elff_mesh = ib_force_coupled_get_current (fc_ptr);

  for (int ni = 0; ni < ib_nodes_count; ni++) {
    IBNode* node = ib_nodes[ni];
    const double nodal_measure =
      elff_mesh.measure ? elff_mesh.measure[ni] : 1.0;
    foreach_dimension () {
      ibval (npos.x) = elff_mesh.position[ni].x;
      ibval (nvel.x) = elff_mesh.velocity[ni].x;
    }
    ibval (nweight) = nodal_measure;
  }

  ib_mesh_free (&elff_mesh);
}

/**
 * @brief
 */
trace void elff_fc_advance (void* ctx, void* mesh, double dt) {
  // Basilisk pointers
  IBMesh* ib_mesh = (IBMesh*) mesh;
  int ib_nodes_count = ib_mesh->nodes.size;
  IBNode** ib_nodes = ib_mesh->nodes.ptrs;

  // Pack nodal forces into vertex_t array
  vertex_t* forces = calloc (ib_nodes_count, sizeof (vertex_t));
  for (int ni = 0; ni < ib_nodes_count; ni++) {
    IBNode* node = ib_nodes[ni];
    foreach_dimension () {
      forces[ni].x = ibval (nforce.x);
    }
  }

  // ELFF Pointers
  ib_force_coupled_t fc_ptr = (ib_force_coupled_t) ctx;
  ib_mesh_t elff_mesh = ib_force_coupled_get_next (fc_ptr, forces, ib_nodes_count, dt);

  for (int ni = 0; ni < ib_nodes_count; ni++) {
    IBNode* node = ib_nodes[ni];
    const double nodal_measure =
      elff_mesh.measure ? elff_mesh.measure[ni] : 1.0;
    foreach_dimension () {
      ibval (npos.x) = elff_mesh.position[ni].x;
      ibval (nvel.x) = elff_mesh.velocity[ni].x;
    }
    ibval (nweight) = nodal_measure;
  }

  free (forces);
  ib_mesh_free (&elff_mesh);
}

static inline IBMeshModel elff_force_coupled_model (void* ctx) {
  IBMeshModel model = ibmeshmodel_force_coupled_init ();
  model.ctx = ctx;
  model.force_ops->node_count = elff_fc_node_count;
  model.force_ops->sync = elff_fc_sync;
  model.force_ops->advance = elff_fc_advance;
  model.force_ops->destroy = NULL;
  return model;
}

/** Install a non-owning Basilisk adapter for one runtime-owned view. */
static inline IBMeshModel elff_model_from_view (ib_model_view_t view) {
  if (view.kind == IB_MODEL_VIEW_FORCE_COUPLED)
    return elff_force_coupled_model (view.handle);
  if (view.kind == IB_MODEL_VIEW_VELOCITY_COUPLED)
    return elff_velocity_coupled_model (view.handle);
  IBMeshModel invalid = {0};
  return invalid;
}

// ============================================================================
// ELFF EulerBernoulli Beam
// ============================================================================

IBMeshModel elff_euler_beam_new (double length,
                                 double EI,
                                 double mu,
                                 int nodes,
                                 double r_penalty,
                                 ib_euler_beam_bcs_t bcs,
                                 coord s0);

IBMeshModel elff_euler_beam_new_theta (double length,
                                       double EI,
                                       double mu,
                                       int nodes,
                                       double r_penalty,
                                       double theta,
                                       ib_euler_beam_bcs_t bcs,
                                       coord s0);

IBMeshModel elff_euler_beam_penalty_new (double length,
                                         double EI,
                                         double mu,
                                         int nodes,
                                         double r_penalty,
                                         ib_euler_beam_bcs_t bcs,
                                         coord s0);
IBMeshModel elff_euler_beam_penalty_new_theta (double length,
                                               double EI,
                                               double mu,
                                               int nodes,
                                               double r_penalty,
                                               double theta,
                                               ib_euler_beam_bcs_t bcs,
                                               coord s0);
IBMeshModel elff_euler_beam_ggl_new (double length,
                                     double EI,
                                     double mu,
                                     int nodes,
                                     double r_penalty,
                                     ib_euler_beam_bcs_t bcs,
                                     coord s0);
IBMeshModel elff_euler_beam_ggl_new_theta (double length,
                                           double EI,
                                           double mu,
                                           int nodes,
                                           double r_penalty,
                                           double theta,
                                           ib_euler_beam_bcs_t bcs,
                                           coord s0);
IBMeshModel elff_euler_beam_addm_new (double length,
                                      double EI,
                                      double mu,
                                      int nodes,
                                      double r_penalty,
                                      ib_euler_beam_bcs_t bcs,
                                      coord s0);
IBMeshModel elff_euler_beam_addm_new_theta (double length,
                                            double EI,
                                            double mu,
                                            int nodes,
                                            double r_penalty,
                                            double theta,
                                            ib_euler_beam_bcs_t bcs,
                                            coord s0);
IBMeshModel elff_euler_beam_huang_new (double length,
                                       double EI,
                                       double mu,
                                       int nodes,
                                       ib_euler_beam_bcs_t bcs,
                                       coord s0);
IBMeshModel elff_euler_beam_huang_new_theta (double length,
                                             double EI,
                                             double mu,
                                             int nodes,
                                             double theta,
                                             ib_euler_beam_bcs_t bcs,
                                             coord s0);
IBMeshModel elff_euler_beam_huang_new_direction (double length,
                                                 double EI,
                                                 double mu,
                                                 int nodes,
                                                 coord direction,
                                                 ib_euler_beam_bcs_t bcs,
                                                 coord s0);

IBMeshModel elff_capsule_sphere_new (double radius,
                                     coord center,
                                     int refinements,
                                     double elastic_modulus);

IBMeshModel elff_capsule_ellipsoid_new (coord radii,
                                        coord center,
                                        int refinements,
                                        double elastic_modulus);

IBMeshModel elff_capsule_biconcave_new (double radius,
                                        coord center,
                                        int refinements,
                                        double elastic_modulus);

IBMeshModel elff_rigid_body_circle_new (double radius,
                                        int point_count,
                                        double density,
                                        coord center,
                                        coord velocity,
                                        double angle,
                                        double angular_velocity_z);

IBMeshModel elff_rigid_body_fibre_new (double length,
                                       double diameter,
                                       int nodes,
                                       double linear_density,
                                       coord center,
                                       coord velocity,
                                       double q_w,
                                       coord q_xyz,
                                       coord angular_velocity_body);

IBMeshModel elff_pinned_rigid_body_circle_new (double radius,
                                               int point_count,
                                               double density,
                                               coord center,
                                               double angle);

IBMeshModel elff_pinned_rigid_body_fibre_new (double length,
                                              double diameter,
                                              int nodes,
                                              double linear_density,
                                              coord center,
                                              double q_w,
                                              coord q_xyz);

IBMeshModel elff_rigid_body_sphere_new (double radius,
                                        int point_count,
                                        double density,
                                        coord center,
                                        coord velocity,
                                        double q_w,
                                        coord q_xyz,
                                        coord angular_velocity_body);

IBMeshModel elff_pinned_rigid_body_sphere_new (double radius,
                                               int point_count,
                                               double density,
                                               coord center,
                                               double q_w,
                                               coord q_xyz);

/**
 * @brief
 */
IBMeshModel elff_euler_beam_new (double length,
                                 double EI,
                                 double mu,
                                 int nodes,
                                 double r_penalty,
                                 ib_euler_beam_bcs_t bcs = elff_euler_beam_bcs_default (),
                                 coord s0 = {0}) {
  vertex_t v0 = {s0.x, s0.y, s0.z};
  ib_euler_beam_t beam_ptr =
    ib_euler_beam_new (elff_runtime_get (), v0, bcs, length, EI, mu, nodes, r_penalty);

  IBMeshModel ib_model = ibmeshmodel_force_coupled_init ();

  ib_model.ctx = beam_ptr;
  ib_model.force_ops->node_count = elff_fc_node_count;
  ib_model.force_ops->sync = elff_fc_sync;
  ib_model.force_ops->advance = elff_fc_advance;
  ib_model.force_ops->destroy = NULL;

  return ib_model;
}

/**
 * @brief
 */
IBMeshModel elff_euler_beam_new_theta (double length,
                                       double EI,
                                       double mu,
                                       int nodes,
                                       double r_penalty,
                                       double theta,
                                       ib_euler_beam_bcs_t bcs = elff_euler_beam_bcs_default (),
                                       coord s0 = {0}) {
  vertex_t v0 = {s0.x, s0.y, s0.z};
  ib_euler_beam_t beam_ptr = ib_euler_beam_new_theta (
    elff_runtime_get (), v0, bcs, length, EI, mu, nodes, r_penalty, theta);

  IBMeshModel ib_model = ibmeshmodel_force_coupled_init ();

  ib_model.ctx = beam_ptr;
  ib_model.force_ops->node_count = elff_fc_node_count;
  ib_model.force_ops->sync = elff_fc_sync;
  ib_model.force_ops->advance = elff_fc_advance;
  ib_model.force_ops->destroy = NULL;

  return ib_model;
}

/**
 * @brief
 */
IBMeshModel elff_euler_beam_penalty_new (double length,
                                         double EI,
                                         double mu,
                                         int nodes,
                                         double r_penalty,
                                         ib_euler_beam_bcs_t bcs = elff_euler_beam_bcs_default (),
                                         coord s0 = {0}) {
  vertex_t v0 = {s0.x, s0.y, s0.z};
  ib_euler_beam_penalty_t beam_ptr = ib_euler_beam_penalty_new (
    elff_runtime_get (), v0, bcs, length, EI, mu, nodes, r_penalty);

  IBMeshModel ib_model = ibmeshmodel_force_coupled_init ();

  ib_model.ctx = beam_ptr;
  ib_model.force_ops->node_count = elff_fc_node_count;
  ib_model.force_ops->sync = elff_fc_sync;
  ib_model.force_ops->advance = elff_fc_advance;
  ib_model.force_ops->destroy = NULL;

  return ib_model;
}

/**
 * @brief
 */
IBMeshModel elff_euler_beam_penalty_new_theta (double length,
                                               double EI,
                                               double mu,
                                               int nodes,
                                               double r_penalty,
                                               double theta,
                                               ib_euler_beam_bcs_t bcs = elff_euler_beam_bcs_default (),
                                               coord s0 = {0}) {
  vertex_t v0 = {s0.x, s0.y, s0.z};
  ib_euler_beam_penalty_t beam_ptr = ib_euler_beam_penalty_new_theta (
    elff_runtime_get (), v0, bcs, length, EI, mu, nodes, r_penalty, theta);

  IBMeshModel ib_model = ibmeshmodel_force_coupled_init ();

  ib_model.ctx = beam_ptr;
  ib_model.force_ops->node_count = elff_fc_node_count;
  ib_model.force_ops->sync = elff_fc_sync;
  ib_model.force_ops->advance = elff_fc_advance;
  ib_model.force_ops->destroy = NULL;

  return ib_model;
}

/**
 * @brief
 */
IBMeshModel elff_euler_beam_ggl_new (double length,
                                     double EI,
                                     double mu,
                                     int nodes,
                                     double r_penalty,
                                     ib_euler_beam_bcs_t bcs = elff_euler_beam_bcs_default (),
                                     coord s0 = {0}) {
  vertex_t v0 = {s0.x, s0.y, s0.z};
  ib_euler_beam_ggl_t beam_ptr =
    ib_euler_beam_ggl_new (elff_runtime_get (), v0, bcs, length, EI, mu, nodes, r_penalty);

  IBMeshModel ib_model = ibmeshmodel_force_coupled_init ();

  ib_model.ctx = beam_ptr;
  ib_model.force_ops->node_count = elff_fc_node_count;
  ib_model.force_ops->sync = elff_fc_sync;
  ib_model.force_ops->advance = elff_fc_advance;
  ib_model.force_ops->destroy = NULL;

  return ib_model;
}

/**
 * @brief
 */
IBMeshModel elff_euler_beam_ggl_new_theta (double length,
                                           double EI,
                                           double mu,
                                           int nodes,
                                           double r_penalty,
                                           double theta,
                                           ib_euler_beam_bcs_t bcs = elff_euler_beam_bcs_default (),
                                           coord s0 = {0}) {
  vertex_t v0 = {s0.x, s0.y, s0.z};
  ib_euler_beam_ggl_t beam_ptr = ib_euler_beam_ggl_new_theta (
    elff_runtime_get (), v0, bcs, length, EI, mu, nodes, r_penalty, theta);

  IBMeshModel ib_model = ibmeshmodel_force_coupled_init ();

  ib_model.ctx = beam_ptr;
  ib_model.force_ops->node_count = elff_fc_node_count;
  ib_model.force_ops->sync = elff_fc_sync;
  ib_model.force_ops->advance = elff_fc_advance;
  ib_model.force_ops->destroy = NULL;

  return ib_model;
}

/**
 * @brief
 */
IBMeshModel elff_euler_beam_addm_new (double length,
                                      double EI,
                                      double mu,
                                      int nodes,
                                      double r_penalty,
                                      ib_euler_beam_bcs_t bcs = elff_euler_beam_bcs_default (),
                                      coord s0 = {0}) {
  vertex_t v0 = {s0.x, s0.y, s0.z};
  ib_euler_beam_addm_t beam_ptr = ib_euler_beam_addm_new (
    elff_runtime_get (), v0, bcs, length, EI, mu, nodes, r_penalty);

  IBMeshModel ib_model = ibmeshmodel_force_coupled_init ();

  ib_model.ctx = beam_ptr;
  ib_model.force_ops->node_count = elff_fc_node_count;
  ib_model.force_ops->sync = elff_fc_sync;
  ib_model.force_ops->advance = elff_fc_advance;
  ib_model.force_ops->destroy = NULL;

  return ib_model;
}

/**
 * @brief
 */
IBMeshModel elff_euler_beam_addm_new_theta (double length,
                                            double EI,
                                            double mu,
                                            int nodes,
                                            double r_penalty,
                                            double theta,
                                            ib_euler_beam_bcs_t bcs = elff_euler_beam_bcs_default (),
                                            coord s0 = {0}) {
  vertex_t v0 = {s0.x, s0.y, s0.z};
  ib_euler_beam_addm_t beam_ptr = ib_euler_beam_addm_new_theta (
    elff_runtime_get (), v0, bcs, length, EI, mu, nodes, r_penalty, theta);

  IBMeshModel ib_model = ibmeshmodel_force_coupled_init ();

  ib_model.ctx = beam_ptr;
  ib_model.force_ops->node_count = elff_fc_node_count;
  ib_model.force_ops->sync = elff_fc_sync;
  ib_model.force_ops->advance = elff_fc_advance;
  ib_model.force_ops->destroy = NULL;

  return ib_model;
}

/**
 * @brief
 */
IBMeshModel elff_euler_beam_huang_new (double length,
                                       double EI,
                                       double mu,
                                       int nodes,
                                       ib_euler_beam_bcs_t bcs = elff_euler_beam_bcs_default (),
                                       coord s0 = {0}) {
  vertex_t v0 = {s0.x, s0.y, s0.z};
  ib_euler_beam_huang_t beam_ptr =
    ib_euler_beam_huang_new (elff_runtime_get (), v0, bcs, length, EI, mu, nodes);

  IBMeshModel ib_model = ibmeshmodel_force_coupled_init ();

  ib_model.ctx = beam_ptr;
  ib_model.force_ops->node_count = elff_fc_node_count;
  ib_model.force_ops->sync = elff_fc_sync;
  ib_model.force_ops->advance = elff_fc_advance;
  ib_model.force_ops->destroy = NULL;

  return ib_model;
}

/**
 * @brief
 */
IBMeshModel elff_euler_beam_huang_new_theta (double length,
                                             double EI,
                                             double mu,
                                             int nodes,
                                             double theta,
                                             ib_euler_beam_bcs_t bcs = elff_euler_beam_bcs_default (),
                                             coord s0 = {0}) {
  vertex_t v0 = {s0.x, s0.y, s0.z};
  ib_euler_beam_huang_t beam_ptr = ib_euler_beam_huang_new_theta (
    elff_runtime_get (), v0, bcs, length, EI, mu, nodes, theta);

  IBMeshModel ib_model = ibmeshmodel_force_coupled_init ();

  ib_model.ctx = beam_ptr;
  ib_model.force_ops->node_count = elff_fc_node_count;
  ib_model.force_ops->sync = elff_fc_sync;
  ib_model.force_ops->advance = elff_fc_advance;
  ib_model.force_ops->destroy = NULL;

  return ib_model;
}

/**
 * @brief
 */
IBMeshModel elff_euler_beam_huang_new_direction (
  double length,
  double EI,
  double mu,
  int nodes,
  coord direction,
  ib_euler_beam_bcs_t bcs = elff_euler_beam_bcs_default (),
  coord s0 = {0}) {
  vertex_t v0 = {s0.x, s0.y, s0.z};
  vertex_t tangent = {direction.x, direction.y, direction.z};
  ib_euler_beam_huang_t beam_ptr = ib_euler_beam_huang_new_direction (
    elff_runtime_get (), v0, bcs, length, EI, mu, nodes, tangent);

  IBMeshModel ib_model = ibmeshmodel_force_coupled_init ();

  ib_model.ctx = beam_ptr;
  ib_model.force_ops->node_count = elff_fc_node_count;
  ib_model.force_ops->sync = elff_fc_sync;
  ib_model.force_ops->advance = elff_fc_advance;
  ib_model.force_ops->destroy = NULL;

  return ib_model;
}

/**
 * @brief
 */
IBMeshModel elff_capsule_sphere_new (double radius,
                                     coord center = {0},
                                     int refinements = 0,
                                     double elastic_modulus = 1.) {
  ib_capsule_t capsule_ptr = ib_capsule_sphere_new (
    elff_runtime_get (), radius, elff_vertex_from_coord (center), refinements, elastic_modulus);

  return elff_velocity_coupled_model (capsule_ptr);
}

/**
 * @brief
 */
IBMeshModel elff_capsule_ellipsoid_new (coord radii,
                                        coord center = {0},
                                        int refinements = 0,
                                        double elastic_modulus = 1.) {
  ib_capsule_t capsule_ptr = ib_capsule_ellipsoid_new (
    elff_runtime_get (),
    elff_vertex_from_coord (radii),
    elff_vertex_from_coord (center),
    refinements,
    elastic_modulus);

  return elff_velocity_coupled_model (capsule_ptr);
}

/**
 * @brief
 */
IBMeshModel elff_capsule_biconcave_new (double radius,
                                        coord center = {0},
                                        int refinements = 0,
                                        double elastic_modulus = 1.) {
  ib_capsule_t capsule_ptr = ib_capsule_biconcave_new (
    elff_runtime_get (), radius, elff_vertex_from_coord (center), refinements, elastic_modulus);

  return elff_velocity_coupled_model (capsule_ptr);
}

/**
 * @brief
 */
IBMeshModel elff_rigid_body_circle_new (double radius,
                                        int point_count,
                                        double density,
                                        coord center = {0},
                                        coord velocity = {0},
                                        double angle = 0.,
                                        double angular_velocity_z = 0.) {
  vertex_t x = {center.x, center.y, center.z};
  vertex_t v = {velocity.x, velocity.y, velocity.z};
  ib_rigid_body_circle_t body_ptr = ib_rigid_body_circle_new (
    elff_runtime_get (), radius, point_count, density, x, v, angle, angular_velocity_z);

  IBMeshModel ib_model = ibmeshmodel_force_coupled_init ();

  ib_model.ctx = body_ptr;
  ib_model.force_ops->node_count = elff_fc_node_count;
  ib_model.force_ops->sync = elff_fc_sync;
  ib_model.force_ops->advance = elff_fc_advance;
  ib_model.force_ops->destroy = NULL;

  return ib_model;
}

/**
 * @brief
 */
IBMeshModel elff_rigid_body_fibre_new (double length,
                                       double diameter,
                                       int nodes,
                                       double linear_density,
                                       coord center = {0},
                                       coord velocity = {0},
                                       double q_w = 1.,
                                       coord q_xyz = {0},
                                       coord angular_velocity_body = {0}) {
  vertex_t x = {center.x, center.y, center.z};
  vertex_t v = {velocity.x, velocity.y, velocity.z};
  quaternion_t q = {q_w, q_xyz.x, q_xyz.y, q_xyz.z};
  vertex_t omega = {angular_velocity_body.x,
                    angular_velocity_body.y,
                    angular_velocity_body.z};

  ib_rigid_body_fibre_t body_ptr = ib_rigid_body_fibre_new (
    elff_runtime_get (), length, diameter, nodes, linear_density, x, v, q, omega);

  IBMeshModel ib_model = ibmeshmodel_force_coupled_init ();

  ib_model.ctx = body_ptr;
  ib_model.force_ops->node_count = elff_fc_node_count;
  ib_model.force_ops->sync = elff_fc_sync;
  ib_model.force_ops->advance = elff_fc_advance;
  ib_model.force_ops->destroy = NULL;

  return ib_model;
}

/**
 * @brief
 */
IBMeshModel elff_pinned_rigid_body_circle_new (double radius,
                                               int point_count,
                                               double density,
                                               coord center = {0},
                                               double angle = 0.) {
  vertex_t x = {center.x, center.y, center.z};
  ib_pinned_rigid_body_circle_t body_ptr =
    ib_pinned_rigid_body_circle_new (elff_runtime_get (), radius, point_count, density, x, angle);

  IBMeshModel ib_model = ibmeshmodel_force_coupled_init ();

  ib_model.ctx = body_ptr;
  ib_model.force_ops->node_count = elff_fc_node_count;
  ib_model.force_ops->sync = elff_fc_sync;
  ib_model.force_ops->advance = elff_fc_advance;
  ib_model.force_ops->destroy = NULL;

  return ib_model;
}

/**
 * @brief
 */
IBMeshModel elff_pinned_rigid_body_fibre_new (double length,
                                              double diameter,
                                              int nodes,
                                              double linear_density,
                                              coord center = {0},
                                              double q_w = 1.,
                                              coord q_xyz = {0}) {
  vertex_t x = {center.x, center.y, center.z};
  quaternion_t q = {q_w, q_xyz.x, q_xyz.y, q_xyz.z};

  ib_pinned_rigid_body_fibre_t body_ptr = ib_pinned_rigid_body_fibre_new (
    elff_runtime_get (), length, diameter, nodes, linear_density, x, q);

  IBMeshModel ib_model = ibmeshmodel_force_coupled_init ();

  ib_model.ctx = body_ptr;
  ib_model.force_ops->node_count = elff_fc_node_count;
  ib_model.force_ops->sync = elff_fc_sync;
  ib_model.force_ops->advance = elff_fc_advance;
  ib_model.force_ops->destroy = NULL;

  return ib_model;
}

/**
 * @brief
 */
IBMeshModel elff_rigid_body_sphere_new (double radius,
                                        int point_count,
                                        double density,
                                        coord center = {0},
                                        coord velocity = {0},
                                        double q_w = 1.,
                                        coord q_xyz = {0},
                                        coord angular_velocity_body = {0}) {
  vertex_t x = {center.x, center.y, center.z};
  vertex_t v = {velocity.x, velocity.y, velocity.z};
  quaternion_t q = {q_w, q_xyz.x, q_xyz.y, q_xyz.z};
  vertex_t omega = {angular_velocity_body.x,
                    angular_velocity_body.y,
                    angular_velocity_body.z};

  ib_rigid_body_sphere_t body_ptr =
    ib_rigid_body_sphere_new (elff_runtime_get (), radius, point_count, density, x, v, q, omega);

  IBMeshModel ib_model = ibmeshmodel_force_coupled_init ();

  ib_model.ctx = body_ptr;
  ib_model.force_ops->node_count = elff_fc_node_count;
  ib_model.force_ops->sync = elff_fc_sync;
  ib_model.force_ops->advance = elff_fc_advance;
  ib_model.force_ops->destroy = NULL;

  return ib_model;
}

/**
 * @brief
 */
IBMeshModel elff_pinned_rigid_body_sphere_new (double radius,
                                               int point_count,
                                               double density,
                                               coord center = {0},
                                               double q_w = 1.,
                                               coord q_xyz = {0}) {
  vertex_t x = {center.x, center.y, center.z};
  quaternion_t q = {q_w, q_xyz.x, q_xyz.y, q_xyz.z};

  ib_pinned_rigid_body_sphere_t body_ptr =
    ib_pinned_rigid_body_sphere_new (elff_runtime_get (), radius, point_count, density, x, q);

  IBMeshModel ib_model = ibmeshmodel_force_coupled_init ();

  ib_model.ctx = body_ptr;
  ib_model.force_ops->node_count = elff_fc_node_count;
  ib_model.force_ops->sync = elff_fc_sync;
  ib_model.force_ops->advance = elff_fc_advance;
  ib_model.force_ops->destroy = NULL;

  return ib_model;
}
