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

#include "utils.h"
#include "library/ibm/IBMeshManager.h"
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

/** Adapter buffers are indexed by model node id, independent of local storage. */
typedef struct {
  ib_force_coupled_t model; /* Borrowed from IBRuntime. */
  int count;
  vertex_t* forces;
  ib_mesh_t output;
} ELFFForceAdapter;

/** Access the borrowed physical model before transferring the adapter to a manager. */
static inline ib_force_coupled_t elff_force_coupled_handle (IBMeshModel model) {
  assert (model.type == IB_MODEL_FORCE_COUPLED && model.ctx);
  return ((ELFFForceAdapter*) model.ctx)->model;
}

int elff_fc_node_count (void* ctx) {
  return ((ELFFForceAdapter*) ctx)->count;
}

int elff_fc_sync_node (void* ctx, void* mesh, int node_lid, void* raw_node) {
  ELFFForceAdapter* adapter = (ELFFForceAdapter*) ctx;
  IBNode* node = (IBNode*) raw_node;
  if (node_lid < 0 || node_lid >= adapter->output.n)
    return -1;
  foreach_dimension () {
    ibval (npos.x) = adapter->output.position[node_lid].x;
    ibval (nvel.x) = adapter->output.velocity[node_lid].x;
  }
  ibval (nweight) = adapter->output.measure ? adapter->output.measure[node_lid] : 1.;
  return 0;
}

int elff_fc_input_node (void* ctx, void* mesh, int node_lid, void* raw_node) {
  ELFFForceAdapter* adapter = (ELFFForceAdapter*) ctx;
  IBNode* node = (IBNode*) raw_node;
  if (node_lid < 0 || node_lid >= adapter->count)
    return -1;
  foreach_dimension ()
    adapter->forces[node_lid].x = ibval (nforce.x);
  return 0;
}

static void elff_fc_publish_local (ELFFForceAdapter* adapter, IBMesh* mesh) {
  for (size_t ni = 0; ni < mesh->nodes.size; ni++) {
    IBNode* node = mesh->nodes.ptrs[ni];
    int rc = elff_fc_sync_node (adapter, mesh, node->node_lid, node);
    assert (rc == 0);
  }
}

trace void elff_fc_sync (void* ctx, void* raw_mesh) {
  ELFFForceAdapter* adapter = (ELFFForceAdapter*) ctx;
  ib_mesh_free (&adapter->output);
  adapter->output = ib_force_coupled_get_current (adapter->model);
  assert (adapter->output.n == adapter->count);
  elff_fc_publish_local (adapter, (IBMesh*) raw_mesh);
}

trace void elff_fc_advance (void* ctx, void* raw_mesh, double dt) {
  ELFFForceAdapter* adapter = (ELFFForceAdapter*) ctx;
  IBMesh* mesh = (IBMesh*) raw_mesh;
  /* Sparse inputs were already delivered for every logical node by the manager. */
  if (!mesh->sparse_managed)
    for (size_t ni = 0; ni < mesh->nodes.size; ni++) {
      IBNode* node = mesh->nodes.ptrs[ni];
      int rc = elff_fc_input_node (adapter, mesh, node->node_lid, node);
      assert (rc == 0);
    }
  ib_mesh_free (&adapter->output);
  adapter->output = ib_force_coupled_get_next (
    adapter->model, adapter->forces, adapter->count, dt);
  assert (adapter->output.n == adapter->count);
  elff_fc_publish_local (adapter, mesh);
}

static void elff_fc_destroy (void* ctx) {
  ELFFForceAdapter* adapter = (ELFFForceAdapter*) ctx;
  if (!adapter)
    return;
  ib_mesh_free (&adapter->output);
  free (adapter->forces);
  free (adapter);
}

static inline IBMeshModel elff_force_coupled_model (void* ctx) {
  ELFFForceAdapter* adapter = calloc (1, sizeof (ELFFForceAdapter));
  assert (adapter);
  adapter->model = (ib_force_coupled_t) ctx;
  adapter->count = ib_force_coupled_get_number_of_nodes (adapter->model);
  adapter->forces = calloc (adapter->count, sizeof (vertex_t));
  assert (adapter->forces || adapter->count == 0);
  IBMeshModel model = ibmeshmodel_force_coupled_init ();
  model.ctx = adapter;
  model.force_ops->node_count = elff_fc_node_count;
  model.force_ops->sync = elff_fc_sync;
  model.force_ops->sync_node = elff_fc_sync_node;
  model.force_ops->input_node = elff_fc_input_node;
  model.force_ops->advance = elff_fc_advance;
  model.force_ops->destroy = elff_fc_destroy;
  return model;
}

/** Install a non-owning Basilisk adapter for one runtime-owned view. */
static inline IBMeshModel elff_model_from_view (ib_model_view_t view) {
  if (!view.handle) {
    IBMeshModel invalid = {0};
    return invalid;
  }
  if (view.kind == IB_MODEL_VIEW_FORCE_COUPLED)
    return elff_force_coupled_model (view.handle);
  if (view.kind == IB_MODEL_VIEW_VELOCITY_COUPLED)
    return elff_velocity_coupled_model (view.handle);
  IBMeshModel invalid = {0};
  return invalid;
}

/** Register an XML model using its borrowed descriptor handle.
 * Call on every rank, in model order, with the same mesh_gid.
 * The runtime must outlive the mesh manager. Legacy views use fallback_level.
 */
static inline int elff_register_model (int mesh_id, int mesh_gid,
                                       ib_model_view_t view,
                                       int fallback_level) {
  ib_model_info_t info;
  if (ib_model_descriptor_get_info (view.descriptor, &info) != 0)
    return -1;
  IBMeshModel model = elff_model_from_view (view);
  if (info.sparse) {
    int type = view.kind == IB_MODEL_VIEW_FORCE_COUPLED ?
      IB_MODEL_FORCE_COUPLED : IB_MODEL_VELOCITY_COUPLED;
    ibmeshmanager_set_model_sparse_descriptor (mesh_id, mesh_gid, model,
      type, info.node_count, info.owner_rank,
      info.level >= 0 ? info.level : fallback_level);
  }
  else {
    ibmm.meshes[mesh_id].depth = fallback_level;
    ibmeshmanager_set_model (mesh_id, model, info.owner_rank);
  }
  return 0;
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

  return elff_force_coupled_model (beam_ptr);
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

  return elff_force_coupled_model (beam_ptr);
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

  return elff_force_coupled_model (beam_ptr);
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

  return elff_force_coupled_model (beam_ptr);
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

  return elff_force_coupled_model (beam_ptr);
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

  return elff_force_coupled_model (beam_ptr);
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

  return elff_force_coupled_model (beam_ptr);
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

  return elff_force_coupled_model (beam_ptr);
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

  return elff_force_coupled_model (beam_ptr);
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

  return elff_force_coupled_model (beam_ptr);
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

  return elff_force_coupled_model (beam_ptr);
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

  return elff_force_coupled_model (body_ptr);
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

  return elff_force_coupled_model (body_ptr);
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

  return elff_force_coupled_model (body_ptr);
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

  return elff_force_coupled_model (body_ptr);
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

  return elff_force_coupled_model (body_ptr);
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

  return elff_force_coupled_model (body_ptr);
}
