#ifndef ELFF_FLUID_DLM_TYPES_H_
#define ELFF_FLUID_DLM_TYPES_H_

#include <elff/fluid/LevelData.H>
#include <elff/fluid/Options.H>

#include <AMReX_Array.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>

#include <array>
#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace ELFF::Fluid::DLM {

/** @brief One Eulerian stencil index and interpolation weight. */
struct StencilEntry {
  amrex::IntVect cell;
  amrex::Real weight = 0.0;
};

/** @brief Geometric multiplier point with vector work arrays and quadrature. */
struct MarkerPoint {
  std::size_t body_index = 0;
  std::array<amrex::Real, 3> position{0.0, 0.0, 0.0};
  std::array<amrex::Real, 3> center_vector{0.0, 0.0, 0.0};
  std::array<amrex::Real, 3> interpolated_velocity{0.0, 0.0, 0.0};
  std::array<amrex::Real, 3> residual{0.0, 0.0, 0.0};
  std::array<amrex::Real, 3> lambda{0.0, 0.0, 0.0};
  std::array<amrex::Real, 3> search{0.0, 0.0, 0.0};
  std::array<amrex::Real, 3> response{0.0, 0.0, 0.0};
  amrex::Real quadrature_weight = 0.0;
  bool boundary = false;
  std::vector<StencilEntry> stencil;
};

/** @brief Scalar constraint on one staggered velocity component.
 * Boundary markers can produce three constraints; interior points constrain
 * individual face degrees of freedom. Q2 stencil weights may be negative. */
struct ConstraintPoint {
  std::size_t body_index = 0;
  int component = 0;
  amrex::IntVect dof;
  std::array<amrex::Real, 3> position{0.0, 0.0, 0.0};
  std::array<amrex::Real, 3> center_vector{0.0, 0.0, 0.0};
  amrex::Real interpolated_value = 0.0;
  amrex::Real residual = 0.0;
  amrex::Real lambda = 0.0;
  amrex::Real search = 0.0;
  amrex::Real response = 0.0;
  amrex::Real quadrature_weight = 0.0;
  bool boundary = false;
  std::vector<StencilEntry> stencil;
};

/** @brief Rigid-body kinematics, diagonal inertia, and decomposed body loads.
 * Hydrodynamic loads combine raw multiplier and fictitious-fluid corrections.
 */
struct BodyState {
  bool initialized = false;
  std::size_t body_index = 0;
  int id = 1;
  int geometry_id = 0;
  bool particle = true;
  bool prescribed = false;
  amrex::Real fluid_density = 1.0;
  amrex::Real solid_density = 2.0;
  amrex::Real coupling_factor = 1.0;
  amrex::Real mass = 0.0;
  amrex::Real radius = 0.0;
  std::array<amrex::Real, 3> center{0.0, 0.0, 0.0};
  std::array<amrex::Real, 3> inertia{0.0, 0.0, 0.0};
  std::array<amrex::Real, 3> velocity{0.0, 0.0, 0.0};
  std::array<amrex::Real, 3> angular_velocity{0.0, 0.0, 0.0};
  std::array<amrex::Real, 3> previous_velocity{0.0, 0.0, 0.0};
  std::array<amrex::Real, 3> previous_angular_velocity{0.0, 0.0, 0.0};
  std::array<amrex::Real, 3> multiplier_force{0.0, 0.0, 0.0};
  std::array<amrex::Real, 3> multiplier_torque{0.0, 0.0, 0.0};
  std::array<amrex::Real, 3> fictitious_force{0.0, 0.0, 0.0};
  std::array<amrex::Real, 3> fictitious_torque{0.0, 0.0, 0.0};
  std::array<amrex::Real, 3> hydrodynamic_force{0.0, 0.0, 0.0};
  std::array<amrex::Real, 3> hydrodynamic_torque{0.0, 0.0, 0.0};

  void clear() {
    initialized = false;
    body_index = 0;
    id = 1;
    geometry_id = 0;
    particle = true;
    prescribed = false;
    radius = 0.0;
    center = {0.0, 0.0, 0.0};
    mass = 0.0;
    inertia = {0.0, 0.0, 0.0};
    velocity = {0.0, 0.0, 0.0};
    angular_velocity = {0.0, 0.0, 0.0};
    previous_velocity = {0.0, 0.0, 0.0};
    previous_angular_velocity = {0.0, 0.0, 0.0};
    multiplier_force = {0.0, 0.0, 0.0};
    multiplier_torque = {0.0, 0.0, 0.0};
    fictitious_force = {0.0, 0.0, 0.0};
    fictitious_torque = {0.0, 0.0, 0.0};
    hydrodynamic_force = {0.0, 0.0, 0.0};
    hydrodynamic_torque = {0.0, 0.0, 0.0};
  }
};

/** @brief Input geometry and kinematics for one spherical body. */
struct RigidSphere {
  int id = 1;
  int geometry_id = 0;
  bool particle = true;
  bool prescribed = false;
  amrex::Real radius = 0.0;
  amrex::Real density = 0.0;
  amrex::Real mass = 0.0;
  std::array<amrex::Real, 3> inertia{0.0, 0.0, 0.0};
  std::array<amrex::Real, 3> center{0.0, 0.0, 0.0};
  std::array<amrex::Real, 3> velocity{0.0, 0.0, 0.0};
  std::array<amrex::Real, 3> angular_velocity{0.0, 0.0, 0.0};
};

/** @brief Body load identified by id; force/torque include inertia correction.
 */
struct HydrodynamicLoad {
  int id = 1;
  std::array<amrex::Real, 3> multiplier_force{0.0, 0.0, 0.0};
  std::array<amrex::Real, 3> multiplier_torque{0.0, 0.0, 0.0};
  std::array<amrex::Real, 3> fictitious_force{0.0, 0.0, 0.0};
  std::array<amrex::Real, 3> fictitious_torque{0.0, 0.0, 0.0};
  std::array<amrex::Real, 3> force{0.0, 0.0, 0.0};
  std::array<amrex::Real, 3> torque{0.0, 0.0, 0.0};
};

/** @brief Persistent level-local constraint catalog, body state, and solve
 * diagnostics. Warm starts require a compatible constraint catalog;
 * reinitialization clears multiplier history. Algebraic and realized residuals
 * are recorded separately. */
struct State {
  bool initialized = false;
  bool converged = true;
  int iterations_performed = 0;
  amrex::Real initial_residual_norm = 0.0;
  amrex::Real final_residual_norm = 0.0;
  amrex::Real realized_residual_norm = 0.0;
  amrex::Real realized_relative_residual = 0.0;
  amrex::Real realized_residual_linf = 0.0;
  std::string termination_reason = "not_started";
  std::vector<MarkerPoint> markers;
  std::vector<ConstraintPoint> constraints;
  std::vector<BodyState> bodies;
  std::vector<HydrodynamicLoad> hydrodynamic_loads;
  BodyState body;

  void clear() {
    initialized = false;
    converged = true;
    iterations_performed = 0;
    initial_residual_norm = 0.0;
    final_residual_norm = 0.0;
    realized_residual_norm = 0.0;
    realized_relative_residual = 0.0;
    realized_residual_linf = 0.0;
    termination_reason = "not_started";
    markers.clear();
    constraints.clear();
    bodies.clear();
    hydrodynamic_loads.clear();
    body.clear();
  }
};

/** @brief One scalar face-centered MultiFab per spatial direction. */
using FaceFields = amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>;
/** @brief Map a staggered multiplier RHS increment to a velocity increment.
 * The caller provides allocated input/output fields. The selected response
 * supplies mass scaling or a projected Helmholtz solve during DLM iteration. */
using StaggeredFluidResponse = std::function<void(FaceFields const &rhs_delta,
                                                  FaceFields &velocity_delta)>;

} // namespace ELFF::Fluid::DLM

#endif
