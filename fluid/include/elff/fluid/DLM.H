#ifndef ELFF_FLUID_DLM_H_
#define ELFF_FLUID_DLM_H_

#include <elff/fluid/DLMTypes.H>

#include <AMReX_BoxArray.H>
#include <AMReX_TagBox.H>

#include <string>

namespace ELFF::Fluid::DLM {

/** @brief Build a single-body constraint state from DLMOptions.
 * @pre The fluid level is defined. Existing multiplier history is cleared. */
void initialize(State &dlm_state, LevelData &state, DLMOptions const &options);

/** @brief Build level-local constraints for the supplied spherical bodies.
 * Rebuilds body state, rigid mask, markers, and staggered constraints. */
void initialize(State &dlm_state, LevelData &state,
                std::vector<RigidSphere> const &spheres,
                DLMOptions const &options);

/** @brief Enforce the option-defined body motion using the supplied fluid
 * response. Updates fluid fields, body kinematics, loads, and convergence
 * diagnostics.
 * @param dlm_state Persistent constraints, body state, and diagnostics.
 * @param state Fluid fields updated by the solve.
 * @param options Constraint scheme, body properties, and convergence policy.
 * @param dt Timestep used by the constraint and body response.
 * @param fluid_response Operator mapping multiplier RHS to velocity response.
 */
void applyRigidProjection(State &dlm_state, LevelData &state,
                          DLMOptions const &options, amrex::Real dt,
                          StaggeredFluidResponse const &fluid_response);

/** @brief Apply a level-local constraint solve for multiple bodies.
 * This is not a composite coarse-fine constraint solve. Nonconvergence aborts
 * only when requested by DLMOptions::fail_on_nonconvergence. */
void applyRigidProjection(State &dlm_state, LevelData &state,
                          std::vector<RigidSphere> const &spheres,
                          DLMOptions const &options, amrex::Real dt,
                          StaggeredFluidResponse const &fluid_response);

/// Initialize diagnostics using temporary, nonpersistent DLM state.
void initialize(LevelData &state, DLMOptions const &options);
/// Convenience solve with temporary state and an identity fluid response.
void applyRigidProjection(LevelData &state, DLMOptions const &options,
                          amrex::Real dt);

/** @brief Optionally overwrite interior velocity with the rigid target.
 * This diagnostic/algorithm option differs from enforcing marker constraints
 * and is disabled in the matched Pacific reference case. */
void enforceInteriorVelocity(State &dlm_state, LevelData &state,
                             DLMOptions const &options, amrex::Real dt);

/// Add body-volume, rigid-mask, and stencil-support refinement tags as
/// configured.
void tagRefinement(amrex::TagBoxArray &tags, LevelData const &state,
                   std::vector<RigidSphere> const &spheres,
                   DLMOptions const &options);

/** @brief Check that this level covers every required constraint stencil.
 * @param state Fluid level whose box layout is checked.
 * @param spheres Body geometry used to construct constraint support.
 * @param options Interpolation and constraint-layout settings.
 * @param reason Optional output describing the first coverage failure.
 * @return Whether the level layout contains complete required support. */
[[nodiscard]] bool
hasFullStencilSupport(LevelData const &state,
                      std::vector<RigidSphere> const &spheres,
                      DLMOptions const &options, std::string *reason = nullptr);

/// Abort with a coverage diagnostic if hasFullStencilSupport() is false.
void requireFullStencilSupport(LevelData const &state,
                               std::vector<RigidSphere> const &spheres,
                               DLMOptions const &options);

} // namespace ELFF::Fluid::DLM

#endif
