#include "dlm/DLMInternal.H"

namespace ELFF::Fluid::DLM {

void
initialize(State& dlm_state, LevelData& state, DLMOptions const& options)
{
    initialize(dlm_state, state, std::vector<RigidSphere>{Detail::optionSphere(options)}, options);
}

void
initialize(State& dlm_state,
           LevelData& state,
           std::vector<RigidSphere> const& spheres,
           DLMOptions const& options)
{
    AMREX_ALWAYS_ASSERT(state.isDefined());
    dlm_state.clear();
    if (options.enabled && options.scheme == "marker") {
        dlm_state.bodies = Detail::buildBodyStates(spheres, options);
        Detail::fillRigidMask(state, dlm_state.bodies, options);
        dlm_state.markers = Detail::buildMarkers(state.geom(), spheres, options);
        dlm_state.constraints =
            Detail::buildConstraints(state.geom(), dlm_state.markers, spheres, options);
        Detail::writeConstraintCatalog(dlm_state.constraints, options.constraint_output);
        dlm_state.initialized = true;
        Detail::syncPrimaryBody(dlm_state);
    } else {
        Detail::fillRigidMask(state, options);
    }
}

void
applyRigidProjection(State& dlm_state,
                     LevelData& state,
                     DLMOptions const& options,
                     amrex::Real dt,
                     StaggeredFluidResponse const& fluid_response)
{
    applyRigidProjection(dlm_state,
                         state,
                         std::vector<RigidSphere>{Detail::optionSphere(options)},
                         options,
                         dt,
                         fluid_response);
}

void
applyRigidProjection(State& dlm_state,
                     LevelData& state,
                     std::vector<RigidSphere> const& spheres,
                     DLMOptions const& options,
                     amrex::Real dt,
                     StaggeredFluidResponse const& fluid_response)
{
    AMREX_ALWAYS_ASSERT(state.isDefined());
    if (!options.enabled) {
        state.rigidVelocity().setVal(0.0);
        state.dlmForce().setVal(0.0);
        dlm_state.clear();
        return;
    }

    if (options.scheme == "cell_mask") {
        if (spheres.size() != 1) {
            amrex::Abort("dlm.scheme = cell_mask supports only one sphere");
        }
        Detail::fillRigidMask(state, options);
        dlm_state.clear();
        Detail::applyCellMaskProjection(state, options, dt);
    } else if (options.scheme == "marker") {
        Detail::ensureMarkerState(dlm_state, state, spheres, options);
        Detail::fillRigidMask(state, dlm_state.bodies, options);
        Detail::applyMarkerProjection(dlm_state, state, spheres, options, dt, fluid_response);
    } else {
        amrex::Abort("dlm.scheme must be either 'marker' or 'cell_mask'");
    }
}

void
initialize(LevelData& state, DLMOptions const& options)
{
    State dlm_state;
    initialize(dlm_state, state, options);
}

void
applyRigidProjection(LevelData& state, DLMOptions const& options, amrex::Real dt)
{
    State dlm_state;
    initialize(dlm_state, state, options);
    StaggeredFluidResponse identity_response = [] (FaceFields const& rhs_delta,
                                                   FaceFields& velocity_delta) {
        for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
            amrex::MultiFab::Copy(velocity_delta[comp], rhs_delta[comp], 0, 0, 1, 0);
        }
    };
    applyRigidProjection(dlm_state, state, options, dt, identity_response);
}

void
enforceInteriorVelocity(State& dlm_state,
                        LevelData& state,
                        DLMOptions const& options,
                        amrex::Real dt)
{
    AMREX_ALWAYS_ASSERT(state.isDefined());
    if (!options.enabled ||
        options.scheme != "marker" ||
        !options.enforce_interior_velocity ||
        dlm_state.bodies.empty())
    {
        return;
    }

    amrex::MultiFab previous_force(state.dlmForce().boxArray(),
                                   state.dlmForce().DistributionMap(),
                                   AMREX_SPACEDIM,
                                   state.dlmForce().nGrowVect(),
                                   amrex::MFInfo(),
                                   state.dlmForce().Factory());
    amrex::MultiFab::Copy(previous_force,
                          state.dlmForce(),
                          0,
                          0,
                          AMREX_SPACEDIM,
                          state.dlmForce().nGrowVect());
    Detail::fillRigidMask(state, dlm_state.bodies, options);
    amrex::MultiFab::Copy(state.dlmForce(),
                          previous_force,
                          0,
                          0,
                          AMREX_SPACEDIM,
                          state.dlmForce().nGrowVect());
    amrex::MultiFab cumulative_delta(state.velocity().boxArray(),
                                     state.velocity().DistributionMap(),
                                     AMREX_SPACEDIM,
                                     state.velocity().nGrowVect(),
                                     amrex::MFInfo(),
                                     state.velocity().Factory());
    cumulative_delta.setVal(0.0);
    Detail::enforceInteriorRigidVelocityCells(state, dlm_state.bodies, cumulative_delta);
    Detail::fillRigidVelocityDiagnostic(state, dlm_state.bodies, options);
    Detail::addDlmForceFromDelta(state, cumulative_delta, dt);
}

} // namespace ELFF::Fluid::DLM
