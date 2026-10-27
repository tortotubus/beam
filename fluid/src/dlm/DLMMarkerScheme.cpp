#include "DLMInternal.H"

#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <limits>

namespace ELFF::Fluid::DLM::Detail {

namespace {

bool
sameConstraintMaterialPoint(ConstraintPoint const& lhs, ConstraintPoint const& rhs)
{
    if (lhs.body_index != rhs.body_index ||
        lhs.component != rhs.component ||
        lhs.boundary != rhs.boundary)
    {
        return false;
    }
    const amrex::Real scale = std::max({amrex::Real(1.0),
                                        std::abs(lhs.quadrature_weight),
                                        std::abs(rhs.quadrature_weight)});
    if (std::abs(lhs.quadrature_weight - rhs.quadrature_weight)
        > amrex::Real(1.0e-12) * scale) {
        return false;
    }
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const amrex::Real coord_scale = std::max({amrex::Real(1.0),
                                                  std::abs(lhs.position[dir]),
                                                  std::abs(rhs.position[dir]),
                                                  std::abs(lhs.center_vector[dir]),
                                                  std::abs(rhs.center_vector[dir])});
        if (std::abs(lhs.position[dir] - rhs.position[dir])
            > amrex::Real(1.0e-12) * coord_scale) {
            return false;
        }
        if (std::abs(lhs.center_vector[dir] - rhs.center_vector[dir])
            > amrex::Real(1.0e-12) * coord_scale) {
            return false;
        }
    }
    if (lhs.stencil.size() != rhs.stencil.size()) {
        return false;
    }
    for (std::size_t i = 0; i < lhs.stencil.size(); ++i) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            if (lhs.stencil[i].cell[dir] != rhs.stencil[i].cell[dir]) {
                return false;
            }
        }
        const amrex::Real weight_scale = std::max({amrex::Real(1.0),
                                                   std::abs(lhs.stencil[i].weight),
                                                   std::abs(rhs.stencil[i].weight)});
        if (std::abs(lhs.stencil[i].weight - rhs.stencil[i].weight)
            > amrex::Real(1.0e-12) * weight_scale) {
            return false;
        }
    }
    return true;
}

void
copyCompatibleLambda(std::vector<ConstraintPoint> const& old_constraints,
                     std::vector<ConstraintPoint>& new_constraints)
{
    if (old_constraints.size() != new_constraints.size()) {
        return;
    }
    for (std::size_t i = 0; i < new_constraints.size(); ++i) {
        if (!sameConstraintMaterialPoint(old_constraints[i], new_constraints[i])) {
            return;
        }
    }
    for (std::size_t i = 0; i < new_constraints.size(); ++i) {
        new_constraints[i].lambda = old_constraints[i].lambda;
    }
}

} // namespace

void
ensureMarkerState(State& dlm_state, LevelData const& state, DLMOptions const& options)
{
    ensureMarkerState(dlm_state, state, std::vector<RigidSphere>{optionSphere(options)}, options);
}

void
ensureMarkerState(State& dlm_state,
                  LevelData const& state,
                  std::vector<RigidSphere> const& spheres,
                  DLMOptions const& options)
{
    std::vector<ConstraintPoint> old_constraints = std::move(dlm_state.constraints);
    std::vector<BodyState> old_bodies = std::move(dlm_state.bodies);
    dlm_state.bodies = buildBodyStates(spheres, options);
    for (BodyState& body : dlm_state.bodies) {
        for (BodyState const& old_body : old_bodies) {
            if (body.id == old_body.id) {
                body.previous_velocity = old_body.velocity;
                body.previous_angular_velocity = old_body.angular_velocity;
                break;
            }
        }
    }
    dlm_state.markers = buildMarkers(state.geom(), spheres, options);
    dlm_state.constraints = buildConstraints(state.geom(), dlm_state.markers, spheres, options);
    copyCompatibleLambda(old_constraints, dlm_state.constraints);
    dlm_state.initialized = true;
    syncPrimaryBody(dlm_state);
}

void
applyMarkerProjection(State& dlm_state,
                      LevelData& state,
                      DLMOptions const& options,
                      amrex::Real dt,
                      StaggeredFluidResponse const& fluid_response)
{
    applyMarkerProjection(dlm_state,
                          state,
                          std::vector<RigidSphere>{optionSphere(options)},
                          options,
                          dt,
                          fluid_response);
}

void
applyMarkerProjection(State& dlm_state,
                      LevelData& state,
                      std::vector<RigidSphere> const& spheres,
                      DLMOptions const& options,
                      amrex::Real dt,
                      StaggeredFluidResponse const& fluid_response)
{
    BL_PROFILE("ELFF::Fluid::DLM::applyMarkerProjection");
    const amrex::Real profile_start = amrex::ParallelDescriptor::second();
    amrex::Real phase_start = profile_start;
    amrex::Real setup_time = 0.0;
    amrex::Real cell_to_face_time = 0.0;
    amrex::Real initial_body_time = 0.0;
    amrex::Real warm_start_time = 0.0;
    amrex::Real residual_time = 0.0;
    amrex::Real search_body_time = 0.0;
    amrex::Real spread_time = 0.0;
    amrex::Real fluid_response_time = 0.0;
    amrex::Real response_interp_time = 0.0;
    amrex::Real update_time = 0.0;
    amrex::Real final_diagnostics_time = 0.0;

    ensureMarkerState(dlm_state, state, spheres, options);
    std::vector<ConstraintPoint>& constraints = dlm_state.constraints;
    resetConstraintUzawaFields(constraints);
    const bool staggered_primary =
        options.multiplier_layout == "pacific" && options.fluid_response == "mass";

    FaceFields face_values;
    FaceFields face_rhs;
    FaceFields face_delta;
    FaceFields realized_face_delta;
    defineFaceFieldsImpl(state, face_values);
    defineFaceFieldsImpl(state, face_rhs);
    defineFaceFieldsImpl(state, face_delta);
    defineFaceFieldsImpl(state, realized_face_delta);
    FaceOwnerMasks face_owner_masks = makeFaceOwnerMasks(face_values, state.geom());
    ConstraintStencilOperator constraint_operator(face_values, constraints);

    amrex::MultiFab cumulative_delta(state.velocity().boxArray(),
                                     state.velocity().DistributionMap(),
                                     AMREX_SPACEDIM,
                                     state.velocity().nGrowVect(),
                                     amrex::MFInfo(),
                                     state.velocity().Factory());
    cumulative_delta.setVal(0.0);

    amrex::MultiFab cell_delta(state.velocity().boxArray(),
                               state.velocity().DistributionMap(),
                               AMREX_SPACEDIM,
                               state.velocity().nGrowVect(),
                               amrex::MFInfo(),
                               state.velocity().Factory());
    cell_delta.setVal(0.0);
    setup_time = amrex::ParallelDescriptor::second() - phase_start;

    phase_start = amrex::ParallelDescriptor::second();
    if (staggered_primary) {
        for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
            amrex::MultiFab::Copy(face_values[comp],
                                  state.faceVelocity()[comp],
                                  0,
                                  0,
                                  1,
                                  0);
        }
    } else {
        cellVectorToComponentFacesImpl(state, state.velocity(), face_values);
    }
    cell_to_face_time += amrex::ParallelDescriptor::second() - phase_start;

    phase_start = amrex::ParallelDescriptor::second();
    std::vector<RigidMotion> motions(dlm_state.bodies.size());
    for (std::size_t body_index = 0; body_index < dlm_state.bodies.size(); ++body_index) {
        motions[body_index] =
            solveInitialBodyMotion(dlm_state.bodies[body_index], options, constraints, dt);
        setBodyMotion(dlm_state.bodies[body_index], motions[body_index]);
    }
    syncPrimaryBody(dlm_state);
    initial_body_time = amrex::ParallelDescriptor::second() - phase_start;

    phase_start = amrex::ParallelDescriptor::second();
    if (hasNonzeroConstraintField(constraints, MarkerField::Lambda)) {
        constraint_operator.spread(
            state, constraints, MarkerField::Lambda, -1.0, face_rhs);
        fluid_response(face_rhs, face_delta);
        if (staggered_primary) {
            for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
                amrex::MultiFab::Add(face_values[comp], face_delta[comp], 0, 0, 1, 0);
            }
        } else {
            addFaceDeltaToVelocity(state, face_delta, cell_delta, cumulative_delta, 1.0);
            state.velocity().FillBoundary(state.geom().periodicity());
            cellVectorToComponentFacesImpl(state, state.velocity(), face_values);
        }
    }
    warm_start_time = amrex::ParallelDescriptor::second() - phase_start;

    phase_start = amrex::ParallelDescriptor::second();
    constraint_operator.interpolate(face_values, face_owner_masks, constraints, false);
    fillConstraintResiduals(constraints, dlm_state.bodies);
    copyConstraintField(constraints, MarkerField::Search, MarkerField::Residual);
    residual_time = amrex::ParallelDescriptor::second() - phase_start;

    const int niter = std::max(0, options.constraint_iterations);
    const amrex::Real relaxation = std::clamp(options.relaxation,
                                             amrex::Real(0.0),
                                             amrex::Real(1.0));
    amrex::Real rdotr = constraintDot(constraints, MarkerField::Residual,
                                      MarkerField::Residual);
    const amrex::Real initial_rdotr = std::max(rdotr, amrex::Real(0.0));
    const amrex::Real initial_residual = std::sqrt(initial_rdotr);
    const amrex::Real abs_tol = std::max(options.constraint_abs_tol, amrex::Real(0.0));
    const amrex::Real rel_tol = std::max(options.constraint_rel_tol, amrex::Real(0.0));
    const amrex::Real tolerance =
        std::max(abs_tol, rel_tol * initial_residual);
    const amrex::Real tolerance_sq =
        std::max(tiny, tolerance * tolerance);
    const amrex::Real effective_tolerance = std::sqrt(tolerance_sq);

    dlm_state.converged = (rdotr <= tolerance_sq);
    dlm_state.termination_reason = dlm_state.converged
        ? "initial_residual"
        : "iteration_limit";
    dlm_state.iterations_performed = 0;
    dlm_state.initial_residual_norm = initial_residual;
    dlm_state.final_residual_norm = initial_residual;

    for (int iter = 0; iter < niter && !dlm_state.converged; ++iter) {
        BL_PROFILE("ELFF::Fluid::DLM::UzawaIteration");
        if (rdotr <= tolerance_sq) {
            dlm_state.converged = true;
            break;
        }
        dlm_state.iterations_performed = iter + 1;

        phase_start = amrex::ParallelDescriptor::second();
        std::vector<RigidMotion> body_responses(dlm_state.bodies.size());
        for (std::size_t body_index = 0; body_index < dlm_state.bodies.size(); ++body_index) {
            body_responses[body_index] =
                solveSearchBodyMotion(dlm_state.bodies[body_index], options, constraints, dt);
        }
        search_body_time += amrex::ParallelDescriptor::second() - phase_start;

        phase_start = amrex::ParallelDescriptor::second();
        constraint_operator.spread(
            state, constraints, MarkerField::Search, 1.0, face_rhs);
        spread_time += amrex::ParallelDescriptor::second() - phase_start;
        phase_start = amrex::ParallelDescriptor::second();
        fluid_response(face_rhs, face_delta);
        fluid_response_time += amrex::ParallelDescriptor::second() - phase_start;
        phase_start = amrex::ParallelDescriptor::second();
        if (staggered_primary) {
            constraint_operator.interpolate(
                face_delta, face_owner_masks, constraints, true);
        } else {
            componentFacesToCellVectorImpl(state, face_delta, cell_delta);
            cellVectorToComponentFacesImpl(state, cell_delta, realized_face_delta);
            constraint_operator.interpolate(
                realized_face_delta, face_owner_masks, constraints, true);
        }
        subtractRigidResponseFromConstraints(constraints, body_responses, dlm_state.bodies);
        response_interp_time += amrex::ParallelDescriptor::second() - phase_start;

        phase_start = amrex::ParallelDescriptor::second();
        const amrex::Real denom = constraintDot(constraints, MarkerField::Search,
                                                MarkerField::Response);
        if (!std::isfinite(denom)) {
            amrex::Abort("DLM marker Uzawa denominator is non-finite");
        }
        if (std::abs(denom) <= tiny) {
            dlm_state.termination_reason = "denominator_breakdown";
            break;
        }

        const amrex::Real alpha = relaxation * rdotr / denom;
        if (!std::isfinite(alpha)) {
            amrex::Abort("DLM marker Uzawa alpha is non-finite");
        }
        updateLambdaAndResidual(constraints, alpha);
        if (staggered_primary) {
            for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
                amrex::MultiFab::Saxpy(
                    face_values[comp], alpha, face_delta[comp], 0, 0, 1, 0);
            }
        } else {
            addDeltaToVelocity(state, cell_delta, cumulative_delta, alpha);
        }
        for (std::size_t body_index = 0; body_index < dlm_state.bodies.size(); ++body_index) {
            addBodyMotion(dlm_state.bodies[body_index], body_responses[body_index], alpha);
        }
        syncPrimaryBody(dlm_state);
        if (!staggered_primary) {
            state.velocity().FillBoundary(state.geom().periodicity());
        }

        const amrex::Real new_rdotr = constraintDot(constraints, MarkerField::Residual,
                                                    MarkerField::Residual);
        if (new_rdotr <= tolerance_sq) {
            rdotr = new_rdotr;
            dlm_state.converged = true;
            dlm_state.termination_reason = "converged";
            break;
        }

        const amrex::Real beta = new_rdotr / rdotr;
        updateSearch(constraints, beta);
        rdotr = new_rdotr;
        update_time += amrex::ParallelDescriptor::second() - phase_start;
    }
    dlm_state.final_residual_norm = std::sqrt(std::max(rdotr, amrex::Real(0.0)));

    phase_start = amrex::ParallelDescriptor::second();
    if (staggered_primary) {
        componentFacesToCellVectorImpl(state, face_values, cell_delta);
        amrex::MultiFab::Copy(cumulative_delta,
                              cell_delta,
                              0,
                              0,
                              AMREX_SPACEDIM,
                              0);
        amrex::MultiFab::Subtract(cumulative_delta,
                                  state.velocity(),
                                  0,
                                  0,
                                  AMREX_SPACEDIM,
                                  0);
        amrex::MultiFab::Copy(state.velocity(),
                              cell_delta,
                              0,
                              0,
                              AMREX_SPACEDIM,
                              0);
        state.velocity().FillBoundary(state.geom().periodicity());
    } else {
        cellVectorToComponentFacesImpl(state, state.velocity(), face_values);
    }
    constraint_operator.interpolate(face_values, face_owner_masks, constraints, false);
    if (options.enforce_interior_velocity) {
        enforceInteriorRigidVelocityCells(state, dlm_state.bodies, cumulative_delta);
        cellVectorToComponentFacesImpl(state, state.velocity(), face_values);
        constraint_operator.interpolate(face_values, face_owner_masks, constraints, false);
    }
    fillConstraintResiduals(constraints, dlm_state.bodies);
    const amrex::Real realized_rdotr =
        std::max(amrex::Real(0.0),
                 constraintDot(constraints, MarkerField::Residual, MarkerField::Residual));
    dlm_state.realized_residual_norm = std::sqrt(realized_rdotr);
    dlm_state.realized_relative_residual = initial_residual > tiny
        ? dlm_state.realized_residual_norm / initial_residual
        : dlm_state.realized_residual_norm;
    dlm_state.realized_residual_linf = 0.0;
    for (ConstraintPoint const& constraint : constraints) {
        dlm_state.realized_residual_linf =
            std::max(dlm_state.realized_residual_linf, std::abs(constraint.residual));
    }
    // Reconstructing J*u from the realized float field is not bitwise
    // identical to the algebraic residual update.  Bound that roundoff in
    // units of amrex::Real and of the constrained velocity signal, rather
    // than introducing a backend-specific DLM tolerance.
    amrex::Real constrained_signal_sq = amrex::Real(0.0);
    for (ConstraintPoint const& constraint : constraints) {
        const amrex::Real rigid_value =
            constraint.interpolated_value + constraint.residual;
        const amrex::Real signal =
            std::max(std::abs(rigid_value),
                     std::abs(constraint.interpolated_value));
        constrained_signal_sq += std::abs(constraint.quadrature_weight)
                               * signal * signal;
    }
    const amrex::Real realized_roundoff_allowance =
        amrex::Real(2.0) * std::numeric_limits<amrex::Real>::epsilon()
        * std::sqrt(std::max(constrained_signal_sq, amrex::Real(0.0)));
    const amrex::Real realized_tolerance =
        effective_tolerance + realized_roundoff_allowance;
    const bool realized_converged =
        dlm_state.realized_residual_norm <= realized_tolerance;
    if (dlm_state.converged && !realized_converged) {
        dlm_state.termination_reason = "realized_residual_mismatch";
    } else if (realized_converged) {
        dlm_state.termination_reason = "converged";
    }
    dlm_state.converged = realized_converged;
    for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
        amrex::MultiFab::Copy(state.faceVelocity()[comp],
                              face_values[comp],
                              0,
                              0,
                              1,
                              0);
    }
    fillRigidVelocityDiagnostic(state, dlm_state.bodies, options);
    setDlmForceFromDelta(state, cumulative_delta, dt);
    setHydrodynamicLoads(dlm_state, dt);
    final_diagnostics_time = amrex::ParallelDescriptor::second() - phase_start;

    if (options.profile_phases) {
        amrex::Print() << "ELFF_FLUID_PROFILE dlm_marker"
                       << " total "
                       << amrex::ParallelDescriptor::second() - profile_start
                       << " bodies " << dlm_state.bodies.size()
                       << " markers " << dlm_state.markers.size()
                       << " constraints " << constraints.size()
                       << " iterations " << dlm_state.iterations_performed
                       << " algebraic_initial " << dlm_state.initial_residual_norm
                       << " algebraic_final " << dlm_state.final_residual_norm
                       << " realized_l2 " << dlm_state.realized_residual_norm
                       << " realized_tolerance " << realized_tolerance
                       << " realized_relative " << dlm_state.realized_relative_residual
                       << " realized_linf " << dlm_state.realized_residual_linf
                       << " termination " << dlm_state.termination_reason
                       << " setup " << setup_time
                       << " cell_to_face " << cell_to_face_time
                       << " initial_body " << initial_body_time
                       << " warm_start " << warm_start_time
                       << " residual " << residual_time
                       << " search_body " << search_body_time
                       << " spread " << spread_time
                       << " fluid_response " << fluid_response_time
                       << " response_interp " << response_interp_time
                       << " update " << update_time
                       << " final_diagnostics " << final_diagnostics_time
                       << "\n";
    }

    if (!dlm_state.converged && options.fail_on_nonconvergence) {
        amrex::Abort("DLM marker Uzawa iteration failed to converge");
    }
}

} // namespace ELFF::Fluid::DLM::Detail
