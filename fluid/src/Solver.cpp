#include "SolverInternals.H"

#include <elff/fluid/Boundary.H>
#include <elff/fluid/DLM.H>
#include <elff/fluid/Solver.H>

#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

namespace ELFF::Fluid {

using namespace SolverDetail;

Solver::Solver(SolverOptions options)
    : options_(options)
{
    standalone_dlm_spheres_ = {optionSphere(options_.dlm)};
    options_.print();
}

void
Solver::initData()
{
    BL_PROFILE("ELFF::Fluid::Solver::initData");
    const amrex::Real profile_start = amrex::ParallelDescriptor::second();
    InitFromScratch(time_);
    defineFluxRegisters();
    defineVelocitySyncCorrections();
    defineSyncProjectionResiduals();
    synchronizeHierarchy();
    if (options_.profile_phases) {
        amrex::Print() << "ELFF_FLUID_PROFILE init total "
                       << amrex::ParallelDescriptor::second() - profile_start
                       << "\n";
    }
}

void
Solver::evolve()
{
    BL_PROFILE("ELFF::Fluid::Solver::evolve");
    if (options_.plot_initial && options_.plot_int > 0) {
        writePlotFile();
    }

    while (step_ < options_.max_step &&
           (options_.stop_time < 0.0 || time_ < options_.stop_time))
    {
        BL_PROFILE("ELFF::Fluid::Solver::evolve::step");
        const amrex::Real dt = options_.stop_time > 0.0
            ? std::min(options_.fixed_dt, options_.stop_time - time_)
            : options_.fixed_dt;

        const amrex::Real step_profile_start = amrex::ParallelDescriptor::second();
        amrex::Real phase_start = step_profile_start;
        amrex::Real setup_time = 0.0;
        amrex::Real explicit_time = 0.0;
        amrex::Real viscosity_time = 0.0;
        amrex::Real constraint_time = 0.0;
        amrex::Real persist_time = 0.0;
        amrex::Real reflux_time = 0.0;
        amrex::Real sync_time = 0.0;
        amrex::Real regrid_time = 0.0;
        amrex::Real dlm_reassert_time = 0.0;
        amrex::Real plot_time = 0.0;

        snapshotOldStates(time_);
        resetFluxRegisters();
        resetVelocitySyncCorrections();
        resetSyncProjectionResiduals();
        setup_time = amrex::ParallelDescriptor::second() - phase_start;

        phase_start = amrex::ParallelDescriptor::second();
        for (int lev = 0; lev <= finestLevel(); ++lev) {
            advanceLevelExplicit(lev, time_, dt);
        }
        explicit_time = amrex::ParallelDescriptor::second() - phase_start;

        phase_start = amrex::ParallelDescriptor::second();
        applyViscosityHierarchy(dt);
        viscosity_time = amrex::ParallelDescriptor::second() - phase_start;

        const amrex::Real new_time = time_ + dt;
        for (int lev = 0; lev <= finestLevel(); ++lev) {
            if (lev < static_cast<int>(levels_.size()) && levels_[lev] != nullptr) {
                levels_[lev]->setNewTime(new_time);
            }
        }

        const bool dlm_enabled_this_step = dlmOptionsForCurrentStep().enabled;
        if (dlm_enabled_this_step) {
            phase_start = amrex::ParallelDescriptor::second();
            refluxHierarchy();
            reflux_time = amrex::ParallelDescriptor::second() - phase_start;
            phase_start = amrex::ParallelDescriptor::second();
            const DLMOptions dlm_options = dlmOptionsForCurrentStep();
            const bool retain_single_level_staggered_velocity =
                finestLevel() == 0 &&
                dlm_options.multiplier_layout == "pacific" &&
                dlm_options.fluid_response == "mass";
            if (!retain_single_level_staggered_velocity) {
                synchronizeHierarchy();
            }
            sync_time = amrex::ParallelDescriptor::second() - phase_start;
        }

        phase_start = amrex::ParallelDescriptor::second();
        for (int lev = 0; lev <= finestLevel(); ++lev) {
            projectAndApplyConstraints(lev, dt);
        }
        constraint_time = amrex::ParallelDescriptor::second() - phase_start;
        writeDlmHistory(new_time);

        phase_start = amrex::ParallelDescriptor::second();
        persistStandaloneDlmSpheres(dt);
        persist_time = amrex::ParallelDescriptor::second() - phase_start;

        if (!dlm_enabled_this_step) {
            phase_start = amrex::ParallelDescriptor::second();
            refluxHierarchy();
            reflux_time = amrex::ParallelDescriptor::second() - phase_start;
            phase_start = amrex::ParallelDescriptor::second();
            synchronizeHierarchy();
            sync_time = amrex::ParallelDescriptor::second() - phase_start;
        }
        phase_start = amrex::ParallelDescriptor::second();
        enforceDlmInteriorVelocity(dt);
        dlm_reassert_time += amrex::ParallelDescriptor::second() - phase_start;

        time_ += dt;
        ++step_;

        phase_start = amrex::ParallelDescriptor::second();
        regridIfNeeded();
        regrid_time = amrex::ParallelDescriptor::second() - phase_start;
        phase_start = amrex::ParallelDescriptor::second();
        enforceDlmInteriorVelocity(dt);
        dlm_reassert_time += amrex::ParallelDescriptor::second() - phase_start;

        if (options_.verbose > 0) {
            amrex::Print() << "Advanced step " << step_
                           << " to time " << time_ << "\n";
        }

        phase_start = amrex::ParallelDescriptor::second();
        if (options_.plot_int > 0 && step_ % options_.plot_int == 0) {
            writePlotFile();
        }
        plot_time = amrex::ParallelDescriptor::second() - phase_start;

        if (options_.profile_phases) {
            const amrex::Real total_time =
                amrex::ParallelDescriptor::second() - step_profile_start;
            amrex::Print() << "ELFF_FLUID_PROFILE step " << step_
                           << " total " << total_time
                           << " setup " << setup_time
                           << " explicit " << explicit_time
                           << " viscosity " << viscosity_time
                           << " project_dlm " << constraint_time
                           << " persist " << persist_time
                           << " reflux " << reflux_time
                           << " sync " << sync_time
                           << " regrid " << regrid_time
                           << " dlm_reassert " << dlm_reassert_time
                           << " plot " << plot_time
                           << "\n";
        }
    }
}

void
Solver::MakeNewLevelFromScratch(int lev,
                                amrex::Real time,
                                amrex::BoxArray const& ba,
                                amrex::DistributionMapping const& dm)
{
    ensureLevelCapacity(lev);
    levels_[lev] = std::make_unique<LevelData>();
    levels_[lev]->define(lev, Geom(lev), ba, dm);
    levels_[lev]->setInitialState(time);
    if (lev > 0) {
        fillLevelFromCoarse(lev);
    }
    Boundary::fillVelocityPhysicalBoundary(*levels_[lev], options_);
    refreshFaceVelocityFromCell(*levels_[lev], options_);
    const DLMOptions dlm_options = dlmOptionsForCurrentStep();
    DLM::initialize(dlm_states_[lev],
                    *levels_[lev],
                    currentDlmSpheres(),
                    dlm_options);
}

void
Solver::MakeNewLevelFromCoarse(int lev,
                               amrex::Real time,
                               amrex::BoxArray const& ba,
                               amrex::DistributionMapping const& dm)
{
    ensureLevelCapacity(lev);
    levels_[lev] = std::make_unique<LevelData>();
    levels_[lev]->define(lev, Geom(lev), ba, dm);
    levels_[lev]->setInitialState(time);
    fillLevelFromCoarse(lev);
    Boundary::fillVelocityPhysicalBoundary(*levels_[lev], options_);
    refreshFaceVelocityFromCell(*levels_[lev], options_);
    const DLMOptions dlm_options = dlmOptionsForCurrentStep();
    DLM::initialize(dlm_states_[lev],
                    *levels_[lev],
                    currentDlmSpheres(),
                    dlm_options);
}

void
Solver::RemakeLevel(int lev,
                    amrex::Real time,
                    amrex::BoxArray const& ba,
                    amrex::DistributionMapping const& dm)
{
    ensureLevelCapacity(lev);
    auto old_level = std::move(levels_[lev]);

    levels_[lev] = std::make_unique<LevelData>();
    levels_[lev]->define(lev, Geom(lev), ba, dm);
    levels_[lev]->setInitialState(time);
    if (lev > 0) {
        fillLevelFromCoarse(lev);
    }

    if (old_level != nullptr && old_level->isDefined()) {
        levels_[lev]->copyStateFrom(*old_level);
    }
    Boundary::fillVelocityPhysicalBoundary(*levels_[lev], options_);
    refreshFaceVelocityFromCell(*levels_[lev], options_);
    const DLMOptions dlm_options = dlmOptionsForCurrentStep();
    DLM::initialize(dlm_states_[lev],
                    *levels_[lev],
                    currentDlmSpheres(),
                    dlm_options);
}

void
Solver::ClearLevel(int lev)
{
    if (lev >= 0 && lev < static_cast<int>(levels_.size())) {
        levels_[lev].reset();
        dlm_states_[lev].clear();
    }
    if (lev >= 0 && lev < static_cast<int>(adv_flux_registers_.size())) {
        adv_flux_registers_[lev].reset();
    }
    if (lev >= 0 && lev < static_cast<int>(visc_flux_registers_.size())) {
        visc_flux_registers_[lev].reset();
    }
    if (lev >= 0 && lev < static_cast<int>(velocity_sync_corrections_.size())) {
        velocity_sync_corrections_[lev].reset();
    }
    if (lev >= 0 && lev < static_cast<int>(sync_projection_residuals_.size())) {
        sync_projection_residuals_[lev].reset();
    }
}

void
Solver::ErrorEst(int lev,
                 amrex::TagBoxArray& tags,
                 amrex::Real /*time*/,
                 int /*ngrow*/)
{
    tags.setVal(amrex::TagBox::CLEAR);

    if (lev < 0 ||
        lev >= static_cast<int>(levels_.size()) ||
        levels_[lev] == nullptr ||
        !levels_[lev]->isDefined())
    {
        return;
    }

    const DLMOptions dlm_options = dlmOptionsForCurrentStep();
    if (!dlm_options.enabled || !dlm_options.refine) {
        return;
    }
    if (dlm_options.refinement_max_level >= 0 &&
        lev >= dlm_options.refinement_max_level)
    {
        return;
    }

    DLM::tagRefinement(tags,
                       *levels_[lev],
                       currentDlmSpheres(),
                       dlm_options);
}

void
Solver::ensureLevelCapacity(int lev)
{
    if (lev >= static_cast<int>(levels_.size())) {
        levels_.resize(lev + 1);
        dlm_states_.resize(lev + 1);
        adv_flux_registers_.resize(lev + 1);
        visc_flux_registers_.resize(lev + 1);
        velocity_sync_corrections_.resize(lev + 1);
        sync_projection_residuals_.resize(lev + 1);
    }
}

DLMOptions
Solver::dlmOptionsForCurrentStep() const
{
    DLMOptions options = options_.dlm;
    options.profile_phases = options_.profile_phases;
    return options;
}

std::vector<DLM::RigidSphere>
Solver::currentDlmSpheres() const
{
    if (!standalone_dlm_spheres_.empty()) {
        return standalone_dlm_spheres_;
    }
    return {optionSphere(options_.dlm)};
}

void
Solver::persistStandaloneDlmSpheres(amrex::Real dt)
{
    BL_PROFILE("ELFF::Fluid::Solver::persistStandaloneDlmSpheres");
    if (!options_.dlm.enabled || options_.dlm.scheme != "marker") {
        return;
    }

    const DLMOptions dlm_options = dlmOptionsForCurrentStep();
    const int source_lev = dlm_options.apply_on_finest_level_only ? finestLevel() : 0;
    if (source_lev < 0 ||
        source_lev >= static_cast<int>(dlm_states_.size()) ||
        dlm_states_[source_lev].bodies.empty())
    {
        return;
    }

    standalone_dlm_spheres_.clear();
    standalone_dlm_spheres_.reserve(dlm_states_[source_lev].bodies.size());
    for (DLM::BodyState const& body : dlm_states_[source_lev].bodies) {
        DLM::RigidSphere sphere;
        sphere.id = body.id;
        sphere.geometry_id = body.geometry_id;
        sphere.particle = body.particle;
        sphere.prescribed = body.prescribed;
        sphere.radius = body.radius;
        sphere.density = body.solid_density;
        sphere.mass = body.mass;
        sphere.inertia = body.inertia;
        sphere.center = body.center;
        sphere.velocity = body.velocity;
        sphere.angular_velocity = body.angular_velocity;
        if (!body.prescribed && dlm_options.mode == "free") {
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                sphere.center[dir] += dt * sphere.velocity[dir];
            }
        }
        standalone_dlm_spheres_.push_back(sphere);
    }
}

void
Solver::enforceDlmInteriorVelocity(amrex::Real dt)
{
    BL_PROFILE("ELFF::Fluid::Solver::enforceDlmInteriorVelocity");
    const DLMOptions dlm_options = dlmOptionsForCurrentStep();
    if (!dlm_options.enabled ||
        dlm_options.scheme != "marker" ||
        !dlm_options.enforce_interior_velocity)
    {
        return;
    }

    const int finest = finestLevel();
    for (int lev = 0; lev <= finest; ++lev) {
        if (lev >= static_cast<int>(levels_.size()) ||
            lev >= static_cast<int>(dlm_states_.size()) ||
            levels_[lev] == nullptr)
        {
            continue;
        }
        DLM::enforceInteriorVelocity(dlm_states_[lev], *levels_[lev], dlm_options, dt);
    }
    averageDownHierarchy();
}

} // namespace ELFF::Fluid
