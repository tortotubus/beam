#ifndef ELFF_FLUID_SOLVER_H_
#define ELFF_FLUID_SOLVER_H_

#include <elff/fluid/DLM.H>
#include <elff/fluid/LevelData.H>
#include <elff/fluid/Options.H>

#include <AMReX_AmrCore.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_REAL.H>
#include <AMReX_TagBox.H>

#include <memory>
#include <string>
#include <vector>

namespace ELFF::Fluid {

/** @brief Owns and advances the block-structured incompressible-flow hierarchy.
 *
 * All levels share one timestep. Fluid synchronization uses refluxing,
 * averaging down, and optional composite projections. DLM constraints remain
 * level-local; see @ref architecture for ordering and limitations. */
class Solver final : public amrex::AmrCore
{
public:
    /// Construct the solver and initialize standalone DLM bodies.
    explicit Solver(SolverOptions options);

    /// Allocate the initial hierarchy, initialize fields, and synchronize levels.
    void initData();
    /// Advance until max_step or stop_time, writing configured diagnostics.
    void evolve();

private:
    void MakeNewLevelFromScratch(int lev,
                                 amrex::Real time,
                                 amrex::BoxArray const& ba,
                                 amrex::DistributionMapping const& dm) override;

    void MakeNewLevelFromCoarse(int lev,
                                amrex::Real time,
                                amrex::BoxArray const& ba,
                                amrex::DistributionMapping const& dm) override;

    void RemakeLevel(int lev,
                     amrex::Real time,
                     amrex::BoxArray const& ba,
                     amrex::DistributionMapping const& dm) override;

    void ClearLevel(int lev) override;

    void ErrorEst(int lev,
                  amrex::TagBoxArray& tags,
                  amrex::Real time,
                  int ngrow) override;

    void ensureLevelCapacity(int lev);
    void defineFluxRegisters();
    void resetFluxRegisters();
    void defineVelocitySyncCorrections();
    void resetVelocitySyncCorrections();
    void defineSyncProjectionResiduals();
    void resetSyncProjectionResiduals();
    void snapshotOldStates(amrex::Real time);
    void refluxHierarchy();
    void interpolateVelocitySyncCorrections();
    void fillPatchLevel(int lev, amrex::Real time);
    void fillLevelFromCoarse(int lev);
    void averageDownLevel(int fine_lev);
    void averageDownHierarchy();
    void projectFaceVelocityHierarchy();
    void projectVelocityHierarchy();
    void refreshFaceVelocityHierarchy();

    /// Interpolate reflux corrections, average down, and project cells and faces.
    void synchronizeHierarchy();

    /// Rebuild eligible levels and synchronize at the configured step interval.
    void regridIfNeeded();

    void advanceLevelExplicit(int lev, amrex::Real time, amrex::Real dt);
    void applyViscosityLevel(int lev, amrex::Real dt);
    void applyViscosityHierarchy(amrex::Real dt);

    /// Project one level and apply its DLM solve.
    void projectAndApplyConstraints(int lev, amrex::Real dt);

    void advanceLevel(int lev, amrex::Real time, amrex::Real dt);
    [[nodiscard]] DLMOptions dlmOptionsForCurrentStep() const;
    [[nodiscard]] std::vector<DLM::RigidSphere> currentDlmSpheres() const;
    void persistStandaloneDlmSpheres(amrex::Real dt);
    void enforceDlmInteriorVelocity(amrex::Real dt);
    void writeDlmHistory(amrex::Real output_time) const;
    void writePlotFile() const;
    [[nodiscard]] std::string plotFileName() const;

    SolverOptions options_;
    std::vector<std::unique_ptr<LevelData>> levels_;
    std::vector<DLM::State> dlm_states_;
    std::vector<std::unique_ptr<amrex::FluxRegister>> adv_flux_registers_;
    std::vector<std::unique_ptr<amrex::FluxRegister>> visc_flux_registers_;
    std::vector<std::unique_ptr<amrex::MultiFab>> velocity_sync_corrections_;
    std::vector<std::unique_ptr<amrex::MultiFab>> sync_projection_residuals_;
    std::vector<DLM::RigidSphere> standalone_dlm_spheres_;
    amrex::Real time_ = 0.0;
    int step_ = 0;
};

} // namespace ELFF::Fluid

#endif
