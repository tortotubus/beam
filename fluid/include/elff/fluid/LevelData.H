#ifndef ELFF_FLUID_LEVELDATA_H_
#define ELFF_FLUID_LEVELDATA_H_

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_SPACE.H>

namespace ELFF::Fluid {

/** @brief Owns fluid fields and diagnostics on one AMR level.
 *
 * Velocity, density, and gradp are cell-centered; pressure is a nodal
 * projection potential, not physical pressure. Face velocity is staggered
 * by component. Accessors return borrowed references owned by this object.
 * Plot output converts pressure using density / timestep.
 * @pre Call define() before initializing or accessing field data. */
class LevelData
{
public:
    LevelData() = default;

    /// Allocate fields for the supplied level geometry, grids, and MPI ownership.
    void define(int lev,
                amrex::Geometry geom,
                amrex::BoxArray const& grids,
                amrex::DistributionMapping const& dmap);

    [[nodiscard]] bool isDefined() const noexcept { return defined_; }
    [[nodiscard]] int level() const noexcept { return level_; }
    [[nodiscard]] amrex::Geometry const& geom() const noexcept { return geom_; }
    [[nodiscard]] amrex::Real oldTime() const noexcept { return old_time_; }
    [[nodiscard]] amrex::Real newTime() const noexcept { return new_time_; }

    amrex::MultiFab& velocity() noexcept { return velocity_; }
    amrex::MultiFab& density() noexcept { return density_; }
    amrex::MultiFab& pressure() noexcept { return pressure_; }
    amrex::MultiFab& gradp() noexcept { return gradp_; }
    amrex::MultiFab& oldVelocity() noexcept { return old_velocity_; }
    amrex::MultiFab& oldDensity() noexcept { return old_density_; }
    amrex::MultiFab& oldPressure() noexcept { return old_pressure_; }
    amrex::MultiFab& oldGradp() noexcept { return old_gradp_; }
    amrex::MultiFab& divu() noexcept { return divu_; }
    amrex::MultiFab& rigidMask() noexcept { return rigid_mask_; }
    amrex::MultiFab& rigidVelocity() noexcept { return rigid_velocity_; }
    /// DLM velocity increment divided by dt; plot output multiplies by density.
    amrex::MultiFab& dlmForce() noexcept { return dlm_force_; }
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& faceVelocity() noexcept { return face_velocity_; }
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& faceBeta() noexcept { return face_beta_; }

    amrex::MultiFab const& velocity() const noexcept { return velocity_; }
    amrex::MultiFab const& density() const noexcept { return density_; }
    amrex::MultiFab const& pressure() const noexcept { return pressure_; }
    amrex::MultiFab const& gradp() const noexcept { return gradp_; }
    amrex::MultiFab const& oldVelocity() const noexcept { return old_velocity_; }
    amrex::MultiFab const& oldDensity() const noexcept { return old_density_; }
    amrex::MultiFab const& oldPressure() const noexcept { return old_pressure_; }
    amrex::MultiFab const& oldGradp() const noexcept { return old_gradp_; }
    amrex::MultiFab const& divu() const noexcept { return divu_; }
    amrex::MultiFab const& rigidMask() const noexcept { return rigid_mask_; }
    amrex::MultiFab const& rigidVelocity() const noexcept { return rigid_velocity_; }
    amrex::MultiFab const& dlmForce() const noexcept { return dlm_force_; }
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> const& faceVelocity() const noexcept { return face_velocity_; }
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> const& faceBeta() const noexcept { return face_beta_; }

    /// Initialize fields from prob parameters and set both state timestamps.
    void setInitialState(amrex::Real time);
    /// Copy current cell/nodal state into old fields and set state timestamps.
    void snapshotOldState(amrex::Real time);
    void setNewTime(amrex::Real time) noexcept { new_time_ = time; }
    /// Copy overlapping fields from another layout using ParallelCopy.
    void copyStateFrom(LevelData const& src);

private:
    int level_ = -1;
    bool defined_ = false;
    amrex::Real old_time_ = 0.0;
    amrex::Real new_time_ = 0.0;
    amrex::Geometry geom_;
    amrex::MultiFab velocity_;
    amrex::MultiFab density_;
    amrex::MultiFab pressure_;
    amrex::MultiFab gradp_;
    amrex::MultiFab old_velocity_;
    amrex::MultiFab old_density_;
    amrex::MultiFab old_pressure_;
    amrex::MultiFab old_gradp_;
    amrex::MultiFab divu_;
    amrex::MultiFab rigid_mask_;
    amrex::MultiFab rigid_velocity_;
    amrex::MultiFab dlm_force_;
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> face_velocity_;
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> face_beta_;
};

} // namespace ELFF::Fluid

#endif
