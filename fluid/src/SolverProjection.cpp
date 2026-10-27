#include "SolverInternals.H"

#include <elff/fluid/Boundary.H>
#include <elff/fluid/Solver.H>

#include <AMReX_MultiFabUtil.H>

#include <hydro_MacProjector.H>
#include <hydro_NodalProjector.H>

#include <algorithm>

namespace ELFF::Fluid {

using namespace SolverDetail;

void
Solver::projectVelocityHierarchy()
{
    if (!options_.composite_projection || finestLevel() < 0) {
        return;
    }

    const int nlevels = finestLevel() + 1;
    amrex::Vector<amrex::MultiFab*> velocities;
    amrex::Vector<amrex::MultiFab*> cell_sources;
    amrex::Vector<amrex::MultiFab const*> nodal_sources;
    amrex::Vector<amrex::Geometry> geoms;
    velocities.reserve(nlevels);
    cell_sources.reserve(nlevels);
    nodal_sources.reserve(nlevels);
    geoms.reserve(nlevels);
    const DLMOptions dlm_options = dlmOptionsForCurrentStep();
    const bool use_sync_residuals =
        options_.sync_projection_residuals &&
        !dlm_options.enabled &&
        static_cast<int>(sync_projection_residuals_.size()) >= nlevels;

    for (int lev = 0; lev < nlevels; ++lev) {
        AMREX_ALWAYS_ASSERT(lev < static_cast<int>(levels_.size()));
        AMREX_ALWAYS_ASSERT(levels_[lev] != nullptr);

        Boundary::fillVelocityPhysicalBoundary(*levels_[lev], options_);
        velocities.push_back(&levels_[lev]->velocity());
        if (use_sync_residuals) {
            nodal_sources.push_back(sync_projection_residuals_[lev].get());
        }
        geoms.push_back(Geom(lev));
    }

    amrex::LPInfo lp_info;
    Hydro::NodalProjector nodal_projector(velocities,
                                          amrex::Real(1.0),
                                          geoms,
                                          lp_info,
                                          cell_sources,
                                          nodal_sources);
    auto [lo, hi] = Boundary::projectionBCs(Geom(0), options_);
    nodal_projector.setDomainBC(lo, hi);
    nodal_projector.setVerbose(std::max(0, options_.verbose - 1));
    nodal_projector.project(options_.projection_rel_tol,
                            options_.projection_abs_tol);

    auto phi = nodal_projector.getPhi();
    auto grad_phi = nodal_projector.getGradPhi();
    for (int lev = 0; lev < nlevels; ++lev) {
        if (lev < static_cast<int>(phi.size()) && phi[lev] != nullptr) {
            levels_[lev]->pressure().ParallelCopy(*phi[lev],
                                                  0,
                                                  0,
                                                  1,
                                                  0,
                                                  0,
                                                  Geom(lev).periodicity());
        }
        if (lev < static_cast<int>(grad_phi.size()) && grad_phi[lev] != nullptr) {
            amrex::MultiFab::Copy(levels_[lev]->gradp(),
                                  *grad_phi[lev],
                                  0,
                                  0,
                                  AMREX_SPACEDIM,
                                  0);
        }
        Boundary::fillVelocityPhysicalBoundary(*levels_[lev], options_);
        levels_[lev]->pressure().FillBoundary(Geom(lev).periodicity());
        levels_[lev]->gradp().FillBoundary(Geom(lev).periodicity());
    }
}

void
Solver::projectFaceVelocityHierarchy()
{
    if (!options_.mac_sync || finestLevel() <= 0) {
        for (int lev = 0; lev <= finestLevel(); ++lev) {
            if (lev < static_cast<int>(levels_.size()) && levels_[lev] != nullptr) {
                projectFaceVelocity(*levels_[lev], options_);
            }
        }
        return;
    }

    const int nlevels = finestLevel() + 1;
    amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>> face_velocities;
    amrex::Vector<amrex::Array<amrex::MultiFab const*, AMREX_SPACEDIM>> face_betas;
    amrex::Vector<amrex::Geometry> geoms;
    face_velocities.reserve(nlevels);
    face_betas.reserve(nlevels);
    geoms.reserve(nlevels);

    for (int lev = 0; lev < nlevels; ++lev) {
        AMREX_ALWAYS_ASSERT(lev < static_cast<int>(levels_.size()));
        AMREX_ALWAYS_ASSERT(levels_[lev] != nullptr);

        face_velocities.push_back(amrex::GetArrOfPtrs(levels_[lev]->faceVelocity()));
        face_betas.push_back(amrex::GetArrOfConstPtrs(levels_[lev]->faceBeta()));
        geoms.push_back(Geom(lev));
    }

    amrex::LPInfo lp_info;
    Hydro::MacProjector mac_projector(face_velocities,
                                      face_betas,
                                      geoms,
                                      lp_info);
    auto [lo, hi] = Boundary::projectionBCs(Geom(0), options_);
    mac_projector.setDomainBC(lo, hi);
    mac_projector.setVerbose(std::max(0, options_.verbose - 1));
    mac_projector.project(options_.projection_rel_tol,
                          options_.projection_abs_tol);

    for (int lev = 0; lev < nlevels; ++lev) {
        Boundary::fillFaceVelocityPhysicalBoundary(levels_[lev]->faceVelocity(),
                                                   Geom(lev),
                                                   options_);
        amrex::computeDivergence(levels_[lev]->divu(),
                                 amrex::GetArrOfConstPtrs(levels_[lev]->faceVelocity()),
                                 Geom(lev));
    }
}

void
Solver::refreshFaceVelocityHierarchy()
{
    for (int lev = 0; lev <= finestLevel(); ++lev) {
        if (lev < static_cast<int>(levels_.size()) && levels_[lev] != nullptr) {
            extrapolateFaceVelocityFromCell(*levels_[lev], options_);
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                levels_[lev]->faceBeta()[dir].setVal(1.0);
            }
            Boundary::fillFaceVelocityPhysicalBoundary(levels_[lev]->faceVelocity(),
                                                       Geom(lev),
                                                       options_);
        }
    }
    projectFaceVelocityHierarchy();
}

void
Solver::synchronizeHierarchy()
{
    if (finestLevel() < 0) {
        return;
    }

    interpolateVelocitySyncCorrections();
    averageDownHierarchy();
    projectVelocityHierarchy();
    refreshFaceVelocityHierarchy();
    averageDownHierarchy();
}

} // namespace ELFF::Fluid
