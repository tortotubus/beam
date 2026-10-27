#include "SolverInternals.H"

#include <elff/fluid/Boundary.H>
#include <elff/fluid/Solver.H>

#include <AMReX_FillPatchUtil.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Print.H>

#include <memory>

namespace ELFF::Fluid {

using namespace SolverDetail;

void
Solver::defineFluxRegisters()
{
    const int nlevels = finestLevel() + 1;
    if (nlevels <= 0) {
        adv_flux_registers_.clear();
        visc_flux_registers_.clear();
        return;
    }

    if (static_cast<int>(adv_flux_registers_.size()) < nlevels) {
        adv_flux_registers_.resize(nlevels);
    }
    if (static_cast<int>(visc_flux_registers_.size()) < nlevels) {
        visc_flux_registers_.resize(nlevels);
    }
    adv_flux_registers_[0].reset();
    visc_flux_registers_[0].reset();

    for (int lev = 1; lev < nlevels; ++lev) {
        if (lev >= static_cast<int>(levels_.size()) ||
            levels_[lev] == nullptr ||
            levels_[lev - 1] == nullptr)
        {
            adv_flux_registers_[lev].reset();
            visc_flux_registers_[lev].reset();
            continue;
        }

        adv_flux_registers_[lev] =
            std::make_unique<amrex::FluxRegister>(levels_[lev]->velocity().boxArray(),
                                                  levels_[lev]->velocity().DistributionMap(),
                                                  refRatio(lev - 1),
                                                  lev,
                                                  AMREX_SPACEDIM);
        visc_flux_registers_[lev] =
            std::make_unique<amrex::FluxRegister>(levels_[lev]->velocity().boxArray(),
                                                  levels_[lev]->velocity().DistributionMap(),
                                                  refRatio(lev - 1),
                                                  lev,
                                                  AMREX_SPACEDIM);
    }

    for (int lev = nlevels; lev < static_cast<int>(adv_flux_registers_.size()); ++lev) {
        adv_flux_registers_[lev].reset();
    }
    for (int lev = nlevels; lev < static_cast<int>(visc_flux_registers_.size()); ++lev) {
        visc_flux_registers_[lev].reset();
    }
}

void
Solver::resetFluxRegisters()
{
    for (std::unique_ptr<amrex::FluxRegister>& flux_register : adv_flux_registers_) {
        if (flux_register != nullptr) {
            flux_register->setVal(0.0);
        }
    }
    for (std::unique_ptr<amrex::FluxRegister>& flux_register : visc_flux_registers_) {
        if (flux_register != nullptr) {
            flux_register->setVal(0.0);
        }
    }
}

void
Solver::defineVelocitySyncCorrections()
{
    const int nlevels = finestLevel() + 1;
    if (nlevels <= 0) {
        velocity_sync_corrections_.clear();
        return;
    }

    if (static_cast<int>(velocity_sync_corrections_.size()) < nlevels) {
        velocity_sync_corrections_.resize(nlevels);
    }

    for (int lev = 0; lev < nlevels; ++lev) {
        if (lev >= static_cast<int>(levels_.size()) || levels_[lev] == nullptr) {
            velocity_sync_corrections_[lev].reset();
            continue;
        }

        velocity_sync_corrections_[lev] =
            std::make_unique<amrex::MultiFab>(levels_[lev]->velocity().boxArray(),
                                              levels_[lev]->velocity().DistributionMap(),
                                              AMREX_SPACEDIM,
                                              1,
                                              amrex::MFInfo(),
                                              levels_[lev]->velocity().Factory());
        velocity_sync_corrections_[lev]->setVal(0.0,
                                                0,
                                                AMREX_SPACEDIM,
                                                velocity_sync_corrections_[lev]->nGrowVect());
    }

    for (int lev = nlevels; lev < static_cast<int>(velocity_sync_corrections_.size()); ++lev) {
        velocity_sync_corrections_[lev].reset();
    }
}

void
Solver::resetVelocitySyncCorrections()
{
    for (std::unique_ptr<amrex::MultiFab>& correction : velocity_sync_corrections_) {
        if (correction != nullptr) {
            correction->setVal(0.0, 0, AMREX_SPACEDIM, correction->nGrowVect());
        }
    }
}

void
Solver::defineSyncProjectionResiduals()
{
    const int nlevels = finestLevel() + 1;
    if (nlevels <= 0) {
        sync_projection_residuals_.clear();
        return;
    }

    if (static_cast<int>(sync_projection_residuals_.size()) < nlevels) {
        sync_projection_residuals_.resize(nlevels);
    }

    for (int lev = 0; lev < nlevels; ++lev) {
        if (lev >= static_cast<int>(levels_.size()) || levels_[lev] == nullptr) {
            sync_projection_residuals_[lev].reset();
            continue;
        }

        sync_projection_residuals_[lev] =
            std::make_unique<amrex::MultiFab>(levels_[lev]->pressure().boxArray(),
                                              levels_[lev]->pressure().DistributionMap(),
                                              1,
                                              levels_[lev]->pressure().nGrowVect(),
                                              amrex::MFInfo(),
                                              levels_[lev]->pressure().Factory());
        sync_projection_residuals_[lev]->setVal(0.0,
                                                0,
                                                1,
                                                sync_projection_residuals_[lev]->nGrowVect());
    }

    for (int lev = nlevels; lev < static_cast<int>(sync_projection_residuals_.size()); ++lev) {
        sync_projection_residuals_[lev].reset();
    }
}

void
Solver::resetSyncProjectionResiduals()
{
    for (std::unique_ptr<amrex::MultiFab>& residual : sync_projection_residuals_) {
        if (residual != nullptr) {
            residual->setVal(0.0, 0, 1, residual->nGrowVect());
        }
    }
}

void
Solver::snapshotOldStates(amrex::Real time)
{
    for (int lev = 0; lev <= finestLevel(); ++lev) {
        if (lev < static_cast<int>(levels_.size()) && levels_[lev] != nullptr) {
            levels_[lev]->snapshotOldState(time);
        }
    }
}

void
Solver::refluxHierarchy()
{
    if ((!options_.advective_reflux && !options_.viscous_reflux) ||
        finestLevel() <= 0)
    {
        return;
    }

    for (int fine_lev = finestLevel(); fine_lev > 0; --fine_lev) {
        if (levels_[fine_lev - 1] == nullptr)
        {
            continue;
        }

        const bool have_advective_register =
            options_.advective_reflux &&
            fine_lev < static_cast<int>(adv_flux_registers_.size()) &&
            adv_flux_registers_[fine_lev] != nullptr;
        const bool have_viscous_register =
            options_.viscous_reflux &&
            fine_lev < static_cast<int>(visc_flux_registers_.size()) &&
            visc_flux_registers_[fine_lev] != nullptr;
        if (!have_advective_register && !have_viscous_register) {
            continue;
        }

        LevelData& crse = *levels_[fine_lev - 1];
        amrex::MultiFab before(crse.velocity().boxArray(),
                               crse.velocity().DistributionMap(),
                               AMREX_SPACEDIM,
                               0,
                               amrex::MFInfo(),
                               crse.velocity().Factory());
        amrex::MultiFab::Copy(before, crse.velocity(), 0, 0, AMREX_SPACEDIM, 0);

        if (have_advective_register) {
            adv_flux_registers_[fine_lev]->Reflux(crse.velocity(),
                                                  amrex::Real(1.0),
                                                  0,
                                                  0,
                                                  AMREX_SPACEDIM,
                                                  crse.geom());
        }
        if (have_viscous_register) {
            visc_flux_registers_[fine_lev]->Reflux(crse.velocity(),
                                                   amrex::Real(1.0),
                                                   0,
                                                   0,
                                                   AMREX_SPACEDIM,
                                                   crse.geom());
        }
        Boundary::fillVelocityPhysicalBoundary(crse, options_);

        if (fine_lev - 1 < static_cast<int>(velocity_sync_corrections_.size()) &&
            velocity_sync_corrections_[fine_lev - 1] != nullptr)
        {
            amrex::MultiFab delta(crse.velocity().boxArray(),
                                  crse.velocity().DistributionMap(),
                                  AMREX_SPACEDIM,
                                  0,
                                  amrex::MFInfo(),
                                  crse.velocity().Factory());
            amrex::MultiFab::Copy(delta, crse.velocity(), 0, 0, AMREX_SPACEDIM, 0);
            amrex::MultiFab::Subtract(delta, before, 0, 0, AMREX_SPACEDIM, 0);
            amrex::MultiFab::Add(*velocity_sync_corrections_[fine_lev - 1],
                                 delta,
                                 0,
                                 0,
                                 AMREX_SPACEDIM,
                                 0);
        }
    }
}

void
Solver::interpolateVelocitySyncCorrections()
{
    if (!options_.sync_interpolation || finestLevel() <= 0) {
        return;
    }

    const amrex::IntVect no_grow(0);
    const amrex::IntVect outside_domain_ghosts(0);
    for (int crse_lev = 0; crse_lev < finestLevel(); ++crse_lev) {
        const int fine_lev = crse_lev + 1;
        if (crse_lev >= static_cast<int>(velocity_sync_corrections_.size()) ||
            fine_lev >= static_cast<int>(velocity_sync_corrections_.size()) ||
            velocity_sync_corrections_[crse_lev] == nullptr ||
            velocity_sync_corrections_[fine_lev] == nullptr ||
            levels_[fine_lev] == nullptr)
        {
            continue;
        }

        LevelData& fine = *levels_[fine_lev];
        amrex::MultiFab fine_increment(fine.velocity().boxArray(),
                                       fine.velocity().DistributionMap(),
                                       AMREX_SPACEDIM,
                                       0,
                                       amrex::MFInfo(),
                                       fine.velocity().Factory());
        fine_increment.setVal(0.0);
        velocity_sync_corrections_[crse_lev]->FillBoundary(Geom(crse_lev).periodicity());

        auto const velocity_bcs = Boundary::makeVelocityBCs(Geom(fine_lev), options_);
        amrex::InterpFromCoarseLevel(fine_increment,
                                     no_grow,
                                     outside_domain_ghosts,
                                     *velocity_sync_corrections_[crse_lev],
                                     0,
                                     0,
                                     AMREX_SPACEDIM,
                                     Geom(crse_lev),
                                     Geom(fine_lev),
                                     refRatio(crse_lev),
                                     &amrex::cell_cons_interp,
                                     velocity_bcs,
                                     0);

        amrex::MultiFab::Add(fine.velocity(), fine_increment, 0, 0, AMREX_SPACEDIM, 0);
        amrex::MultiFab::Add(*velocity_sync_corrections_[fine_lev],
                             fine_increment,
                             0,
                             0,
                             AMREX_SPACEDIM,
                             0);
        Boundary::fillVelocityPhysicalBoundary(fine, options_);
    }
}

void
Solver::fillPatchLevel(int lev, amrex::Real time)
{
    AMREX_ALWAYS_ASSERT(lev >= 0);
    AMREX_ALWAYS_ASSERT(lev < static_cast<int>(levels_.size()));
    AMREX_ALWAYS_ASSERT(levels_[lev] != nullptr);

    LevelData& fine = *levels_[lev];
    amrex::Vector<amrex::MultiFab*> fvel{&fine.oldVelocity(), &fine.velocity()};
    amrex::Vector<amrex::MultiFab*> frho{&fine.oldDensity(), &fine.density()};
    amrex::Vector<amrex::MultiFab*> fp{&fine.oldPressure(), &fine.pressure()};
    amrex::Vector<amrex::MultiFab*> fgp{&fine.oldGradp(), &fine.gradp()};
    amrex::Vector<amrex::Real> ftime{fine.oldTime(), fine.newTime()};

    auto const velocity_bcs = Boundary::makeVelocityBCs(Geom(lev), options_);
    auto const scalar_bcs = makeExtrapolationBCs(Geom(lev), 1);
    auto const vector_bcs = makeExtrapolationBCs(Geom(lev), AMREX_SPACEDIM);

    if (lev == 0) {
        amrex::FillPatchSingleLevel(fine.velocity(),
                                    fine.velocity().nGrowVect(),
                                    time,
                                    fvel,
                                    fine.velocity().nGrowVect(),
                                    ftime,
                                    0,
                                    0,
                                    AMREX_SPACEDIM,
                                    Geom(lev));
        amrex::FillPatchSingleLevel(fine.density(),
                                    fine.density().nGrowVect(),
                                    time,
                                    frho,
                                    fine.density().nGrowVect(),
                                    ftime,
                                    0,
                                    0,
                                    1,
                                    Geom(lev));
        amrex::FillPatchSingleLevel(fine.gradp(),
                                    fine.gradp().nGrowVect(),
                                    time,
                                    fgp,
                                    fine.gradp().nGrowVect(),
                                    ftime,
                                    0,
                                    0,
                                    AMREX_SPACEDIM,
                                    Geom(lev));
        amrex::FillPatchSingleLevel(fine.pressure(),
                                    fine.pressure().nGrowVect(),
                                    time,
                                    fp,
                                    fine.pressure().nGrowVect(),
                                    ftime,
                                    0,
                                    0,
                                    1,
                                    Geom(lev));
    } else {
        AMREX_ALWAYS_ASSERT(levels_[lev - 1] != nullptr);
        LevelData& crse = *levels_[lev - 1];
        amrex::Vector<amrex::MultiFab*> cvel{&crse.oldVelocity(), &crse.velocity()};
        amrex::Vector<amrex::MultiFab*> crho{&crse.oldDensity(), &crse.density()};
        amrex::Vector<amrex::MultiFab*> cp{&crse.oldPressure(), &crse.pressure()};
        amrex::Vector<amrex::MultiFab*> cgp{&crse.oldGradp(), &crse.gradp()};
        amrex::Vector<amrex::Real> ctime{crse.oldTime(), crse.newTime()};

        const amrex::IntVect outside_domain_ghosts(0);
        amrex::FillPatchTwoLevels(fine.velocity(),
                                  fine.velocity().nGrowVect(),
                                  outside_domain_ghosts,
                                  time,
                                  cvel,
                                  ctime,
                                  fvel,
                                  ftime,
                                  0,
                                  0,
                                  AMREX_SPACEDIM,
                                  Geom(lev - 1),
                                  Geom(lev),
                                  refRatio(lev - 1),
                                  &amrex::cell_cons_interp,
                                  velocity_bcs,
                                  0);
        amrex::FillPatchTwoLevels(fine.density(),
                                  fine.density().nGrowVect(),
                                  outside_domain_ghosts,
                                  time,
                                  crho,
                                  ctime,
                                  frho,
                                  ftime,
                                  0,
                                  0,
                                  1,
                                  Geom(lev - 1),
                                  Geom(lev),
                                  refRatio(lev - 1),
                                  &amrex::cell_cons_interp,
                                  scalar_bcs,
                                  0);
        amrex::FillPatchTwoLevels(fine.gradp(),
                                  fine.gradp().nGrowVect(),
                                  outside_domain_ghosts,
                                  time,
                                  cgp,
                                  ctime,
                                  fgp,
                                  ftime,
                                  0,
                                  0,
                                  AMREX_SPACEDIM,
                                  Geom(lev - 1),
                                  Geom(lev),
                                  refRatio(lev - 1),
                                  &amrex::cell_cons_interp,
                                  vector_bcs,
                                  0);
        amrex::FillPatchTwoLevels(fine.pressure(),
                                  fine.pressure().nGrowVect(),
                                  outside_domain_ghosts,
                                  time,
                                  cp,
                                  ctime,
                                  fp,
                                  ftime,
                                  0,
                                  0,
                                  1,
                                  Geom(lev - 1),
                                  Geom(lev),
                                  refRatio(lev - 1),
                                  &amrex::node_bilinear_interp,
                                  scalar_bcs,
                                  0);
    }

    Boundary::fillVelocityPhysicalBoundary(fine, options_);
    fine.density().FillBoundary(fine.geom().periodicity());
    fine.pressure().FillBoundary(fine.geom().periodicity());
    fine.gradp().FillBoundary(fine.geom().periodicity());
}

void
Solver::fillLevelFromCoarse(int lev)
{
    AMREX_ALWAYS_ASSERT(lev > 0);
    AMREX_ALWAYS_ASSERT(lev < static_cast<int>(levels_.size()));
    AMREX_ALWAYS_ASSERT(levels_[lev] != nullptr);
    AMREX_ALWAYS_ASSERT(levels_[lev - 1] != nullptr);

    LevelData& fine = *levels_[lev];
    LevelData const& crse = *levels_[lev - 1];
    amrex::Geometry const& cgeom = Geom(lev - 1);
    amrex::Geometry const& fgeom = Geom(lev);
    amrex::IntVect const& ratio = refRatio(lev - 1);
    amrex::IntVect const outside_domain_ghosts(0);

    auto const velocity_bcs = Boundary::makeVelocityBCs(fgeom, options_);
    amrex::InterpFromCoarseLevel(fine.velocity(),
                                 fine.velocity().nGrowVect(),
                                 outside_domain_ghosts,
                                 crse.velocity(),
                                 0,
                                 0,
                                 AMREX_SPACEDIM,
                                 cgeom,
                                 fgeom,
                                 ratio,
                                 &amrex::cell_cons_interp,
                                 velocity_bcs,
                                 0);

    auto const scalar_bcs = makeExtrapolationBCs(fgeom, 1);
    amrex::InterpFromCoarseLevel(fine.density(),
                                 fine.density().nGrowVect(),
                                 outside_domain_ghosts,
                                 crse.density(),
                                 0,
                                 0,
                                 1,
                                 cgeom,
                                 fgeom,
                                 ratio,
                                 &amrex::cell_cons_interp,
                                 scalar_bcs,
                                 0);

    auto const vector_bcs = makeExtrapolationBCs(fgeom, AMREX_SPACEDIM);
    amrex::InterpFromCoarseLevel(fine.gradp(),
                                 fine.gradp().nGrowVect(),
                                 outside_domain_ghosts,
                                 crse.gradp(),
                                 0,
                                 0,
                                 AMREX_SPACEDIM,
                                 cgeom,
                                 fgeom,
                                 ratio,
                                 &amrex::cell_cons_interp,
                                 vector_bcs,
                                 0);

    amrex::InterpFromCoarseLevel(fine.pressure(),
                                 fine.pressure().nGrowVect(),
                                 outside_domain_ghosts,
                                 crse.pressure(),
                                 0,
                                 0,
                                 1,
                                 cgeom,
                                 fgeom,
                                 ratio,
                                 &amrex::node_bilinear_interp,
                                 scalar_bcs,
                                 0);

    fine.snapshotOldState(fine.newTime());
    fine.divu().setVal(0.0);
    fine.rigidMask().setVal(0.0);
    fine.rigidVelocity().setVal(0.0);
    fine.dlmForce().setVal(0.0);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        fine.faceVelocity()[dir].setVal(0.0);
        fine.faceBeta()[dir].setVal(1.0);
    }
}

void
Solver::averageDownLevel(int fine_lev)
{
    if (fine_lev <= 0 || fine_lev > finestLevel()) {
        return;
    }
    if (fine_lev >= static_cast<int>(levels_.size()) ||
        levels_[fine_lev] == nullptr ||
        levels_[fine_lev - 1] == nullptr)
    {
        return;
    }

    LevelData const& fine = *levels_[fine_lev];
    LevelData& crse = *levels_[fine_lev - 1];
    amrex::IntVect const& ratio = refRatio(fine_lev - 1);

    amrex::average_down(fine.velocity(),
                        crse.velocity(),
                        fine.geom(),
                        crse.geom(),
                        0,
                        AMREX_SPACEDIM,
                        ratio);
    amrex::average_down(fine.density(),
                        crse.density(),
                        fine.geom(),
                        crse.geom(),
                        0,
                        1,
                        ratio);
    amrex::average_down(fine.gradp(),
                        crse.gradp(),
                        fine.geom(),
                        crse.geom(),
                        0,
                        AMREX_SPACEDIM,
                        ratio);
    amrex::average_down(fine.divu(),
                        crse.divu(),
                        fine.geom(),
                        crse.geom(),
                        0,
                        1,
                        ratio);
    amrex::average_down(fine.rigidMask(),
                        crse.rigidMask(),
                        fine.geom(),
                        crse.geom(),
                        0,
                        1,
                        ratio);
    amrex::average_down(fine.rigidVelocity(),
                        crse.rigidVelocity(),
                        fine.geom(),
                        crse.geom(),
                        0,
                        AMREX_SPACEDIM,
                        ratio);
    amrex::average_down(fine.dlmForce(),
                        crse.dlmForce(),
                        fine.geom(),
                        crse.geom(),
                        0,
                        AMREX_SPACEDIM,
                        ratio);
    amrex::average_down(fine.pressure(),
                        crse.pressure(),
                        0,
                        1,
                        ratio);
    amrex::average_down_faces(amrex::GetArrOfConstPtrs(fine.faceVelocity()),
                              amrex::GetArrOfPtrs(crse.faceVelocity()),
                              ratio,
                              crse.geom());
    amrex::average_down_faces(amrex::GetArrOfConstPtrs(fine.faceBeta()),
                              amrex::GetArrOfPtrs(crse.faceBeta()),
                              ratio,
                              crse.geom());

    Boundary::fillVelocityPhysicalBoundary(crse, options_);
    crse.density().FillBoundary(crse.geom().periodicity());
    crse.pressure().FillBoundary(crse.geom().periodicity());
    crse.gradp().FillBoundary(crse.geom().periodicity());
}

void
Solver::averageDownHierarchy()
{
    for (int fine_lev = finestLevel(); fine_lev > 0; --fine_lev) {
        averageDownLevel(fine_lev);
    }
}

void
Solver::regridIfNeeded()
{
    if (maxLevel() <= 0 || options_.regrid_int <= 0 ||
        step_ % options_.regrid_int != 0)
    {
        return;
    }

    const int old_finest = finestLevel();
    regrid(0, time_);
    defineFluxRegisters();
    defineVelocitySyncCorrections();
    defineSyncProjectionResiduals();
    synchronizeHierarchy();

    if (options_.verbose > 0 && finestLevel() != old_finest) {
        amrex::Print() << "Regridded hierarchy from finest level "
                       << old_finest << " to " << finestLevel() << "\n";
    }
}

} // namespace ELFF::Fluid
