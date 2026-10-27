#include "SolverInternals.H"

#include <elff/fluid/Boundary.H>
#include <elff/fluid/DLM.H>
#include <elff/fluid/Solver.H>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_MFIter.H>
#include <AMReX_MLLinOp.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLTensorOp.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

#include <hydro_mol.H>
#include <hydro_utils.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <vector>

namespace ELFF::Fluid {

using namespace SolverDetail;

namespace {

void
computeMOLAdvectiveTerm(LevelData& state,
                        SolverOptions const& options,
                        amrex::MultiFab& advective_term,
                        amrex::FluxRegister* crse_flux_register,
                        amrex::FluxRegister* fine_flux_register,
                        amrex::Real dt)
{
    BL_PROFILE("ELFF::Fluid::computeMOLAdvectiveTerm");
    auto velocity_bcs = Boundary::makeVelocityBCs(state.geom(), options);
    auto device_velocity_bcs = Boundary::makeDeviceBCs(velocity_bcs);

    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> edge_state;
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> fluxes;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        edge_state[dir].define(state.faceVelocity()[dir].boxArray(),
                               state.faceVelocity()[dir].DistributionMap(),
                               AMREX_SPACEDIM,
                               0);
        fluxes[dir].define(state.faceVelocity()[dir].boxArray(),
                           state.faceVelocity()[dir].DistributionMap(),
                           AMREX_SPACEDIM,
                           0);
    }

    amrex::Gpu::DeviceVector<int> iconserv(AMREX_SPACEDIM, 1);
    int const* iconserv_ptr = iconserv.dataPtr();
    amrex::BCRec const* bcrec_ptr = device_velocity_bcs.dataPtr();

    advective_term.setVal(0.0);

#ifdef AMREX_USE_GPU
    amrex::Gpu::streamSynchronize();
#endif

    for (amrex::MFIter mfi(state.velocity(), amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();

        MOL::ComputeEdgeState(bx,
                              AMREX_D_DECL(edge_state[0].array(mfi),
                                           edge_state[1].array(mfi),
                                           edge_state[2].array(mfi)),
                              state.velocity().const_array(mfi),
                              AMREX_SPACEDIM,
                              AMREX_D_DECL(state.faceVelocity()[0].const_array(mfi),
                                           state.faceVelocity()[1].const_array(mfi),
                                           state.faceVelocity()[2].const_array(mfi)),
                              state.geom().Domain(),
                              velocity_bcs,
                              bcrec_ptr,
                              true);

        HydroUtils::ComputeFluxes(bx,
                                  AMREX_D_DECL(fluxes[0].array(mfi),
                                               fluxes[1].array(mfi),
                                               fluxes[2].array(mfi)),
                                  AMREX_D_DECL(state.faceVelocity()[0].const_array(mfi),
                                               state.faceVelocity()[1].const_array(mfi),
                                               state.faceVelocity()[2].const_array(mfi)),
                                  AMREX_D_DECL(edge_state[0].const_array(mfi),
                                               edge_state[1].const_array(mfi),
                                               edge_state[2].const_array(mfi)),
                                  state.geom(),
                                  AMREX_SPACEDIM,
                                  true,
                                  iconserv_ptr);

        HydroUtils::ComputeDivergence(bx,
                                      advective_term.array(mfi),
                                      AMREX_D_DECL(fluxes[0].const_array(mfi),
                                                   fluxes[1].const_array(mfi),
                                                   fluxes[2].const_array(mfi)),
                                      AMREX_SPACEDIM,
                                      state.geom(),
                                      amrex::Real(1.0),
                                      true);
    }

    if (crse_flux_register != nullptr) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            crse_flux_register->CrseInit(fluxes[dir],
                                         dir,
                                         0,
                                         0,
                                         AMREX_SPACEDIM,
                                         -dt);
        }
    }
    if (fine_flux_register != nullptr) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            fine_flux_register->FineAdd(fluxes[dir],
                                        dir,
                                        0,
                                        0,
                                        AMREX_SPACEDIM,
                                        dt);
        }
    }
}

void
applyExplicitVelocityUpdate(LevelData& state,
                            amrex::MultiFab const& advective_term,
                            amrex::Real dt)
{
    for (amrex::MFIter mfi(state.velocity(), amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto const vel = state.velocity().array(mfi);
        auto const adv = advective_term.const_array(mfi);

        amrex::ParallelFor(bx, AMREX_SPACEDIM,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
            vel(i,j,k,n) -= dt * adv(i,j,k,n);
        });
    }
}

bool
hasNonzeroDriving(std::array<amrex::Real, 3> const& values)
{
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (std::abs(values[dir]) > 0.0) {
            return true;
        }
    }
    return false;
}

void
applyDrivingForces(LevelData& state, SolverOptions const& options, amrex::Real dt)
{
    if (!hasNonzeroDriving(options.body_force) &&
        !hasNonzeroDriving(options.pressure_gradient))
    {
        return;
    }

    const amrex::Real body_force0 = options.body_force[0];
    const amrex::Real body_force1 = options.body_force[1];
    const amrex::Real body_force2 = options.body_force[2];
    const amrex::Real pressure_gradient0 = options.pressure_gradient[0];
    const amrex::Real pressure_gradient1 = options.pressure_gradient[1];
    const amrex::Real pressure_gradient2 = options.pressure_gradient[2];

    for (amrex::MFIter mfi(state.velocity(), amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto const vel = state.velocity().array(mfi);
        auto const rho = state.density().const_array(mfi);

        amrex::ParallelFor(bx, AMREX_SPACEDIM,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
            const amrex::Real raw_density = rho(i,j,k);
            const amrex::Real density =
                raw_density > amrex::Real(1.0e-30) ? raw_density : amrex::Real(1.0e-30);
            const amrex::Real force =
                (n == 0) ? body_force0 :
                (n == 1) ? body_force1 : body_force2;
            const amrex::Real gradp =
                (n == 0) ? pressure_gradient0 :
                (n == 1) ? pressure_gradient1 : pressure_gradient2;
            vel(i,j,k,n) += dt * (force - gradp / density);
        });
    }
}

void
addCellVelocityIncrementToFaces(LevelData& state)
{
    const amrex::Box& domain = state.geom().Domain();

    for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
        for (amrex::MFIter mfi(state.faceVelocity()[comp], amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            auto const face = state.faceVelocity()[comp].array(mfi);
            auto const velocity = state.velocity().const_array(mfi);
            auto const old_velocity = state.oldVelocity().const_array(mfi);
            const int lo = domain.smallEnd(comp);
            const int hi = domain.bigEnd(comp);
            const bool periodic = state.geom().isPeriodic(comp);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                int li = i;
                int lj = j;
                int lk = k;
                int ri = i;
                int rj = j;
                int rk = k;

                if (comp == 0) {
                    --li;
                    if (!periodic) {
                        li = amrex::max(li, lo);
                        ri = amrex::min(ri, hi);
                    }
                } else if (comp == 1) {
                    --lj;
                    if (!periodic) {
                        lj = amrex::max(lj, lo);
                        rj = amrex::min(rj, hi);
                    }
#if (AMREX_SPACEDIM == 3)
                } else {
                    --lk;
                    if (!periodic) {
                        lk = amrex::max(lk, lo);
                        rk = amrex::min(rk, hi);
                    }
#endif
                }

                const amrex::Real left_increment =
                    velocity(li,lj,lk,comp) - old_velocity(li,lj,lk,comp);
                const amrex::Real right_increment =
                    velocity(ri,rj,rk,comp) - old_velocity(ri,rj,rk,comp);
                face(i,j,k) += amrex::Real(0.5) *
                               (left_increment + right_increment);
            });
        }
    }
}

void
defineVectorFaceFluxesLike(LevelData const& state,
                           DLM::FaceFields& fluxes)
{
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        fluxes[dir].define(state.faceVelocity()[dir].boxArray(),
                           state.faceVelocity()[dir].DistributionMap(),
                           AMREX_SPACEDIM,
                           0,
                           amrex::MFInfo(),
                           state.faceVelocity()[dir].Factory());
        fluxes[dir].setVal(0.0);
    }
}

void
scaleFluxesByFaceArea(amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& fluxes,
                      amrex::Geometry const& geom,
                      amrex::Real scale)
{
    const auto dx = geom.CellSizeArray();
    amrex::GpuArray<amrex::Real, 3> area{1.0, 1.0, 1.0};
#if (AMREX_SPACEDIM == 2)
    area[0] = dx[1];
    area[1] = dx[0];
#else
    area[0] = dx[1] * dx[2];
    area[1] = dx[0] * dx[2];
    area[2] = dx[0] * dx[1];
#endif

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        for (amrex::MFIter mfi(fluxes[dir], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.tilebox();
            auto const flux = fluxes[dir].array(mfi);
            const amrex::Real factor = scale * area[dir];

            amrex::ParallelFor(bx, AMREX_SPACEDIM,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                flux(i,j,k,n) *= factor;
            });
        }
    }
}

void
computeExtensiveViscousFluxes(amrex::MLMG& mlmg,
                              LevelData& state,
                              DLM::FaceFields& fluxes,
                              amrex::Real scale)
{
    defineVectorFaceFluxesLike(state, fluxes);
    mlmg.getFluxes({amrex::GetArrOfPtrs(fluxes)},
                   {&state.velocity()},
                   amrex::MLLinOp::Location::FaceCenter);
    scaleFluxesByFaceArea(fluxes, state.geom(), scale);
}

void
registerFluxes(amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> const& fluxes,
               amrex::FluxRegister* crse_flux_register,
               amrex::FluxRegister* fine_flux_register,
               amrex::Real dt)
{
    if (crse_flux_register != nullptr) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            crse_flux_register->CrseInit(fluxes[dir],
                                         dir,
                                         0,
                                         0,
                                         AMREX_SPACEDIM,
                                         -dt);
        }
    }
    if (fine_flux_register != nullptr) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            fine_flux_register->FineAdd(fluxes[dir],
                                        dir,
                                        0,
                                        0,
                                        AMREX_SPACEDIM,
                                        dt);
        }
    }
}

void
applyImplicitViscosity(LevelData& state,
                       SolverOptions const& options,
                       amrex::FluxRegister* crse_flux_register,
                       amrex::FluxRegister* fine_flux_register,
                       amrex::Real dt)
{
    BL_PROFILE("ELFF::Fluid::applyImplicitViscosity");
    if (options.kinematic_viscosity <= 0.0) {
        return;
    }

    Boundary::fillVelocityPhysicalBoundary(state, options);

    amrex::MultiFab rhs(state.velocity().boxArray(),
                        state.velocity().DistributionMap(),
                        AMREX_SPACEDIM,
                        0,
                        amrex::MFInfo(),
                        state.velocity().Factory());
    amrex::MultiFab::Copy(rhs, state.velocity(), 0, 0, AMREX_SPACEDIM, 0);

    amrex::LPInfo lp_info;
    amrex::MLTensorOp tensor_op({state.geom()},
                                {state.velocity().boxArray()},
                                {state.velocity().DistributionMap()},
                                lp_info);
    auto [lo, hi] = Boundary::velocityOperatorBCs(state.geom(), options);
    tensor_op.setDomainBC(lo, hi);
    tensor_op.setLevelBC(0, &state.velocity());
    tensor_op.setScalars(1.0, dt);
    tensor_op.setACoeffs(0, 1.0);
    tensor_op.setShearViscosity(0, options.kinematic_viscosity);
    tensor_op.setBulkViscosity(0, 0.0);

    amrex::MLMG mlmg(tensor_op);
    mlmg.setVerbose(std::max(0, options.verbose - 1));
    mlmg.setFinalFillBC(true);
    mlmg.solve({&state.velocity()}, {&rhs},
               options.viscosity_rel_tol, options.viscosity_abs_tol);

    Boundary::fillVelocityPhysicalBoundary(state, options);

    if (options.viscous_reflux &&
        (crse_flux_register != nullptr || fine_flux_register != nullptr))
    {
        DLM::FaceFields viscous_fluxes;
        computeExtensiveViscousFluxes(mlmg, state, viscous_fluxes, amrex::Real(1.0));
        registerFluxes(viscous_fluxes, crse_flux_register, fine_flux_register, dt);
    }
}

} // namespace

void
Solver::applyViscosityLevel(int lev, amrex::Real dt)
{
    BL_PROFILE("ELFF::Fluid::Solver::applyViscosityLevel");
    AMREX_ALWAYS_ASSERT(lev >= 0);
    AMREX_ALWAYS_ASSERT(lev < static_cast<int>(levels_.size()));
    AMREX_ALWAYS_ASSERT(levels_[lev] != nullptr);

    LevelData& state = *levels_[lev];
    amrex::FluxRegister* crse_visc_flux_register =
        (options_.viscous_reflux &&
         lev + 1 < static_cast<int>(visc_flux_registers_.size()))
        ? visc_flux_registers_[lev + 1].get()
        : nullptr;
    amrex::FluxRegister* fine_visc_flux_register =
        (options_.viscous_reflux &&
         lev > 0 &&
         lev < static_cast<int>(visc_flux_registers_.size()))
        ? visc_flux_registers_[lev].get()
        : nullptr;
    applyImplicitViscosity(state,
                           options_,
                           crse_visc_flux_register,
                           fine_visc_flux_register,
                           dt);
}

void
Solver::applyViscosityHierarchy(amrex::Real dt)
{
    BL_PROFILE("ELFF::Fluid::Solver::applyViscosityHierarchy");
    if (options_.kinematic_viscosity <= 0.0) {
        return;
    }

    if (!options_.composite_viscosity || finestLevel() <= 0) {
        for (int lev = 0; lev <= finestLevel(); ++lev) {
            applyViscosityLevel(lev, dt);
        }
        return;
    }

    const int nlevels = finestLevel() + 1;
    amrex::Vector<amrex::Geometry> geoms;
    amrex::Vector<amrex::BoxArray> grids;
    amrex::Vector<amrex::DistributionMapping> dmaps;
    amrex::Vector<amrex::MultiFab*> velocities;
    amrex::Vector<amrex::MultiFab const*> rhs;
    std::vector<std::unique_ptr<amrex::MultiFab>> rhs_storage;
    geoms.reserve(nlevels);
    grids.reserve(nlevels);
    dmaps.reserve(nlevels);
    velocities.reserve(nlevels);
    rhs.reserve(nlevels);
    rhs_storage.reserve(nlevels);

    for (int lev = 0; lev < nlevels; ++lev) {
        AMREX_ALWAYS_ASSERT(lev < static_cast<int>(levels_.size()));
        AMREX_ALWAYS_ASSERT(levels_[lev] != nullptr);

        LevelData& state = *levels_[lev];
        Boundary::fillVelocityPhysicalBoundary(state, options_);
        geoms.push_back(Geom(lev));
        grids.push_back(state.velocity().boxArray());
        dmaps.push_back(state.velocity().DistributionMap());
        velocities.push_back(&state.velocity());

        auto rhs_level = std::make_unique<amrex::MultiFab>(state.velocity().boxArray(),
                                                           state.velocity().DistributionMap(),
                                                           AMREX_SPACEDIM,
                                                           0,
                                                           amrex::MFInfo(),
                                                           state.velocity().Factory());
        amrex::MultiFab::Copy(*rhs_level, state.velocity(), 0, 0, AMREX_SPACEDIM, 0);
        rhs.push_back(rhs_level.get());
        rhs_storage.push_back(std::move(rhs_level));
    }

    amrex::LPInfo lp_info;
    amrex::MLTensorOp tensor_op(geoms, grids, dmaps, lp_info);
    auto [lo, hi] = Boundary::velocityOperatorBCs(Geom(0), options_);
    tensor_op.setDomainBC(lo, hi);
    tensor_op.setScalars(1.0, dt);
    for (int lev = 0; lev < nlevels; ++lev) {
        tensor_op.setLevelBC(lev, &levels_[lev]->velocity());
        tensor_op.setACoeffs(lev, 1.0);
        tensor_op.setShearViscosity(lev, options_.kinematic_viscosity);
        tensor_op.setBulkViscosity(lev, 0.0);
    }

    amrex::MLMG mlmg(tensor_op);
    mlmg.setVerbose(std::max(0, options_.verbose - 1));
    mlmg.setFinalFillBC(true);
    mlmg.solve(velocities, rhs,
               options_.viscosity_rel_tol, options_.viscosity_abs_tol);

    for (int lev = 0; lev < nlevels; ++lev) {
        Boundary::fillVelocityPhysicalBoundary(*levels_[lev], options_);
    }

    if (options_.viscous_reflux && finestLevel() > 0) {
        std::vector<DLM::FaceFields> viscous_fluxes(nlevels);
        amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>> flux_ptrs;
        flux_ptrs.reserve(nlevels);
        for (int lev = 0; lev < nlevels; ++lev) {
            defineVectorFaceFluxesLike(*levels_[lev], viscous_fluxes[lev]);
            flux_ptrs.push_back(amrex::GetArrOfPtrs(viscous_fluxes[lev]));
        }

        mlmg.getFluxes(flux_ptrs,
                       velocities,
                       amrex::MLLinOp::Location::FaceCenter);

        for (int lev = 0; lev < nlevels; ++lev) {
            scaleFluxesByFaceArea(viscous_fluxes[lev], Geom(lev), amrex::Real(1.0));
            amrex::FluxRegister* crse_visc_flux_register =
                (lev + 1 < static_cast<int>(visc_flux_registers_.size()))
                ? visc_flux_registers_[lev + 1].get()
                : nullptr;
            amrex::FluxRegister* fine_visc_flux_register =
                (lev > 0 && lev < static_cast<int>(visc_flux_registers_.size()))
                ? visc_flux_registers_[lev].get()
                : nullptr;
            registerFluxes(viscous_fluxes[lev],
                           crse_visc_flux_register,
                           fine_visc_flux_register,
                           dt);
        }
    }
}

void
Solver::advanceLevelExplicit(int lev, amrex::Real time, amrex::Real dt)
{
    BL_PROFILE("ELFF::Fluid::Solver::advanceLevelExplicit");
    AMREX_ALWAYS_ASSERT(lev >= 0);
    AMREX_ALWAYS_ASSERT(lev < static_cast<int>(levels_.size()));
    AMREX_ALWAYS_ASSERT(levels_[lev] != nullptr);

    LevelData& state = *levels_[lev];

    if (options_.advection_type != "MOL") {
        amrex::Abort("Only fluid.advection_type = MOL is implemented in the AMReX-Hydro solver path");
    }

    fillPatchLevel(lev, time);

    auto velocity_bcs = Boundary::makeVelocityBCs(state.geom(), options_);
    auto device_velocity_bcs = Boundary::makeDeviceBCs(velocity_bcs);

    const DLMOptions dlm_options = dlmOptionsForCurrentStep();
    const bool staggered_primary =
        dlm_options.enabled &&
        dlm_options.multiplier_layout == "pacific" &&
        dlm_options.fluid_response == "mass";
    const bool retain_staggered_state = staggered_primary && finestLevel() == 0;
    if (!retain_staggered_state) {
        MOL::ExtrapVelToFaces(state.velocity(),
                              AMREX_D_DECL(state.faceVelocity()[0],
                                           state.faceVelocity()[1],
                                           state.faceVelocity()[2]),
                              state.geom(),
                              velocity_bcs,
                              device_velocity_bcs.dataPtr());

        projectFaceVelocity(state, options_);
    }

    amrex::MultiFab advective_term(state.velocity().boxArray(),
                                   state.velocity().DistributionMap(),
                                   AMREX_SPACEDIM,
                                   0,
                                   amrex::MFInfo(),
                                   state.velocity().Factory());
    amrex::FluxRegister* crse_flux_register =
        (options_.advective_reflux &&
         lev + 1 < static_cast<int>(adv_flux_registers_.size()))
        ? adv_flux_registers_[lev + 1].get()
        : nullptr;
    amrex::FluxRegister* fine_flux_register =
        (options_.advective_reflux &&
         lev > 0 &&
         lev < static_cast<int>(adv_flux_registers_.size()))
        ? adv_flux_registers_[lev].get()
        : nullptr;
    computeMOLAdvectiveTerm(state,
                            options_,
                            advective_term,
                            crse_flux_register,
                            fine_flux_register,
                            dt);
    applyExplicitVelocityUpdate(state, advective_term, dt);
    applyDrivingForces(state, options_, dt);
}

void
Solver::projectAndApplyConstraints(int lev, amrex::Real dt)
{
    BL_PROFILE("ELFF::Fluid::Solver::projectAndApplyConstraints");
    AMREX_ALWAYS_ASSERT(lev >= 0);
    AMREX_ALWAYS_ASSERT(lev < static_cast<int>(levels_.size()));
    AMREX_ALWAYS_ASSERT(levels_[lev] != nullptr);

    LevelData& state = *levels_[lev];
    const DLMOptions dlm_options = dlmOptionsForCurrentStep();
    const bool staggered_primary =
        dlm_options.enabled &&
        dlm_options.multiplier_layout == "pacific" &&
        dlm_options.fluid_response == "mass";
    const bool retain_staggered_state = staggered_primary && finestLevel() == 0;
    const bool collect_sync_projection_residuals =
        options_.sync_projection_residuals && !dlm_options.enabled;

    const amrex::Real level_profile_start = amrex::ParallelDescriptor::second();
    amrex::Real phase_start = level_profile_start;
    amrex::Real projection_time = 0.0;
    amrex::Real dlm_time = 0.0;
    amrex::Real face_refresh_time = 0.0;
    amrex::Real mac_projection_time = 0.0;

    std::unique_ptr<amrex::MultiFab> sync_resid_crse;
    std::unique_ptr<amrex::MultiFab> sync_resid_fine;
    if (collect_sync_projection_residuals) {
        if (lev < finestLevel()) {
            sync_resid_crse = std::make_unique<amrex::MultiFab>(
                state.pressure().boxArray(),
                state.pressure().DistributionMap(),
                1,
                state.pressure().nGrowVect(),
                amrex::MFInfo(),
                state.pressure().Factory());
            sync_resid_crse->setVal(0.0, 0, 1, sync_resid_crse->nGrowVect());
        }
        if (lev > 0) {
            sync_resid_fine = std::make_unique<amrex::MultiFab>(
                state.pressure().boxArray(),
                state.pressure().DistributionMap(),
                1,
                state.pressure().nGrowVect(),
                amrex::MFInfo(),
                state.pressure().Factory());
            sync_resid_fine->setVal(0.0, 0, 1, sync_resid_fine->nGrowVect());
        }
    }

    if (retain_staggered_state) {
        // Preserve the subcell staggered DLM state while carrying the
        // cell-centered advection/viscosity increment onto its normal faces.
        // Reconstructing the complete face state from cell averages here
        // applies a [1/4, 1/2, 1/4] filter every time step.
        Boundary::fillVelocityPhysicalBoundary(state, options_);
        addCellVelocityIncrementToFaces(state);
    }

    projectVelocityField(state.velocity(),
                         state.geom(),
                         options_,
                         &state.pressure(),
                         &state.gradp(),
                         sync_resid_crse.get(),
                         sync_resid_fine.get(),
                         lev < finestLevel()
                             ? refRatio(lev)
                             : amrex::IntVect(1),
                         lev < finestLevel()
                             ? levels_[lev + 1]->velocity().boxArray()
                             : amrex::BoxArray());
    projection_time = amrex::ParallelDescriptor::second() - phase_start;

    if (collect_sync_projection_residuals &&
        lev < static_cast<int>(sync_projection_residuals_.size()) &&
        sync_projection_residuals_[lev] != nullptr)
    {
        sync_projection_residuals_[lev]->setVal(0.0,
                                                0,
                                                1,
                                                sync_projection_residuals_[lev]->nGrowVect());
        if (sync_resid_crse != nullptr) {
            amrex::MultiFab::Add(*sync_projection_residuals_[lev],
                                 *sync_resid_crse,
                                 0,
                                 0,
                                 1,
                                 0);
        }
        if (sync_resid_fine != nullptr) {
            amrex::MultiFab::Add(*sync_projection_residuals_[lev],
                                 *sync_resid_fine,
                                 0,
                                 0,
                                 1,
                                 0);
        }
        sync_projection_residuals_[lev]->FillBoundary(Geom(lev).periodicity());
    }

    if (dlm_options.enabled) {
        const std::vector<DLM::RigidSphere> spheres = currentDlmSpheres();
        const bool apply_dlm_here =
            !dlm_options.apply_on_finest_level_only || lev == finestLevel();
        if (!apply_dlm_here) {
            phase_start = amrex::ParallelDescriptor::second();
            DLM::initialize(dlm_states_[lev], state, spheres, dlm_options);
            dlm_time = amrex::ParallelDescriptor::second() - phase_start;
            if (options_.profile_phases) {
                amrex::Print() << "ELFF_FLUID_PROFILE project_dlm_level step "
                               << step_ + 1
                               << " lev " << lev
                               << " total "
                               << amrex::ParallelDescriptor::second() - level_profile_start
                               << " projection " << projection_time
                               << " dlm_initialize " << dlm_time
                               << " face_refresh " << face_refresh_time
                               << " mac_projection " << mac_projection_time
                               << "\n";
            }
            return;
        }
        if (dlm_options.require_finest_stencil_support) {
            DLM::requireFullStencilSupport(state, spheres, dlm_options);
        }

        if (staggered_primary) {
            phase_start = amrex::ParallelDescriptor::second();
            if (retain_staggered_state) {
                projectFaceVelocity(state, options_);
            } else {
                refreshFaceVelocityFromCell(state, options_);
            }
            face_refresh_time += amrex::ParallelDescriptor::second() - phase_start;
        }

        StaggeredResponseWorkspace staggered_workspace;
        defineStaggeredResponseWorkspace(state, staggered_workspace);
        DLM::StaggeredFluidResponse fluid_response =
            [this, &state, dt, &staggered_workspace]
            (DLM::FaceFields const& rhs_delta, DLM::FaceFields& velocity_delta) {
                BL_PROFILE("ELFF::Fluid::DLM::fluidResponse");
                solveStaggeredVelocityResponse(state, options_, rhs_delta,
                                               velocity_delta, staggered_workspace, dt);
            };
        phase_start = amrex::ParallelDescriptor::second();
        DLM::applyRigidProjection(dlm_states_[lev],
                                  state,
                                  spheres,
                                  dlm_options,
                                  dt,
                                  fluid_response);
        dlm_time = amrex::ParallelDescriptor::second() - phase_start;

        phase_start = amrex::ParallelDescriptor::second();
        Boundary::fillVelocityPhysicalBoundary(state, options_);
        if (dlm_options.fluid_response == "mass") {
            Boundary::fillFaceVelocityPhysicalBoundary(
                state.faceVelocity(), state.geom(), options_);
            amrex::computeDivergence(
                state.divu(),
                amrex::GetArrOfConstPtrs(state.faceVelocity()),
                state.geom());
        } else {
            auto post_dlm_velocity_bcs = Boundary::makeVelocityBCs(state.geom(), options_);
            auto post_dlm_device_velocity_bcs = Boundary::makeDeviceBCs(
                post_dlm_velocity_bcs);
            MOL::ExtrapVelToFaces(state.velocity(),
                                  AMREX_D_DECL(state.faceVelocity()[0],
                                               state.faceVelocity()[1],
                                               state.faceVelocity()[2]),
                                  state.geom(),
                                  post_dlm_velocity_bcs,
                                  post_dlm_device_velocity_bcs.dataPtr());
        }
        face_refresh_time += amrex::ParallelDescriptor::second() - phase_start;
        if (dlm_options.fluid_response != "mass") {
            phase_start = amrex::ParallelDescriptor::second();
            projectFaceVelocity(state, options_);
            mac_projection_time = amrex::ParallelDescriptor::second() - phase_start;
        }
    }

    if (options_.profile_phases) {
        amrex::Print() << "ELFF_FLUID_PROFILE project_dlm_level step "
                       << step_ + 1
                       << " lev " << lev
                       << " total "
                       << amrex::ParallelDescriptor::second() - level_profile_start
                       << " projection " << projection_time
                       << " dlm_apply " << dlm_time
                       << " face_refresh " << face_refresh_time
                       << " mac_projection " << mac_projection_time
                       << "\n";
    }
}

void
Solver::advanceLevel(int lev, amrex::Real time, amrex::Real dt)
{
    advanceLevelExplicit(lev, time, dt);
    if (!options_.composite_viscosity) {
        applyViscosityLevel(lev, dt);
    }
    projectAndApplyConstraints(lev, dt);
}

} // namespace ELFF::Fluid
