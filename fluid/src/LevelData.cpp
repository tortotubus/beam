#include <elff/fluid/LevelData.H>

#include <AMReX_BLassert.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_MFIter.H>
#include <AMReX_Math.H>
#include <AMReX_ParmParse.H>
#include <AMReX_SPACE.H>

#include <string>
#include <utility>

namespace ELFF::Fluid {

void
LevelData::define(int lev,
                  amrex::Geometry geom,
                  amrex::BoxArray const& grids,
                  amrex::DistributionMapping const& dmap)
{
    level_ = lev;
    geom_ = std::move(geom);

    constexpr int ng_state = 2;
    constexpr int ng_pressure = 1;
    constexpr int ng_gradp = 1;

    velocity_.define(grids, dmap, AMREX_SPACEDIM, ng_state);
    density_.define(grids, dmap, 1, ng_state);
    gradp_.define(grids, dmap, AMREX_SPACEDIM, ng_gradp);
    old_velocity_.define(grids, dmap, AMREX_SPACEDIM, ng_state);
    old_density_.define(grids, dmap, 1, ng_state);
    old_gradp_.define(grids, dmap, AMREX_SPACEDIM, ng_gradp);
    divu_.define(grids, dmap, 1, 0);
    rigid_mask_.define(grids, dmap, 1, 0);
    rigid_velocity_.define(grids, dmap, AMREX_SPACEDIM, 0);
    dlm_force_.define(grids, dmap, AMREX_SPACEDIM, 0);

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        amrex::BoxArray face_grids(grids);
        face_grids.surroundingNodes(dir);
        face_velocity_[dir].define(face_grids, dmap, 1, 0);
        face_beta_[dir].define(face_grids, dmap, 1, 0);
    }

    amrex::BoxArray nodal_grids(grids);
    nodal_grids.surroundingNodes();
    pressure_.define(nodal_grids, dmap, 1, ng_pressure);
    old_pressure_.define(nodal_grids, dmap, 1, ng_pressure);

    defined_ = true;
}

void
LevelData::setInitialState(amrex::Real time)
{
    AMREX_ALWAYS_ASSERT(defined_);

    std::string init_type = "zero";
    amrex::Real velocity_scale = 1.0;
    amrex::Real rho0 = 1.0;
    {
        amrex::ParmParse pp("prob");
        pp.query("init_type", init_type);
        pp.query("velocity_scale", velocity_scale);
        pp.query("rho0", rho0);
    }

    velocity_.setVal(0.0);
    density_.setVal(rho0);
    pressure_.setVal(0.0);
    gradp_.setVal(0.0);
    old_velocity_.setVal(0.0);
    old_density_.setVal(rho0);
    old_pressure_.setVal(0.0);
    old_gradp_.setVal(0.0);
    divu_.setVal(0.0);
    rigid_mask_.setVal(0.0);
    rigid_velocity_.setVal(0.0);
    dlm_force_.setVal(0.0);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        face_velocity_[dir].setVal(0.0);
        face_beta_[dir].setVal(1.0);
    }

    if (init_type == "taylor_green") {
        const auto geom = geom_.data();
        for (amrex::MFIter mfi(velocity_, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.tilebox();
            auto const vel = velocity_.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const auto* dx = geom.CellSize();
                const auto* plo = geom.ProbLo();
                const amrex::Real x = plo[0] + (amrex::Real(i) + amrex::Real(0.5)) * dx[0];
                const amrex::Real y = plo[1] + (amrex::Real(j) + amrex::Real(0.5)) * dx[1];
#if (AMREX_SPACEDIM == 2)
                vel(i,j,k,0) = velocity_scale * amrex::Math::sinpi(amrex::Real(2.0) * x)
                                                * amrex::Math::cospi(amrex::Real(2.0) * y);
                vel(i,j,k,1) = -velocity_scale * amrex::Math::cospi(amrex::Real(2.0) * x)
                                                 * amrex::Math::sinpi(amrex::Real(2.0) * y);
#else
                const amrex::Real z = plo[2] + (amrex::Real(k) + amrex::Real(0.5)) * dx[2];
                vel(i,j,k,0) = velocity_scale * amrex::Math::sinpi(amrex::Real(2.0) * x)
                                                * amrex::Math::cospi(amrex::Real(2.0) * y)
                                                * amrex::Math::cospi(amrex::Real(2.0) * z);
                vel(i,j,k,1) = -velocity_scale * amrex::Math::cospi(amrex::Real(2.0) * x)
                                                 * amrex::Math::sinpi(amrex::Real(2.0) * y)
                                                 * amrex::Math::cospi(amrex::Real(2.0) * z);
                vel(i,j,k,2) = amrex::Real(0.0);
#endif
            });
        }
    } else if (init_type != "zero") {
        amrex::Abort("Unknown prob.init_type: " + init_type);
    }

    velocity_.FillBoundary(geom_.periodicity());
    snapshotOldState(time);
    new_time_ = time;
}

void
LevelData::snapshotOldState(amrex::Real time)
{
    AMREX_ALWAYS_ASSERT(defined_);

    amrex::MultiFab::Copy(old_velocity_, velocity_, 0, 0, AMREX_SPACEDIM,
                          old_velocity_.nGrowVect());
    amrex::MultiFab::Copy(old_density_, density_, 0, 0, 1,
                          old_density_.nGrowVect());
    amrex::MultiFab::Copy(old_pressure_, pressure_, 0, 0, 1,
                          old_pressure_.nGrowVect());
    amrex::MultiFab::Copy(old_gradp_, gradp_, 0, 0, AMREX_SPACEDIM,
                          old_gradp_.nGrowVect());
    old_time_ = time;
    new_time_ = time;
}

void
LevelData::copyStateFrom(LevelData const& src)
{
    AMREX_ALWAYS_ASSERT(defined_);
    AMREX_ALWAYS_ASSERT(src.isDefined());

    velocity_.ParallelCopy(src.velocity(), 0, 0, AMREX_SPACEDIM, 0, 0, geom_.periodicity());
    density_.ParallelCopy(src.density(), 0, 0, 1, 0, 0, geom_.periodicity());
    pressure_.ParallelCopy(src.pressure(), 0, 0, 1, 0, 0, geom_.periodicity());
    gradp_.ParallelCopy(src.gradp(), 0, 0, AMREX_SPACEDIM, 0, 0, geom_.periodicity());
    old_velocity_.ParallelCopy(src.oldVelocity(), 0, 0, AMREX_SPACEDIM, 0, 0,
                               geom_.periodicity());
    old_density_.ParallelCopy(src.oldDensity(), 0, 0, 1, 0, 0, geom_.periodicity());
    old_pressure_.ParallelCopy(src.oldPressure(), 0, 0, 1, 0, 0, geom_.periodicity());
    old_gradp_.ParallelCopy(src.oldGradp(), 0, 0, AMREX_SPACEDIM, 0, 0,
                            geom_.periodicity());
    divu_.ParallelCopy(src.divu(), 0, 0, 1, 0, 0, geom_.periodicity());
    rigid_mask_.ParallelCopy(src.rigidMask(), 0, 0, 1, 0, 0, geom_.periodicity());
    rigid_velocity_.ParallelCopy(src.rigidVelocity(), 0, 0, AMREX_SPACEDIM, 0, 0,
                                geom_.periodicity());
    dlm_force_.ParallelCopy(src.dlmForce(), 0, 0, AMREX_SPACEDIM, 0, 0,
                            geom_.periodicity());
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        face_velocity_[dir].ParallelCopy(src.faceVelocity()[dir], 0, 0, 1, 0, 0,
                                         geom_.periodicity());
        face_beta_[dir].ParallelCopy(src.faceBeta()[dir], 0, 0, 1, 0, 0,
                                     geom_.periodicity());
    }
    old_time_ = src.oldTime();
    new_time_ = src.newTime();
}

} // namespace ELFF::Fluid
