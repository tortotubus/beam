#include "DLMInternal.H"

namespace ELFF::Fluid::DLM::Detail {

void
applyRigidVelocityCells(LevelData& state,
                        RigidMotion const& motion,
                        DLMOptions const& options,
                        amrex::Real dt)
{
    const auto geom = state.geom().data();
    const amrex::Real cx = options.center[0];
    const amrex::Real cy = options.center[1];
    const amrex::Real cz = options.center[2];
    const amrex::Real ux = motion.velocity[0];
    const amrex::Real uy = motion.velocity[1];
    const amrex::Real uz = motion.velocity[2];
    const amrex::Real wx = motion.angular_velocity[0];
    const amrex::Real wy = motion.angular_velocity[1];
    const amrex::Real wz = motion.angular_velocity[2];
    const amrex::Real inv_dt = (dt > 0.0) ? amrex::Real(1.0) / dt : amrex::Real(0.0);

    for (amrex::MFIter mfi(state.velocity(), amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto const vel = state.velocity().array(mfi);
        auto const mask = state.rigidMask().const_array(mfi);
        auto const rigid_vel = state.rigidVelocity().array(mfi);
        auto const force = state.dlmForce().array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const amrex::Real m = mask(i,j,k);
            if (m <= 0.0) {
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                    rigid_vel(i,j,k,dir) = 0.0;
                    force(i,j,k,dir) = 0.0;
                }
                return;
            }

            const auto* dx = geom.CellSize();
            const auto* plo = geom.ProbLo();
            const amrex::Real x = plo[0] + (amrex::Real(i) + amrex::Real(0.5)) * dx[0];
            const amrex::Real y = plo[1] + (amrex::Real(j) + amrex::Real(0.5)) * dx[1];
            const amrex::Real rx = x - cx;
            const amrex::Real ry = y - cy;

#if (AMREX_SPACEDIM == 2)
            const amrex::Real rvx = ux - wz * ry;
            const amrex::Real rvy = uy + wz * rx;
            const amrex::Real old_x = vel(i,j,k,0);
            const amrex::Real old_y = vel(i,j,k,1);
            vel(i,j,k,0) = rvx;
            vel(i,j,k,1) = rvy;
            rigid_vel(i,j,k,0) = rvx;
            rigid_vel(i,j,k,1) = rvy;
            force(i,j,k,0) = (rvx - old_x) * inv_dt;
            force(i,j,k,1) = (rvy - old_y) * inv_dt;
            static_cast<void>(cz);
            static_cast<void>(uz);
            static_cast<void>(wx);
            static_cast<void>(wy);
#else
            const amrex::Real z = plo[2] + (amrex::Real(k) + amrex::Real(0.5)) * dx[2];
            const amrex::Real rz = z - cz;
            const amrex::Real rvx = ux + wy * rz - wz * ry;
            const amrex::Real rvy = uy + wz * rx - wx * rz;
            const amrex::Real rvz = uz + wx * ry - wy * rx;
            const amrex::Real old_x = vel(i,j,k,0);
            const amrex::Real old_y = vel(i,j,k,1);
            const amrex::Real old_z = vel(i,j,k,2);
            vel(i,j,k,0) = rvx;
            vel(i,j,k,1) = rvy;
            vel(i,j,k,2) = rvz;
            rigid_vel(i,j,k,0) = rvx;
            rigid_vel(i,j,k,1) = rvy;
            rigid_vel(i,j,k,2) = rvz;
            force(i,j,k,0) = (rvx - old_x) * inv_dt;
            force(i,j,k,1) = (rvy - old_y) * inv_dt;
            force(i,j,k,2) = (rvz - old_z) * inv_dt;
#endif
        });
    }

    state.velocity().FillBoundary(state.geom().periodicity());
}

void
applyCellMaskProjection(LevelData& state, DLMOptions const& options, amrex::Real dt)
{
    const RigidMotion motion = rigidMotionCells(state, options);
    if (options.mode == "free" && motion.mass <= 0.0) {
        return;
    }
    applyRigidVelocityCells(state, motion, options, dt);
}

} // namespace ELFF::Fluid::DLM::Detail
