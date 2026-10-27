#include "DLMInternal.H"

#include <AMReX_GpuContainers.H>

namespace ELFF::Fluid::DLM::Detail {

namespace {

struct DeviceBody
{
    amrex::Real center[3]{0.0, 0.0, 0.0};
    amrex::Real velocity[3]{0.0, 0.0, 0.0};
    amrex::Real angular_velocity[3]{0.0, 0.0, 0.0};
    amrex::Real radius_sq = 0.0;
};

amrex::Gpu::DeviceVector<DeviceBody>
makeDeviceBodies(std::vector<BodyState> const& bodies)
{
    amrex::Vector<DeviceBody> host_bodies(bodies.size());
    for (std::size_t i = 0; i < bodies.size(); ++i) {
        for (int dir = 0; dir < 3; ++dir) {
            host_bodies[i].center[dir] = bodies[i].center[dir];
            host_bodies[i].velocity[dir] = bodies[i].velocity[dir];
            host_bodies[i].angular_velocity[dir] = bodies[i].angular_velocity[dir];
        }
        host_bodies[i].radius_sq = bodies[i].radius * bodies[i].radius;
    }

    amrex::Gpu::DeviceVector<DeviceBody> device_bodies(host_bodies.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                     host_bodies.begin(),
                     host_bodies.end(),
                     device_bodies.begin());
    return device_bodies;
}

} // namespace

void
fillRigidMask(LevelData& state, DLMOptions const& options)
{
    if (!options.enabled) {
        state.rigidMask().setVal(0.0);
        state.rigidVelocity().setVal(0.0);
        state.dlmForce().setVal(0.0);
        return;
    }
    BodyState body;
    initializeBodyState(body, options);
    fillRigidMask(state, std::vector<BodyState>{body}, options);
}

void
fillRigidMask(LevelData& state,
              std::vector<BodyState> const& bodies,
              DLMOptions const& options)
{
    state.rigidMask().setVal(0.0);
    state.rigidVelocity().setVal(0.0);
    state.dlmForce().setVal(0.0);

    if (!options.enabled || bodies.empty()) {
        return;
    }

    const auto geom = state.geom().data();
    auto device_bodies = makeDeviceBodies(bodies);
    DeviceBody const* body_data = device_bodies.dataPtr();
    const int nbodies = static_cast<int>(bodies.size());

    for (amrex::MFIter mfi(state.rigidMask(), amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto const mask = state.rigidMask().array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const auto* dx = geom.CellSize();
            const auto* plo = geom.ProbLo();
            const amrex::Real x = plo[0] + (amrex::Real(i) + amrex::Real(0.5)) * dx[0];
            const amrex::Real y = plo[1] + (amrex::Real(j) + amrex::Real(0.5)) * dx[1];
            amrex::Real z = 0.0;
#if (AMREX_SPACEDIM == 3)
            z = plo[2] + (amrex::Real(k) + amrex::Real(0.5)) * dx[2];
#endif
            amrex::Real value = 0.0;
            for (int body_index = 0; body_index < nbodies; ++body_index) {
                DeviceBody const& body = body_data[body_index];
                amrex::Real dist_sq = (x - body.center[0]) * (x - body.center[0])
                                    + (y - body.center[1]) * (y - body.center[1]);
#if (AMREX_SPACEDIM == 3)
                dist_sq += (z - body.center[2]) * (z - body.center[2]);
#else
                static_cast<void>(z);
#endif
                if (dist_sq <= body.radius_sq) {
                    value = 1.0;
                }
            }
            mask(i,j,k) = value;
        });
    }
}

void
fillRigidVelocityDiagnostic(LevelData& state,
                            RigidMotion const& motion,
                            DLMOptions const& options)
{
    BodyState body;
    initializeBodyState(body, options);
    setBodyMotion(body, motion);
    fillRigidVelocityDiagnostic(state, std::vector<BodyState>{body}, options);
}

void
fillRigidVelocityDiagnostic(LevelData& state,
                            std::vector<BodyState> const& bodies,
                            DLMOptions const& options)
{
    if (!options.enabled || bodies.empty()) {
        state.rigidVelocity().setVal(0.0);
        return;
    }

    const auto geom = state.geom().data();
    auto device_bodies = makeDeviceBodies(bodies);
    DeviceBody const* body_data = device_bodies.dataPtr();
    const int nbodies = static_cast<int>(bodies.size());

    for (amrex::MFIter mfi(state.rigidVelocity(), amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto const mask = state.rigidMask().const_array(mfi);
        auto const rigid_vel = state.rigidVelocity().array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (mask(i,j,k) <= 0.0) {
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                    rigid_vel(i,j,k,dir) = 0.0;
                }
                return;
            }

            const auto* dx = geom.CellSize();
            const auto* plo = geom.ProbLo();
            const amrex::Real x = plo[0] + (amrex::Real(i) + amrex::Real(0.5)) * dx[0];
            const amrex::Real y = plo[1] + (amrex::Real(j) + amrex::Real(0.5)) * dx[1];
            amrex::Real z = 0.0;
#if (AMREX_SPACEDIM == 3)
            z = plo[2] + (amrex::Real(k) + amrex::Real(0.5)) * dx[2];
#endif

            int owner = -1;
            for (int body_index = 0; body_index < nbodies; ++body_index) {
                DeviceBody const& body = body_data[body_index];
                amrex::Real dist_sq = (x - body.center[0]) * (x - body.center[0])
                                    + (y - body.center[1]) * (y - body.center[1]);
#if (AMREX_SPACEDIM == 3)
                dist_sq += (z - body.center[2]) * (z - body.center[2]);
#else
                static_cast<void>(z);
#endif
                if (dist_sq <= body.radius_sq) {
                    owner = body_index;
                    break;
                }
            }
            if (owner < 0) {
                return;
            }

            DeviceBody const& body = body_data[owner];
            const amrex::Real rx = x - body.center[0];
            const amrex::Real ry = y - body.center[1];

#if (AMREX_SPACEDIM == 2)
            rigid_vel(i,j,k,0) = body.velocity[0] - body.angular_velocity[2] * ry;
            rigid_vel(i,j,k,1) = body.velocity[1] + body.angular_velocity[2] * rx;
#else
            const amrex::Real rz = z - body.center[2];
            rigid_vel(i,j,k,0) = body.velocity[0]
                               + body.angular_velocity[1] * rz
                               - body.angular_velocity[2] * ry;
            rigid_vel(i,j,k,1) = body.velocity[1]
                               + body.angular_velocity[2] * rx
                               - body.angular_velocity[0] * rz;
            rigid_vel(i,j,k,2) = body.velocity[2]
                               + body.angular_velocity[0] * ry
                               - body.angular_velocity[1] * rx;
#endif
        });
    }
}

void
enforceInteriorRigidVelocityCells(LevelData& state,
                                  std::vector<BodyState> const& bodies,
                                  amrex::MultiFab& cumulative_delta)
{
    if (bodies.empty()) {
        return;
    }

    const auto geom = state.geom().data();
    auto device_bodies = makeDeviceBodies(bodies);
    DeviceBody const* body_data = device_bodies.dataPtr();
    const int nbodies = static_cast<int>(bodies.size());

    for (amrex::MFIter mfi(state.velocity(), amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto const mask = state.rigidMask().const_array(mfi);
        auto const vel = state.velocity().array(mfi);
        auto const cumulative = cumulative_delta.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (mask(i,j,k) <= 0.0) {
                return;
            }

            const auto* dx = geom.CellSize();
            const auto* plo = geom.ProbLo();
            const amrex::Real x = plo[0] + (amrex::Real(i) + amrex::Real(0.5)) * dx[0];
            const amrex::Real y = plo[1] + (amrex::Real(j) + amrex::Real(0.5)) * dx[1];
            amrex::Real z = 0.0;
#if (AMREX_SPACEDIM == 3)
            z = plo[2] + (amrex::Real(k) + amrex::Real(0.5)) * dx[2];
#endif

            int owner = -1;
            for (int body_index = 0; body_index < nbodies; ++body_index) {
                DeviceBody const& body = body_data[body_index];
                amrex::Real dist_sq = (x - body.center[0]) * (x - body.center[0])
                                    + (y - body.center[1]) * (y - body.center[1]);
#if (AMREX_SPACEDIM == 3)
                dist_sq += (z - body.center[2]) * (z - body.center[2]);
#else
                static_cast<void>(z);
#endif
                if (dist_sq <= body.radius_sq) {
                    owner = body_index;
                    break;
                }
            }
            if (owner < 0) {
                return;
            }

            DeviceBody const& body = body_data[owner];
            const amrex::Real rx = x - body.center[0];
            const amrex::Real ry = y - body.center[1];

#if (AMREX_SPACEDIM == 2)
            const amrex::Real rigid_x = body.velocity[0] - body.angular_velocity[2] * ry;
            const amrex::Real rigid_y = body.velocity[1] + body.angular_velocity[2] * rx;
            const amrex::Real old_x = vel(i,j,k,0);
            const amrex::Real old_y = vel(i,j,k,1);
            vel(i,j,k,0) = rigid_x;
            vel(i,j,k,1) = rigid_y;
            cumulative(i,j,k,0) += rigid_x - old_x;
            cumulative(i,j,k,1) += rigid_y - old_y;
#else
            const amrex::Real rz = z - body.center[2];
            const amrex::Real rigid_x = body.velocity[0]
                                      + body.angular_velocity[1] * rz
                                      - body.angular_velocity[2] * ry;
            const amrex::Real rigid_y = body.velocity[1]
                                      + body.angular_velocity[2] * rx
                                      - body.angular_velocity[0] * rz;
            const amrex::Real rigid_z = body.velocity[2]
                                      + body.angular_velocity[0] * ry
                                      - body.angular_velocity[1] * rx;
            const amrex::Real old_x = vel(i,j,k,0);
            const amrex::Real old_y = vel(i,j,k,1);
            const amrex::Real old_z = vel(i,j,k,2);
            vel(i,j,k,0) = rigid_x;
            vel(i,j,k,1) = rigid_y;
            vel(i,j,k,2) = rigid_z;
            cumulative(i,j,k,0) += rigid_x - old_x;
            cumulative(i,j,k,1) += rigid_y - old_y;
            cumulative(i,j,k,2) += rigid_z - old_z;
#endif
        });
    }

    state.velocity().FillBoundary(state.geom().periodicity());
}

void
setDlmForceFromDelta(LevelData& state, amrex::MultiFab const& cumulative_delta, amrex::Real dt)
{
    const amrex::Real inv_dt = (dt > 0.0) ? amrex::Real(1.0) / dt : amrex::Real(0.0);
    for (amrex::MFIter mfi(state.dlmForce(), amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto const force = state.dlmForce().array(mfi);
        auto const delta = cumulative_delta.const_array(mfi);
        amrex::ParallelFor(bx, AMREX_SPACEDIM,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
            force(i,j,k,n) = delta(i,j,k,n) * inv_dt;
        });
    }
}

void
addDlmForceFromDelta(LevelData& state, amrex::MultiFab const& cumulative_delta, amrex::Real dt)
{
    const amrex::Real inv_dt = (dt > 0.0) ? amrex::Real(1.0) / dt : amrex::Real(0.0);
    for (amrex::MFIter mfi(state.dlmForce(), amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto const force = state.dlmForce().array(mfi);
        auto const delta = cumulative_delta.const_array(mfi);
        amrex::ParallelFor(bx, AMREX_SPACEDIM,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
            force(i,j,k,n) += delta(i,j,k,n) * inv_dt;
        });
    }
}

void
setHydrodynamicLoads(State& dlm_state, amrex::Real dt)
{
    dlm_state.hydrodynamic_loads.clear();
    dlm_state.hydrodynamic_loads.reserve(dlm_state.bodies.size());
    for (BodyState& body : dlm_state.bodies) {
        const GeneralizedForce fluid_reaction =
            constraintGeneralizedForce(dlm_state.constraints,
                                       MarkerField::Lambda,
                                       body.body_index);
        HydrodynamicLoad load;
        load.id = body.id;
        const amrex::Real fluid_mass = body.fluid_density * sphereVolume(body.radius);
        const amrex::Real angular_inertia_factor =
            body.solid_density > 0.0 ? body.fluid_density / body.solid_density : 0.0;
        for (int dir = 0; dir < 3; ++dir) {
            load.multiplier_force[dir] = -fluid_reaction.force[dir];
            load.multiplier_torque[dir] = -fluid_reaction.torque[dir];
            load.fictitious_force[dir] = dt > 0.0
                ? fluid_mass * (body.velocity[dir] - body.previous_velocity[dir]) / dt
                : 0.0;
            load.fictitious_torque[dir] = dt > 0.0
                ? angular_inertia_factor * body.inertia[dir]
                    * (body.angular_velocity[dir] - body.previous_angular_velocity[dir]) / dt
                : 0.0;
            load.force[dir] = load.multiplier_force[dir] + load.fictitious_force[dir];
            load.torque[dir] = load.multiplier_torque[dir] + load.fictitious_torque[dir];
        }
        body.multiplier_force = load.multiplier_force;
        body.multiplier_torque = load.multiplier_torque;
        body.fictitious_force = load.fictitious_force;
        body.fictitious_torque = load.fictitious_torque;
        body.hydrodynamic_force = load.force;
        body.hydrodynamic_torque = load.torque;
        dlm_state.hydrodynamic_loads.push_back(load);
    }
    syncPrimaryBody(dlm_state);
}

} // namespace ELFF::Fluid::DLM::Detail
