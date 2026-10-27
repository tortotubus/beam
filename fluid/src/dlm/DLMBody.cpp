#include "DLMInternal.H"

namespace ELFF::Fluid::DLM::Detail {

void
initializeBodyState(BodyState& body, DLMOptions const& options)
{
    initializeBodyState(body, optionSphere(options), options);
}

void
initializeBodyState(BodyState& body,
                    RigidSphere const& sphere,
                    DLMOptions const& options)
{
    if (options.fluid_density <= 0.0) {
        amrex::Abort("dlm.fluid_density must be positive");
    }
    const amrex::Real solid_density =
        (sphere.density > 0.0) ? sphere.density : options.solid_density;
    if (solid_density <= 0.0) {
        amrex::Abort("DLM sphere density must be positive");
    }
    if (sphere.radius <= 0.0) {
        amrex::Abort("DLM sphere radius must be positive");
    }

    body.initialized = true;
    body.id = sphere.id;
    body.geometry_id = sphere.geometry_id;
    body.particle = sphere.particle;
    body.prescribed = sphere.prescribed || options.mode == "prescribed";
    body.fluid_density = options.fluid_density;
    body.solid_density = solid_density;
    body.coupling_factor = options.explicit_added_mass
                         ? amrex::Real(1.0)
                         : amrex::Real(1.0) - options.fluid_density / solid_density;
    if (!body.prescribed && options.mode == "free" && std::abs(body.coupling_factor) <= tiny) {
        amrex::Abort("free DLM bodies require dlm.solid_density != dlm.fluid_density "
                     "unless dlm.explicit_added_mass = true");
    }

    body.radius = sphere.radius;
    body.center = sphere.center;
    body.mass = (sphere.mass > 0.0) ? sphere.mass : solid_density * sphereVolume(sphere.radius);
    body.velocity = sphere.velocity;
    body.angular_velocity = sphere.angular_velocity;
    body.previous_velocity = sphere.velocity;
    body.previous_angular_velocity = sphere.angular_velocity;
    body.multiplier_force = {0.0, 0.0, 0.0};
    body.multiplier_torque = {0.0, 0.0, 0.0};
    body.fictitious_force = {0.0, 0.0, 0.0};
    body.fictitious_torque = {0.0, 0.0, 0.0};
    body.hydrodynamic_force = {0.0, 0.0, 0.0};
    body.hydrodynamic_torque = {0.0, 0.0, 0.0};

#if (AMREX_SPACEDIM == 2)
    body.inertia = {0.0, 0.0, amrex::Real(0.5) * body.mass * sphere.radius * sphere.radius};
#else
    if (sphere.inertia[0] > 0.0 && sphere.inertia[1] > 0.0 && sphere.inertia[2] > 0.0) {
        body.inertia = sphere.inertia;
    } else {
        const amrex::Real isotropic_inertia =
            amrex::Real(0.4) * body.mass * sphere.radius * sphere.radius;
        body.inertia = {isotropic_inertia, isotropic_inertia, isotropic_inertia};
    }
#endif
}

std::vector<BodyState>
buildBodyStates(std::vector<RigidSphere> const& spheres, DLMOptions const& options)
{
    std::vector<BodyState> bodies;
    bodies.reserve(spheres.size());
    for (std::size_t body_index = 0; body_index < spheres.size(); ++body_index) {
        BodyState body;
        initializeBodyState(body, spheres[body_index], options);
        body.body_index = body_index;
        bodies.push_back(body);
    }
    return bodies;
}

void
ensureBodyState(BodyState& body, DLMOptions const& options)
{
    if (!body.initialized) {
        initializeBodyState(body, options);
    }
}

void
syncPrimaryBody(State& dlm_state)
{
    if (!dlm_state.bodies.empty()) {
        dlm_state.body = dlm_state.bodies.front();
    } else {
        dlm_state.body.clear();
    }
}

RigidMotion
bodyMotion(BodyState const& body)
{
    RigidMotion motion;
    motion.mass = body.mass;
    motion.velocity = body.velocity;
    motion.angular_velocity = body.angular_velocity;
    return motion;
}

RigidMotion
prescribedRigidMotion(DLMOptions const& options)
{
    RigidMotion motion;
    motion.velocity = options.velocity;
    motion.angular_velocity = options.angular_velocity;
    return motion;
}

RigidMotion
prescribedRigidMotion(BodyState const& body)
{
    return bodyMotion(body);
}

RigidMotion
bestFitRigidMotionCells(LevelData const& state, DLMOptions const& options)
{
    RigidMotion motion;
    const amrex::Real vol = cellVolume(state.geom());

    {
        amrex::ReduceOps<amrex::ReduceOpSum,
                         amrex::ReduceOpSum,
                         amrex::ReduceOpSum,
                         amrex::ReduceOpSum> reduce_op;
        amrex::ReduceData<amrex::Real,
                          amrex::Real,
                          amrex::Real,
                          amrex::Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (amrex::MFIter mfi(state.velocity(), amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.tilebox();
            auto const vel = state.velocity().const_array(mfi);
            auto const mask = state.rigidMask().const_array(mfi);

            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple {
                const amrex::Real weight = vol * mask(i,j,k);
                amrex::Real mom_z = 0.0;
#if (AMREX_SPACEDIM == 3)
                mom_z = weight * vel(i,j,k,2);
#else
                static_cast<void>(k);
#endif
                return {weight,
                        weight * vel(i,j,k,0),
                        weight * vel(i,j,k,1),
                        mom_z};
            });
        }

        const ReduceTuple hv = reduce_data.value(reduce_op);
        std::array<amrex::Real, 4> sums{
            amrex::get<0>(hv),
            amrex::get<1>(hv),
            amrex::get<2>(hv),
            amrex::get<3>(hv)
        };
        amrex::ParallelDescriptor::ReduceRealSum(sums.data(), static_cast<int>(sums.size()));

        motion.mass = sums[0];
        if (motion.mass <= 0.0) {
            return motion;
        }

        motion.velocity[0] = sums[1] / motion.mass;
        motion.velocity[1] = sums[2] / motion.mass;
        motion.velocity[2] = sums[3] / motion.mass;
    }

    {
        const auto geom = state.geom().data();
        const amrex::Real cx = options.center[0];
        const amrex::Real cy = options.center[1];
        const amrex::Real cz = options.center[2];
        const amrex::Real ux = motion.velocity[0];
        const amrex::Real uy = motion.velocity[1];
        const amrex::Real uz = motion.velocity[2];
#if (AMREX_SPACEDIM == 2)
        static_cast<void>(cz);
        static_cast<void>(uz);
#endif

        amrex::ReduceOps<amrex::ReduceOpSum,
                         amrex::ReduceOpSum,
                         amrex::ReduceOpSum,
                         amrex::ReduceOpSum,
                         amrex::ReduceOpSum,
                         amrex::ReduceOpSum,
                         amrex::ReduceOpSum,
                         amrex::ReduceOpSum,
                         amrex::ReduceOpSum> reduce_op;
        amrex::ReduceData<amrex::Real,
                          amrex::Real,
                          amrex::Real,
                          amrex::Real,
                          amrex::Real,
                          amrex::Real,
                          amrex::Real,
                          amrex::Real,
                          amrex::Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (amrex::MFIter mfi(state.velocity(), amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.tilebox();
            auto const vel = state.velocity().const_array(mfi);
            auto const mask = state.rigidMask().const_array(mfi);

            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple {
                const auto* dx = geom.CellSize();
                const auto* plo = geom.ProbLo();
                const amrex::Real x = plo[0] + (amrex::Real(i) + amrex::Real(0.5)) * dx[0];
                const amrex::Real y = plo[1] + (amrex::Real(j) + amrex::Real(0.5)) * dx[1];
                const amrex::Real rx = x - cx;
                const amrex::Real ry = y - cy;
                const amrex::Real qx = vel(i,j,k,0) - ux;
                const amrex::Real qy = vel(i,j,k,1) - uy;
                const amrex::Real weight = vol * mask(i,j,k);

#if (AMREX_SPACEDIM == 2)
                const amrex::Real lz = weight * (rx * qy - ry * qx);
                const amrex::Real izz = weight * (rx * rx + ry * ry);
                return {amrex::Real(0.0), amrex::Real(0.0), lz,
                        amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0),
                        amrex::Real(0.0), amrex::Real(0.0), izz};
#else
                const amrex::Real z = plo[2] + (amrex::Real(k) + amrex::Real(0.5)) * dx[2];
                const amrex::Real rz = z - cz;
                const amrex::Real qz = vel(i,j,k,2) - uz;

                const amrex::Real lx = weight * (ry * qz - rz * qy);
                const amrex::Real ly = weight * (rz * qx - rx * qz);
                const amrex::Real lz = weight * (rx * qy - ry * qx);
                const amrex::Real ixx = weight * (ry * ry + rz * rz);
                const amrex::Real ixy = -weight * rx * ry;
                const amrex::Real ixz = -weight * rx * rz;
                const amrex::Real iyy = weight * (rx * rx + rz * rz);
                const amrex::Real iyz = -weight * ry * rz;
                const amrex::Real izz = weight * (rx * rx + ry * ry);
                return {lx, ly, lz, ixx, ixy, ixz, iyy, iyz, izz};
#endif
            });
        }

        const ReduceTuple hv = reduce_data.value(reduce_op);
        std::array<amrex::Real, 9> sums{
            amrex::get<0>(hv),
            amrex::get<1>(hv),
            amrex::get<2>(hv),
            amrex::get<3>(hv),
            amrex::get<4>(hv),
            amrex::get<5>(hv),
            amrex::get<6>(hv),
            amrex::get<7>(hv),
            amrex::get<8>(hv)
        };
        amrex::ParallelDescriptor::ReduceRealSum(sums.data(), static_cast<int>(sums.size()));

#if (AMREX_SPACEDIM == 2)
        if (sums[8] > 0.0) {
            motion.angular_velocity[2] = sums[2] / sums[8];
        }
#else
        const amrex::Real ixx = sums[3];
        const amrex::Real ixy = sums[4];
        const amrex::Real ixz = sums[5];
        const amrex::Real iyy = sums[6];
        const amrex::Real iyz = sums[7];
        const amrex::Real izz = sums[8];
        const amrex::Real lx = sums[0];
        const amrex::Real ly = sums[1];
        const amrex::Real lz = sums[2];

        const amrex::Real det = ixx * (iyy * izz - iyz * iyz)
                              - ixy * (ixy * izz - iyz * ixz)
                              + ixz * (ixy * iyz - iyy * ixz);
        if (std::abs(det) > 0.0) {
            const amrex::Real inv00 = (iyy * izz - iyz * iyz) / det;
            const amrex::Real inv01 = (ixz * iyz - ixy * izz) / det;
            const amrex::Real inv02 = (ixy * iyz - ixz * iyy) / det;
            const amrex::Real inv11 = (ixx * izz - ixz * ixz) / det;
            const amrex::Real inv12 = (ixy * ixz - ixx * iyz) / det;
            const amrex::Real inv22 = (ixx * iyy - ixy * ixy) / det;
            motion.angular_velocity[0] = inv00 * lx + inv01 * ly + inv02 * lz;
            motion.angular_velocity[1] = inv01 * lx + inv11 * ly + inv12 * lz;
            motion.angular_velocity[2] = inv02 * lx + inv12 * ly + inv22 * lz;
        }
#endif
    }

    return motion;
}

RigidMotion
bestFitRigidMotionMarkers(std::vector<MarkerPoint> const& markers)
{
    return bestFitRigidMotionMarkers(markers, static_cast<std::size_t>(-1));
}

RigidMotion
bestFitRigidMotionMarkers(std::vector<MarkerPoint> const& markers,
                          std::size_t body_index)
{
    RigidMotion motion;
    for (MarkerPoint const& marker : markers) {
        if (body_index != static_cast<std::size_t>(-1) && marker.body_index != body_index) {
            continue;
        }
        const amrex::Real weight = marker.quadrature_weight;
        motion.mass += weight;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            motion.velocity[dir] += weight * marker.interpolated_velocity[dir];
        }
    }
    if (motion.mass <= 0.0) {
        return motion;
    }
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        motion.velocity[dir] /= motion.mass;
    }

    std::array<amrex::Real, 9> sums{};
    for (MarkerPoint const& marker : markers) {
        if (body_index != static_cast<std::size_t>(-1) && marker.body_index != body_index) {
            continue;
        }
        const amrex::Real weight = marker.quadrature_weight;
        const amrex::Real rx = marker.center_vector[0];
        const amrex::Real ry = marker.center_vector[1];
        const amrex::Real qx = marker.interpolated_velocity[0] - motion.velocity[0];
        const amrex::Real qy = marker.interpolated_velocity[1] - motion.velocity[1];
#if (AMREX_SPACEDIM == 2)
        sums[2] += weight * (rx * qy - ry * qx);
        sums[8] += weight * (rx * rx + ry * ry);
#else
        const amrex::Real rz = marker.center_vector[2];
        const amrex::Real qz = marker.interpolated_velocity[2] - motion.velocity[2];
        sums[0] += weight * (ry * qz - rz * qy);
        sums[1] += weight * (rz * qx - rx * qz);
        sums[2] += weight * (rx * qy - ry * qx);
        sums[3] += weight * (ry * ry + rz * rz);
        sums[4] += -weight * rx * ry;
        sums[5] += -weight * rx * rz;
        sums[6] += weight * (rx * rx + rz * rz);
        sums[7] += -weight * ry * rz;
        sums[8] += weight * (rx * rx + ry * ry);
#endif
    }

#if (AMREX_SPACEDIM == 2)
    if (sums[8] > 0.0) {
        motion.angular_velocity[2] = sums[2] / sums[8];
    }
#else
    const amrex::Real ixx = sums[3];
    const amrex::Real ixy = sums[4];
    const amrex::Real ixz = sums[5];
    const amrex::Real iyy = sums[6];
    const amrex::Real iyz = sums[7];
    const amrex::Real izz = sums[8];
    const amrex::Real lx = sums[0];
    const amrex::Real ly = sums[1];
    const amrex::Real lz = sums[2];

    const amrex::Real det = ixx * (iyy * izz - iyz * iyz)
                          - ixy * (ixy * izz - iyz * ixz)
                          + ixz * (ixy * iyz - iyy * ixz);
    if (std::abs(det) > 0.0) {
        const amrex::Real inv00 = (iyy * izz - iyz * iyz) / det;
        const amrex::Real inv01 = (ixz * iyz - ixy * izz) / det;
        const amrex::Real inv02 = (ixy * iyz - ixz * iyy) / det;
        const amrex::Real inv11 = (ixx * izz - ixz * ixz) / det;
        const amrex::Real inv12 = (ixy * ixz - ixx * iyz) / det;
        const amrex::Real inv22 = (ixx * iyy - ixy * ixy) / det;
        motion.angular_velocity[0] = inv00 * lx + inv01 * ly + inv02 * lz;
        motion.angular_velocity[1] = inv01 * lx + inv11 * ly + inv12 * lz;
        motion.angular_velocity[2] = inv02 * lx + inv12 * ly + inv22 * lz;
    }
#endif
    return motion;
}

RigidMotion
rigidMotionCells(LevelData const& state, DLMOptions const& options)
{
    if (options.mode == "prescribed") {
        return prescribedRigidMotion(options);
    }
    if (options.mode == "free") {
        return bestFitRigidMotionCells(state, options);
    }
    amrex::Abort("dlm.mode must be either 'free' or 'prescribed'");
    return {};
}

RigidMotion
rigidMotionMarkers(std::vector<MarkerPoint> const& markers, DLMOptions const& options)
{
    if (options.mode == "prescribed") {
        return prescribedRigidMotion(options);
    }
    if (options.mode == "free") {
        return bestFitRigidMotionMarkers(markers);
    }
    amrex::Abort("dlm.mode must be either 'free' or 'prescribed'");
    return {};
}

RigidMotion
rigidMotionMarkers(std::vector<MarkerPoint> const& markers,
                   BodyState const& body,
                   DLMOptions const& options)
{
    if (body.prescribed || options.mode == "prescribed") {
        return prescribedRigidMotion(body);
    }
    if (options.mode == "free") {
        return bestFitRigidMotionMarkers(markers, body.body_index);
    }
    amrex::Abort("dlm.mode must be either 'free' or 'prescribed'");
    return {};
}

std::array<amrex::Real, 3>
rigidVelocityAt(std::array<amrex::Real, 3> const& position,
                RigidMotion const& motion,
                DLMOptions const& options)
{
    const amrex::Real rx = position[0] - options.center[0];
    const amrex::Real ry = position[1] - options.center[1];
#if (AMREX_SPACEDIM == 2)
    return {motion.velocity[0] - motion.angular_velocity[2] * ry,
            motion.velocity[1] + motion.angular_velocity[2] * rx,
            0.0};
#else
    const amrex::Real rz = position[2] - options.center[2];
    return {motion.velocity[0] + motion.angular_velocity[1] * rz
                               - motion.angular_velocity[2] * ry,
            motion.velocity[1] + motion.angular_velocity[2] * rx
                               - motion.angular_velocity[0] * rz,
            motion.velocity[2] + motion.angular_velocity[0] * ry
                               - motion.angular_velocity[1] * rx};
#endif
}

std::array<amrex::Real, 3>
rigidVelocityAt(std::array<amrex::Real, 3> const& position,
                RigidMotion const& motion,
                BodyState const& body)
{
    const amrex::Real rx = position[0] - body.center[0];
    const amrex::Real ry = position[1] - body.center[1];
#if (AMREX_SPACEDIM == 2)
    return {motion.velocity[0] - motion.angular_velocity[2] * ry,
            motion.velocity[1] + motion.angular_velocity[2] * rx,
            0.0};
#else
    const amrex::Real rz = position[2] - body.center[2];
    return {motion.velocity[0] + motion.angular_velocity[1] * rz
                               - motion.angular_velocity[2] * ry,
            motion.velocity[1] + motion.angular_velocity[2] * rx
                               - motion.angular_velocity[0] * rz,
            motion.velocity[2] + motion.angular_velocity[0] * ry
                               - motion.angular_velocity[1] * rx};
#endif
}

GeneralizedForce
initialBodyRhs(BodyState const& body,
               DLMOptions const& options,
               GeneralizedForce marker_rhs,
               amrex::Real dt)
{
    GeneralizedForce rhs;
    const amrex::Real inv_dt = (dt > 0.0) ? amrex::Real(1.0) / dt : amrex::Real(0.0);
    const amrex::Real buoyant_gravity_factor =
        amrex::Real(1.0) - body.fluid_density / body.solid_density;

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        rhs.force[dir] = body.coupling_factor * body.mass * body.velocity[dir] * inv_dt
                       + buoyant_gravity_factor * body.mass * options.gravity[dir]
                       - marker_rhs.force[dir];
        rhs.torque[dir] = body.coupling_factor * body.inertia[dir]
                        * body.angular_velocity[dir] * inv_dt
                        - marker_rhs.torque[dir];
    }
#if (AMREX_SPACEDIM == 2)
    rhs.torque[2] = body.coupling_factor * body.inertia[2]
                  * body.angular_velocity[2] * inv_dt
                  - marker_rhs.torque[2];
#endif
    return rhs;
}

RigidMotion
solveBodyVelocity(BodyState const& body,
                  GeneralizedForce const& rhs,
                  amrex::Real dt)
{
    RigidMotion motion;
    motion.mass = body.mass;

    const amrex::Real force_denom = body.coupling_factor * body.mass;
    if (std::abs(force_denom) > tiny) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            motion.velocity[dir] = dt * rhs.force[dir] / force_denom;
        }
    }

    for (int dir = 0; dir < 3; ++dir) {
#if (AMREX_SPACEDIM == 2)
        if (dir != 2) {
            continue;
        }
#endif
        const amrex::Real torque_denom = body.coupling_factor * body.inertia[dir];
        if (std::abs(torque_denom) > tiny) {
            motion.angular_velocity[dir] = dt * rhs.torque[dir] / torque_denom;
        }
    }

    return motion;
}

RigidMotion
solveInitialBodyMotion(BodyState const& body,
                       DLMOptions const& options,
                       std::vector<ConstraintPoint> const& constraints,
                       amrex::Real dt)
{
    if (body.prescribed || options.mode == "prescribed") {
        return prescribedRigidMotion(body);
    }
    if (options.mode != "free") {
        amrex::Abort("dlm.mode must be either 'free' or 'prescribed'");
    }

    const GeneralizedForce lambda_rhs =
        constraintGeneralizedForce(constraints, MarkerField::Lambda,
                                   body.body_index);
    return solveBodyVelocity(body, initialBodyRhs(body, options, lambda_rhs, dt), dt);
}

RigidMotion
solveSearchBodyMotion(BodyState const& body,
                      DLMOptions const& options,
                      std::vector<ConstraintPoint> const& constraints,
                      amrex::Real dt)
{
    if (body.prescribed || options.mode == "prescribed") {
        return {};
    }
    if (options.mode != "free") {
        amrex::Abort("dlm.mode must be either 'free' or 'prescribed'");
    }

    const GeneralizedForce search_rhs =
        constraintGeneralizedForce(constraints, MarkerField::Search,
                                   body.body_index);
    return solveBodyVelocity(body, search_rhs, dt);
}

void
setBodyMotion(BodyState& body, RigidMotion const& motion)
{
    body.velocity = motion.velocity;
    body.angular_velocity = motion.angular_velocity;
}

void
addBodyMotion(BodyState& body, RigidMotion const& motion, amrex::Real scale)
{
    for (int dir = 0; dir < 3; ++dir) {
        body.velocity[dir] += scale * motion.velocity[dir];
        body.angular_velocity[dir] += scale * motion.angular_velocity[dir];
    }
}

} // namespace ELFF::Fluid::DLM::Detail
