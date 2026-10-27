#include <elff/fluid/Boundary.H>

#include <AMReX_BC_TYPES.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_MFIter.H>

#include <algorithm>
#include <cctype>
#include <cstddef>

namespace ELFF::Fluid::Boundary {

std::string
canonicalBoundaryName(std::string name)
{
    std::string result;
    result.reserve(name.size());
    for (char ch : name) {
        if (ch == '-') {
            result.push_back('_');
        } else {
            result.push_back(static_cast<char>(
                std::tolower(static_cast<unsigned char>(ch))));
        }
    }
    return result;
}

BoundaryKind
parseBoundaryKind(std::string const& name)
{
    const std::string canonical = canonicalBoundaryName(name);
    if (canonical == "natural" || canonical == "interior") {
        return BoundaryKind::Natural;
    }
    if (canonical == "periodic") {
        return BoundaryKind::Periodic;
    }
    if (canonical == "inflow") {
        return BoundaryKind::Inflow;
    }
    if (canonical == "outflow") {
        return BoundaryKind::Outflow;
    }
    if (canonical == "no_slip" || canonical == "noslip" || canonical == "wall") {
        return BoundaryKind::NoSlip;
    }
    if (canonical == "slip" || canonical == "slip_wall" || canonical == "slipwall") {
        return BoundaryKind::Slip;
    }

    amrex::Abort("Unknown bc face type: " + name);
    return BoundaryKind::Natural;
}

namespace {

char const*
axisName(int dir)
{
    static constexpr std::array<char const*, 3> names{"x", "y", "z"};
    return names[static_cast<std::size_t>(dir)];
}

std::string
faceName(int dir, bool lower)
{
    return std::string("bc.") + axisName(dir) + (lower ? "lo" : "hi");
}

int
velocityBCType(BoundaryKind kind, int component, int normal_dir)
{
    switch (kind) {
    case BoundaryKind::Periodic:
        return amrex::BCType::int_dir;
    case BoundaryKind::Inflow:
    case BoundaryKind::NoSlip:
        return amrex::BCType::ext_dir;
    case BoundaryKind::Slip:
        return (component == normal_dir)
            ? amrex::BCType::ext_dir
            : amrex::BCType::hoextrap;
    case BoundaryKind::Outflow:
    case BoundaryKind::Natural:
        return amrex::BCType::foextrap;
    }
    return amrex::BCType::foextrap;
}

} // namespace

BoundaryFaceOptions const&
faceOptions(BoundaryOptions const& options, int dir, bool lower)
{
    return lower ? options.lo[dir] : options.hi[dir];
}

BoundaryKind
effectiveBoundaryKind(amrex::Geometry const& geom,
                      SolverOptions const& options,
                      int dir,
                      bool lower)
{
    const BoundaryFaceOptions& face = faceOptions(options.boundary, dir, lower);
    const BoundaryKind kind = parseBoundaryKind(face.type);
    if (geom.isPeriodic(dir)) {
        if (kind != BoundaryKind::Natural && kind != BoundaryKind::Periodic) {
            amrex::Abort(faceName(dir, lower)
                         + " cannot be non-periodic when geometry.is_periodic is true");
        }
        return BoundaryKind::Periodic;
    }
    if (kind == BoundaryKind::Periodic) {
        amrex::Abort(faceName(dir, lower)
                     + " = periodic requires geometry.is_periodic in that direction");
    }
    return kind;
}

amrex::Vector<amrex::BCRec>
makeVelocityBCs(amrex::Geometry const& geom, SolverOptions const& options)
{
    amrex::Vector<amrex::BCRec> bcs(AMREX_SPACEDIM);
    for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            const BoundaryKind lo_kind =
                effectiveBoundaryKind(geom, options, dir, true);
            const BoundaryKind hi_kind =
                effectiveBoundaryKind(geom, options, dir, false);
            bcs[comp].setLo(dir, velocityBCType(lo_kind, comp, dir));
            bcs[comp].setHi(dir, velocityBCType(hi_kind, comp, dir));
        }
    }
    return bcs;
}

amrex::Gpu::DeviceVector<amrex::BCRec>
makeDeviceBCs(amrex::Vector<amrex::BCRec> const& host_bcs)
{
    amrex::Gpu::DeviceVector<amrex::BCRec> device_bcs(host_bcs.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, host_bcs.begin(), host_bcs.end(),
                     device_bcs.begin());
    return device_bcs;
}

LinOpDomainBC
projectionBCs(amrex::Geometry const& geom, SolverOptions const& options)
{
    amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> lo;
    amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> hi;

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        auto projector_bc = [] (BoundaryKind kind) {
            switch (kind) {
            case BoundaryKind::Periodic:
                return amrex::LinOpBCType::Periodic;
            case BoundaryKind::Inflow:
                return amrex::LinOpBCType::Neumann;
            case BoundaryKind::Outflow:
                return amrex::LinOpBCType::Dirichlet;
            case BoundaryKind::NoSlip:
            case BoundaryKind::Slip:
            case BoundaryKind::Natural:
                return amrex::LinOpBCType::Neumann;
            }
            return amrex::LinOpBCType::Neumann;
        };
        lo[dir] = projector_bc(effectiveBoundaryKind(geom, options, dir, true));
        hi[dir] = projector_bc(effectiveBoundaryKind(geom, options, dir, false));
    }

    return {lo, hi};
}

LinOpDomainBC
velocityOperatorBCs(amrex::Geometry const& geom, SolverOptions const& options)
{
    amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> lo;
    amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> hi;

    auto velocity_bc = [] (BoundaryKind kind) {
        switch (kind) {
        case BoundaryKind::Periodic:
            return amrex::LinOpBCType::Periodic;
        case BoundaryKind::Inflow:
        case BoundaryKind::NoSlip:
            return amrex::LinOpBCType::Dirichlet;
        case BoundaryKind::Outflow:
        case BoundaryKind::Slip:
        case BoundaryKind::Natural:
            return amrex::LinOpBCType::Neumann;
        }
        return amrex::LinOpBCType::Neumann;
    };

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        lo[dir] = velocity_bc(effectiveBoundaryKind(geom, options, dir, true));
        hi[dir] = velocity_bc(effectiveBoundaryKind(geom, options, dir, false));
    }

    return {lo, hi};
}

void
fillVelocityPhysicalBoundary(amrex::MultiFab& velocity,
                             amrex::Geometry const& geom,
                             SolverOptions const& options)
{
    velocity.FillBoundary(geom.periodicity());

    const amrex::Box& domain = geom.Domain();
    for (int normal_dir = 0; normal_dir < AMREX_SPACEDIM; ++normal_dir) {
        for (int side = 0; side < 2; ++side) {
            const bool lower = (side == 0);
            const BoundaryKind kind =
                effectiveBoundaryKind(geom, options, normal_dir, lower);
            if (kind == BoundaryKind::Periodic) {
                continue;
            }

            const BoundaryFaceOptions& face =
                faceOptions(options.boundary, normal_dir, lower);
            amrex::GpuArray<amrex::Real, 3> face_velocity{
                face.velocity[0], face.velocity[1], face.velocity[2]};

            int fill_mode = 0;
            if (kind == BoundaryKind::Inflow || kind == BoundaryKind::NoSlip) {
                fill_mode = 1;
            } else if (kind == BoundaryKind::Slip) {
                fill_mode = 2;
            }

            const int interior_index = lower
                ? domain.smallEnd(normal_dir)
                : domain.bigEnd(normal_dir);
            const int ghost_limit = lower
                ? domain.smallEnd(normal_dir) - 1
                : domain.bigEnd(normal_dir) + 1;

            for (amrex::MFIter mfi(velocity, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                amrex::Box bx = velocity[mfi].box();
                if (lower) {
                    bx.setBig(normal_dir, std::min(bx.bigEnd(normal_dir), ghost_limit));
                } else {
                    bx.setSmall(normal_dir, std::max(bx.smallEnd(normal_dir), ghost_limit));
                }
                if (!bx.ok()) {
                    continue;
                }

                auto const vel = velocity.array(mfi);
                amrex::ParallelFor(bx, AMREX_SPACEDIM,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                    int ci = i;
                    int cj = j;
                    int ck = k;
                    if (normal_dir == 0) {
                        ci = interior_index;
                    } else if (normal_dir == 1) {
                        cj = interior_index;
#if (AMREX_SPACEDIM == 3)
                    } else {
                        ck = interior_index;
#endif
                    }

                    if (fill_mode == 1) {
                        vel(i,j,k,n) = face_velocity[n];
                    } else if (fill_mode == 2) {
                        vel(i,j,k,n) = (n == normal_dir)
                            ? face_velocity[n]
                            : vel(ci,cj,ck,n);
                    } else {
                        vel(i,j,k,n) = vel(ci,cj,ck,n);
                    }
                });
            }
        }
    }
}

void
fillVelocityPhysicalBoundary(LevelData& state, SolverOptions const& options)
{
    fillVelocityPhysicalBoundary(state.velocity(), state.geom(), options);
}

void
fillFaceVelocityPhysicalBoundary(DLM::FaceFields& face_velocity,
                                 amrex::Geometry const& geom,
                                 SolverOptions const& options)
{
    const amrex::Box& domain = geom.Domain();
    for (int normal_dir = 0; normal_dir < AMREX_SPACEDIM; ++normal_dir) {
        for (int side = 0; side < 2; ++side) {
            const bool lower = (side == 0);
            const BoundaryKind kind =
                effectiveBoundaryKind(geom, options, normal_dir, lower);
            if (kind == BoundaryKind::Periodic ||
                kind == BoundaryKind::Outflow ||
                kind == BoundaryKind::Natural)
            {
                continue;
            }

            const BoundaryFaceOptions& face =
                faceOptions(options.boundary, normal_dir, lower);
            const amrex::Real normal_velocity = face.velocity[normal_dir];
            const int face_index = lower
                ? domain.smallEnd(normal_dir)
                : domain.bigEnd(normal_dir) + 1;

            for (amrex::MFIter mfi(face_velocity[normal_dir], amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.validbox();
                bx.setSmall(normal_dir, face_index);
                bx.setBig(normal_dir, face_index);
                bx &= mfi.validbox();
                if (!bx.ok()) {
                    continue;
                }

                auto const vel = face_velocity[normal_dir].array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    vel(i,j,k) = normal_velocity;
                });
            }
        }
    }
}

} // namespace ELFF::Fluid::Boundary
