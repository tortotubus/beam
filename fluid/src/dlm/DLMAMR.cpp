#include "DLMInternal.H"

#include <AMReX_IndexType.H>

#include <sstream>
#include <string>

namespace ELFF::Fluid::DLM {

namespace {

std::string
formatIntVect(amrex::IntVect const& iv)
{
    std::ostringstream os;
    os << "(" << iv[0] << ", " << iv[1];
#if (AMREX_SPACEDIM == 3)
    os << ", " << iv[2];
#endif
    os << ")";
    return os.str();
}

amrex::Box
cellPointBox(amrex::IntVect const& iv)
{
    return amrex::Box(iv, iv);
}

amrex::Box
facePointBox(amrex::IntVect const& iv, int component)
{
    amrex::IndexType type = amrex::IndexType::TheCellType();
    type.setType(component, amrex::IndexType::NODE);
    return amrex::Box(iv, iv, type);
}

amrex::IntVect
domainCell(amrex::Geometry const& geom, amrex::IntVect iv)
{
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        iv[dir] = Detail::domainIndex(geom, dir, iv[dir]);
    }
    return iv;
}

void
tagCell(amrex::TagBoxArray& tags,
        amrex::Geometry const& geom,
        amrex::IntVect const& cell)
{
    tags.setVal(amrex::TagBox::SET, cellPointBox(domainCell(geom, cell)), 0);
}

void
tagFaceAdjacentCells(amrex::TagBoxArray& tags,
                     amrex::Geometry const& geom,
                     int component,
                     amrex::IntVect const& face)
{
    amrex::IntVect lo = face;
    lo[component] -= 1;
    tagCell(tags, geom, lo);
    tagCell(tags, geom, face);
}

void
tagRigidMask(amrex::TagBoxArray& tags, amrex::MultiFab const& rigid_mask)
{
    const char tagval = amrex::TagBox::SET;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(rigid_mask, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto const mask = rigid_mask.const_array(mfi);
        auto const tag = tags.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (mask(i,j,k) > amrex::Real(0.5)) {
                tag(i,j,k) = tagval;
            }
        });
    }
}

void
tagSphereVolumes(amrex::TagBoxArray& tags,
                 LevelData const& state,
                 std::vector<RigidSphere> const& spheres,
                 DLMOptions const& options)
{
    if (spheres.empty()) {
        return;
    }

    amrex::GpuArray<int, 3> periodic{0, 0, 0};
    amrex::GpuArray<amrex::Real, 3> domain_length{0.0, 0.0, 0.0};
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        periodic[dir] = state.geom().isPeriodic(dir) ? 1 : 0;
        domain_length[dir] = state.geom().ProbHi(dir) - state.geom().ProbLo(dir);
    }

    const amrex::Real buffer = std::max(options.refinement_buffer, amrex::Real(0.0));
    const char tagval = amrex::TagBox::SET;
    const auto geom = state.geom().data();

    for (RigidSphere const& sphere : spheres) {
        const amrex::Real inflated_radius = sphere.radius + buffer;
        if (inflated_radius <= 0.0) {
            continue;
        }
        const amrex::Real radius_sq = inflated_radius * inflated_radius;
        const amrex::Real center0 = sphere.center[0];
        const amrex::Real center1 = sphere.center[1];
        const amrex::Real center2 = sphere.center[2];

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(tags, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.tilebox();
            auto const tag = tags.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const auto* dx = geom.CellSize();
                const auto* plo = geom.ProbLo();
                amrex::Real dist_sq = 0.0;

                amrex::Real x0 = plo[0] + (amrex::Real(i) + amrex::Real(0.5)) * dx[0];
                amrex::Real d0 = x0 - center0;
                if (periodic[0] != 0 && domain_length[0] > 0.0) {
                    if (d0 > amrex::Real(0.5) * domain_length[0]) {
                        d0 -= domain_length[0];
                    } else if (d0 < -amrex::Real(0.5) * domain_length[0]) {
                        d0 += domain_length[0];
                    }
                }
                dist_sq += d0 * d0;

                amrex::Real x1 = plo[1] + (amrex::Real(j) + amrex::Real(0.5)) * dx[1];
                amrex::Real d1 = x1 - center1;
                if (periodic[1] != 0 && domain_length[1] > 0.0) {
                    if (d1 > amrex::Real(0.5) * domain_length[1]) {
                        d1 -= domain_length[1];
                    } else if (d1 < -amrex::Real(0.5) * domain_length[1]) {
                        d1 += domain_length[1];
                    }
                }
                dist_sq += d1 * d1;

#if (AMREX_SPACEDIM == 3)
                amrex::Real x2 = plo[2] + (amrex::Real(k) + amrex::Real(0.5)) * dx[2];
                amrex::Real d2 = x2 - center2;
                if (periodic[2] != 0 && domain_length[2] > 0.0) {
                    if (d2 > amrex::Real(0.5) * domain_length[2]) {
                        d2 -= domain_length[2];
                    } else if (d2 < -amrex::Real(0.5) * domain_length[2]) {
                        d2 += domain_length[2];
                    }
                }
                dist_sq += d2 * d2;
#else
                static_cast<void>(k);
                static_cast<void>(center2);
#endif

                if (dist_sq <= radius_sq) {
                    tag(i,j,k) = tagval;
                }
            });
        }
    }
}

void
tagMarkerStencilSupport(amrex::TagBoxArray& tags,
                        LevelData const& state,
                        std::vector<RigidSphere> const& spheres,
                        DLMOptions const& options)
{
    if (spheres.empty()) {
        return;
    }

    std::vector<MarkerPoint> markers = Detail::buildMarkers(state.geom(), spheres, options);
    for (MarkerPoint const& marker : markers) {
        for (StencilEntry const& entry : marker.stencil) {
            tagCell(tags, state.geom(), entry.cell);
        }
        for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
            ConstraintPoint constraint = Detail::makeConstraint(state.geom(), options, marker, comp);
            for (StencilEntry const& entry : constraint.stencil) {
                tagFaceAdjacentCells(tags, state.geom(), comp, entry.cell);
            }
        }
    }
}

bool
hasMarkerCellSupport(LevelData const& state,
                     MarkerPoint const& marker,
                     std::string* reason)
{
    amrex::BoxArray const& cells = state.velocity().boxArray();
    for (StencilEntry const& entry : marker.stencil) {
        if (!cells.contains(cellPointBox(entry.cell))) {
            if (reason != nullptr) {
                *reason = "missing cell-centered DLM marker support at "
                        + formatIntVect(entry.cell);
            }
            return false;
        }
    }
    return true;
}

bool
hasConstraintFaceSupport(LevelData const& state,
                         ConstraintPoint const& constraint,
                         std::string* reason)
{
    amrex::BoxArray const& faces = state.faceVelocity()[constraint.component].boxArray();
    for (StencilEntry const& entry : constraint.stencil) {
        if (!faces.contains(facePointBox(entry.cell, constraint.component))) {
            if (reason != nullptr) {
                *reason = "missing staggered DLM constraint support for component "
                        + std::to_string(constraint.component) + " at "
                        + formatIntVect(entry.cell);
            }
            return false;
        }
    }
    return true;
}

} // namespace

void
tagRefinement(amrex::TagBoxArray& tags,
              LevelData const& state,
              std::vector<RigidSphere> const& spheres,
              DLMOptions const& options)
{
    if (!options.enabled || !options.refine) {
        return;
    }

    tagSphereVolumes(tags, state, spheres, options);
    if (options.refine_rigid_mask) {
        tagRigidMask(tags, state.rigidMask());
    }
    if (options.refine_stencil_support && options.scheme == "marker") {
        tagMarkerStencilSupport(tags, state, spheres, options);
    }
}

bool
hasFullStencilSupport(LevelData const& state,
                      std::vector<RigidSphere> const& spheres,
                      DLMOptions const& options,
                      std::string* reason)
{
    if (!options.enabled || options.scheme != "marker") {
        return true;
    }
    if (spheres.empty()) {
        if (reason != nullptr) {
            *reason = "no DLM spheres are available";
        }
        return false;
    }

    std::vector<MarkerPoint> markers = Detail::buildMarkers(state.geom(), spheres, options);
    for (MarkerPoint const& marker : markers) {
        if (!hasMarkerCellSupport(state, marker, reason)) {
            return false;
        }
        for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
            ConstraintPoint constraint = Detail::makeConstraint(state.geom(), options, marker, comp);
            if (!hasConstraintFaceSupport(state, constraint, reason)) {
                return false;
            }
        }
    }
    return true;
}

void
requireFullStencilSupport(LevelData const& state,
                          std::vector<RigidSphere> const& spheres,
                          DLMOptions const& options)
{
    std::string reason;
    if (!hasFullStencilSupport(state, spheres, options, &reason)) {
        amrex::Abort("DLM stencil support is incomplete on the active level: " + reason);
    }
}

} // namespace ELFF::Fluid::DLM
