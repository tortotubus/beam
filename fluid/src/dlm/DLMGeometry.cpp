#include "DLMInternal.H"

#include <AMReX_ParallelDescriptor.H>

#include <fstream>
#include <iomanip>
#include <set>
#include <tuple>

namespace ELFF::Fluid::DLM::Detail {

amrex::Real
cellVolume(amrex::Geometry const& geom)
{
    amrex::Real volume = 1.0;
    const auto* dx = geom.CellSize();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        volume *= dx[dir];
    }
    return volume;
}

amrex::Real
minCellSize(amrex::Geometry const& geom)
{
    const auto* dx = geom.CellSize();
    amrex::Real min_dx = dx[0];
    for (int dir = 1; dir < AMREX_SPACEDIM; ++dir) {
        min_dx = std::min(min_dx, dx[dir]);
    }
    return min_dx;
}

amrex::IntVect
makeCell(int i, int j, int k)
{
#if (AMREX_SPACEDIM == 2)
    static_cast<void>(k);
    return amrex::IntVect(i, j);
#else
    return amrex::IntVect(i, j, k);
#endif
}

int
cellK(amrex::IntVect const& iv)
{
#if (AMREX_SPACEDIM == 2)
    static_cast<void>(iv);
    return 0;
#else
    return iv[2];
#endif
}

std::array<amrex::Real, 3>
cellCenter(amrex::Geometry const& geom, amrex::IntVect const& iv)
{
    const auto* dx = geom.CellSize();
    const auto* plo = geom.ProbLo();
    std::array<amrex::Real, 3> point{0.0, 0.0, 0.0};
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        point[dir] = plo[dir] + (amrex::Real(iv[dir]) + amrex::Real(0.5)) * dx[dir];
    }
    return point;
}

std::array<amrex::Real, 3>
faceCenter(amrex::Geometry const& geom, amrex::IntVect const& iv, int component)
{
    const auto* dx = geom.CellSize();
    const auto* plo = geom.ProbLo();
    std::array<amrex::Real, 3> point{0.0, 0.0, 0.0};
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const amrex::Real offset = (dir == component) ? amrex::Real(0.0)
                                                      : amrex::Real(0.5);
        point[dir] = plo[dir] + (amrex::Real(iv[dir]) + offset) * dx[dir];
    }
    return point;
}

RigidSphere
optionSphere(DLMOptions const& options)
{
    RigidSphere sphere;
    sphere.id = 1;
    sphere.geometry_id = 0;
    sphere.prescribed = (options.mode == "prescribed");
    sphere.radius = options.radius;
    sphere.density = options.solid_density;
    sphere.center = options.center;
    sphere.velocity = options.velocity;
    sphere.angular_velocity = options.angular_velocity;
    return sphere;
}

amrex::Real
distanceSquared(std::array<amrex::Real, 3> const& point, DLMOptions const& options)
{
    amrex::Real dist_sq = 0.0;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const amrex::Real d = point[dir] - options.center[dir];
        dist_sq += d * d;
    }
    return dist_sq;
}

amrex::Real
distanceSquared(std::array<amrex::Real, 3> const& point, BodyState const& body)
{
    amrex::Real dist_sq = 0.0;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const amrex::Real d = point[dir] - body.center[dir];
        dist_sq += d * d;
    }
    return dist_sq;
}

bool
sameCell(amrex::IntVect const& a, amrex::IntVect const& b)
{
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (a[dir] != b[dir]) {
            return false;
        }
    }
    return true;
}

int
domainIndex(amrex::Geometry const& geom, int dir, int idx)
{
    const amrex::Box& domain = geom.Domain();
    const int lo = domain.smallEnd(dir);
    const int hi = domain.bigEnd(dir);
    if (geom.isPeriodic(dir)) {
        const int n = hi - lo + 1;
        while (idx < lo) {
            idx += n;
        }
        while (idx > hi) {
            idx -= n;
        }
        return idx;
    }
    return std::clamp(idx, lo, hi);
}

void
addStencilEntry(std::vector<StencilEntry>& stencil, amrex::IntVect const& cell,
                amrex::Real weight)
{
    if (weight == 0.0) {
        return;
    }
    for (StencilEntry& entry : stencil) {
        if (sameCell(entry.cell, cell)) {
            entry.weight += weight;
            return;
        }
    }
    stencil.push_back({cell, weight});
}

std::vector<StencilEntry>
q1Stencil(amrex::Geometry const& geom, std::array<amrex::Real, 3> const& point)
{
    std::array<std::array<int, 2>, 3> cells{};
    std::array<std::array<amrex::Real, 2>, 3> weights{};
    const auto* dx = geom.CellSize();
    const auto* plo = geom.ProbLo();
    const amrex::Box& domain = geom.Domain();

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const amrex::Real grid_pos = (point[dir] - plo[dir]) / dx[dir] - amrex::Real(0.5);
        int lower = static_cast<int>(std::floor(grid_pos));
        amrex::Real frac = grid_pos - amrex::Real(lower);
        if (!geom.isPeriodic(dir)) {
            if (lower < domain.smallEnd(dir)) {
                lower = domain.smallEnd(dir);
                frac = 0.0;
            } else if (lower >= domain.bigEnd(dir)) {
                lower = domain.bigEnd(dir);
                frac = 0.0;
            }
        }
        frac = std::clamp(frac, amrex::Real(0.0), amrex::Real(1.0));
        cells[dir][0] = domainIndex(geom, dir, lower);
        cells[dir][1] = domainIndex(geom, dir, lower + 1);
        weights[dir][0] = 1.0 - frac;
        weights[dir][1] = frac;
    }

#if (AMREX_SPACEDIM == 2)
    cells[2][0] = 0;
    cells[2][1] = 0;
    weights[2][0] = 1.0;
    weights[2][1] = 0.0;
#endif

    std::vector<StencilEntry> stencil;
    stencil.reserve(1 << AMREX_SPACEDIM);
    for (int ii = 0; ii < 2; ++ii) {
        for (int jj = 0; jj < 2; ++jj) {
            for (int kk = 0; kk < 2; ++kk) {
                const amrex::Real weight = weights[0][ii] * weights[1][jj] * weights[2][kk];
                addStencilEntry(stencil,
                                makeCell(cells[0][ii], cells[1][jj], cells[2][kk]),
                                weight);
            }
        }
    }

    amrex::Real sum = 0.0;
    for (StencilEntry const& entry : stencil) {
        sum += entry.weight;
    }
    if (sum > 0.0) {
        for (StencilEntry& entry : stencil) {
            entry.weight /= sum;
        }
    }
    return stencil;
}

int
domainStaggeredIndex(amrex::Geometry const& geom, int dir, int component, int idx)
{
    const amrex::Box& domain = geom.Domain();
    const int lo = domain.smallEnd(dir);
    const int hi = domain.bigEnd(dir) + ((dir == component) ? 1 : 0);
    if (geom.isPeriodic(dir)) {
        const int n = domain.length(dir);
        while (idx < lo) {
            idx += n;
        }
        while (idx > hi) {
            idx -= n;
        }
        return idx;
    }
    return std::clamp(idx, lo, hi);
}

std::vector<StencilEntry>
q1StaggeredStencil(amrex::Geometry const& geom,
                   std::array<amrex::Real, 3> const& point,
                   int component)
{
    std::array<std::array<int, 2>, 3> faces{};
    std::array<std::array<amrex::Real, 2>, 3> weights{};
    const auto* dx = geom.CellSize();
    const auto* plo = geom.ProbLo();
    const amrex::Box& domain = geom.Domain();

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const amrex::Real offset = (dir == component) ? amrex::Real(0.0)
                                                      : amrex::Real(0.5);
        const amrex::Real grid_pos = (point[dir] - plo[dir]) / dx[dir] - offset;
        int lower = static_cast<int>(std::floor(grid_pos));
        amrex::Real frac = grid_pos - amrex::Real(lower);
        const int lo = domain.smallEnd(dir);
        const int hi = domain.bigEnd(dir) + ((dir == component) ? 1 : 0);
        if (!geom.isPeriodic(dir)) {
            if (lower < lo) {
                lower = lo;
                frac = 0.0;
            } else if (lower >= hi) {
                lower = hi;
                frac = 0.0;
            }
        }
        frac = std::clamp(frac, amrex::Real(0.0), amrex::Real(1.0));
        faces[dir][0] = domainStaggeredIndex(geom, dir, component, lower);
        faces[dir][1] = domainStaggeredIndex(geom, dir, component, lower + 1);
        weights[dir][0] = 1.0 - frac;
        weights[dir][1] = frac;
    }

#if (AMREX_SPACEDIM == 2)
    faces[2][0] = 0;
    faces[2][1] = 0;
    weights[2][0] = 1.0;
    weights[2][1] = 0.0;
#endif

    std::vector<StencilEntry> stencil;
    stencil.reserve(1 << AMREX_SPACEDIM);
    for (int ii = 0; ii < 2; ++ii) {
        for (int jj = 0; jj < 2; ++jj) {
            for (int kk = 0; kk < 2; ++kk) {
                const amrex::Real weight = weights[0][ii] * weights[1][jj] * weights[2][kk];
                addStencilEntry(stencil,
                                makeCell(faces[0][ii], faces[1][jj], faces[2][kk]),
                                weight);
            }
        }
    }

    amrex::Real sum = 0.0;
    for (StencilEntry const& entry : stencil) {
        sum += entry.weight;
    }
    if (sum > 0.0) {
        for (StencilEntry& entry : stencil) {
            entry.weight /= sum;
        }
    }
    return stencil;
}

amrex::Real
q2Basis(int node, amrex::Real t)
{
    switch (node) {
    case 0:
        return amrex::Real(2.0) * (t - amrex::Real(0.5)) * (t - amrex::Real(1.0));
    case 1:
        return amrex::Real(4.0) * t * (amrex::Real(1.0) - t);
    case 2:
        return amrex::Real(2.0) * t * (t - amrex::Real(0.5));
    default:
        amrex::Abort("Invalid Q2 basis node");
    }
    return 0.0;
}

bool
supportedDlmInterpolation(std::string const& interpolation)
{
    return interpolation == "q1" || interpolation == "q2";
}

bool
supportedMultiplierLayout(std::string const& layout)
{
    return layout == "legacy" || layout == "pacific";
}

std::vector<StencilEntry>
q2StaggeredStencil(amrex::Geometry const& geom,
                   std::array<amrex::Real, 3> const& point,
                   std::array<amrex::Real, 3> const& center_vector,
                   int component)
{
    std::array<std::array<int, 3>, 3> faces{};
    std::array<std::array<amrex::Real, 3>, 3> weights{};
    const auto* dx = geom.CellSize();
    const auto* plo = geom.ProbLo();
    const amrex::Box& domain = geom.Domain();

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const amrex::Real offset = (dir == component) ? amrex::Real(0.0)
                                                      : amrex::Real(0.5);
        const amrex::Real grid_pos = (point[dir] - plo[dir]) / dx[dir] - offset;
        const int lower = static_cast<int>(std::floor(grid_pos));
        int base = (center_vector[dir] >= amrex::Real(0.0)) ? lower : lower - 1;

        const int lo = domain.smallEnd(dir);
        const int hi = domain.bigEnd(dir) + ((dir == component) ? 1 : 0);
        if (!geom.isPeriodic(dir)) {
            if (hi - lo + 1 < 3) {
                return q1StaggeredStencil(geom, point, component);
            }
            base = std::clamp(base, lo, hi - 2);
        }

        const amrex::Real t = (grid_pos - amrex::Real(base)) / amrex::Real(2.0);
        for (int node = 0; node < 3; ++node) {
            faces[dir][node] = domainStaggeredIndex(geom, dir, component, base + node);
            weights[dir][node] = q2Basis(node, t);
        }
    }

#if (AMREX_SPACEDIM == 2)
    faces[2][0] = 0;
    faces[2][1] = 0;
    faces[2][2] = 0;
    weights[2][0] = 1.0;
    weights[2][1] = 0.0;
    weights[2][2] = 0.0;
#endif

    std::vector<StencilEntry> stencil;
    stencil.reserve(AMREX_SPACEDIM == 2 ? 9 : 27);
    for (int ii = 0; ii < 3; ++ii) {
        for (int jj = 0; jj < 3; ++jj) {
            for (int kk = 0; kk < 3; ++kk) {
                const amrex::Real weight = weights[0][ii] * weights[1][jj] * weights[2][kk];
                addStencilEntry(stencil,
                                makeCell(faces[0][ii], faces[1][jj], faces[2][kk]),
                                weight);
            }
        }
    }

    amrex::Real sum = 0.0;
    for (StencilEntry const& entry : stencil) {
        sum += entry.weight;
    }
    if (sum != 0.0) {
        for (StencilEntry& entry : stencil) {
            entry.weight /= sum;
        }
    }
    return stencil;
}

MarkerPoint
makeMarker(amrex::Geometry const& geom,
           DLMOptions const& options,
           RigidSphere const& sphere,
           std::size_t body_index,
           std::array<amrex::Real, 3> position,
           amrex::Real quadrature_weight,
           bool boundary)
{
    MarkerPoint marker;
    marker.body_index = body_index;
    marker.position = position;
    marker.quadrature_weight = quadrature_weight;
    marker.boundary = boundary;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        marker.center_vector[dir] = marker.position[dir] - sphere.center[dir];
    }
    marker.stencil = q1Stencil(geom, marker.position);
    static_cast<void>(options);
    return marker;
}

void
appendBoundaryMarkers(amrex::Geometry const& geom,
                      DLMOptions const& options,
                      RigidSphere const& sphere,
                      std::size_t body_index,
                      std::vector<MarkerPoint>& markers)
{
    if (!options.use_boundary_points) {
        return;
    }

    const bool pacific_layout = options.multiplier_layout == "pacific";
    const amrex::Real critical_distance =
        std::sqrt(amrex::Real(AMREX_SPACEDIM)) * minCellSize(geom);
    const amrex::Real spacing = std::max(
        options.marker_spacing * (pacific_layout ? critical_distance : minCellSize(geom)),
        amrex::Real(1.0e-12));
#if (AMREX_SPACEDIM == 2)
    const amrex::Real circumference = 2.0 * pi * sphere.radius;
    const int nmarkers = pacific_layout
        ? std::max(8, static_cast<int>(circumference / spacing))
        : std::max(8, static_cast<int>(std::ceil(circumference / spacing)));
    const amrex::Real qweight = pacific_layout
        ? amrex::Real(1.0)
        : circumference * minCellSize(geom) / amrex::Real(nmarkers);
    for (int n = 0; n < nmarkers; ++n) {
        const amrex::Real theta = 2.0 * pi * (amrex::Real(n) + amrex::Real(0.5))
                                / amrex::Real(nmarkers);
        const std::array<amrex::Real, 3> point{
            sphere.center[0] + sphere.radius * std::cos(theta),
            sphere.center[1] + sphere.radius * std::sin(theta),
            sphere.center[2]
        };
        markers.push_back(makeMarker(geom, options, sphere, body_index, point, qweight, true));
    }
#else
    const amrex::Real surface_area = 4.0 * pi * sphere.radius * sphere.radius;
    if (pacific_layout) {
        const int nmarkers = std::max(
            12,
            static_cast<int>(std::pow(amrex::Real(3.809) * sphere.radius / spacing,
                                      amrex::Real(2.0))));
        const amrex::Real spiral_correction = spacing / std::sqrt(amrex::Real(3.0));
        const amrex::Real cspiral = amrex::Real(3.6) / std::sqrt(amrex::Real(nmarkers));
        amrex::Real phi_previous = 0.0;
        amrex::Real dphi = 0.0;
        for (int n = 0; n < nmarkers; ++n) {
            const amrex::Real h = -1.0 + 2.0 * amrex::Real(n)
                                           / amrex::Real(nmarkers - 1);
            amrex::Real theta = std::acos(h);
            amrex::Real phi = 0.0;
            if (n == 0) {
                theta = pi - amrex::Real(0.5) * spiral_correction / sphere.radius;
            } else if (n == nmarkers - 1) {
                phi = phi_previous + dphi;
                if (phi > amrex::Real(2.0) * pi) {
                    phi -= amrex::Real(2.0) * pi;
                }
                theta = amrex::Real(0.5) * spiral_correction / sphere.radius;
            } else {
                dphi = cspiral / std::sqrt(amrex::Real(1.0) - h * h);
                phi = phi_previous + dphi;
                if (phi > amrex::Real(2.0) * pi) {
                    phi -= amrex::Real(2.0) * pi;
                }
            }
            phi_previous = phi;
            if (n == 1) {
                theta -= amrex::Real(0.4) * spiral_correction / sphere.radius;
            }
            if (n == nmarkers - 2) {
                phi -= amrex::Real(0.1) * dphi;
                theta += amrex::Real(0.25) * spiral_correction / sphere.radius;
            }
            const std::array<amrex::Real, 3> point{
                sphere.center[0] + sphere.radius * std::cos(phi) * std::sin(theta),
                sphere.center[1] + sphere.radius * std::sin(phi) * std::sin(theta),
                sphere.center[2] + sphere.radius * std::cos(theta)
            };
            markers.push_back(
                makeMarker(geom, options, sphere, body_index, point, amrex::Real(1.0), true));
        }
    } else {
        const int nmarkers = std::max(
            12, static_cast<int>(std::ceil(surface_area / (spacing * spacing))));
        const amrex::Real qweight =
            surface_area * minCellSize(geom) / amrex::Real(nmarkers);
        const amrex::Real golden =
            pi * (amrex::Real(3.0) - std::sqrt(amrex::Real(5.0)));
        for (int n = 0; n < nmarkers; ++n) {
            const amrex::Real z =
                amrex::Real(1.0)
                - amrex::Real(2.0) * (amrex::Real(n) + amrex::Real(0.5))
                      / amrex::Real(nmarkers);
            const amrex::Real radius_xy =
                std::sqrt(std::max(amrex::Real(0.0), amrex::Real(1.0) - z * z));
            const amrex::Real theta = golden * amrex::Real(n);
            const std::array<amrex::Real, 3> point{
                sphere.center[0] + sphere.radius * radius_xy * std::cos(theta),
                sphere.center[1] + sphere.radius * radius_xy * std::sin(theta),
                sphere.center[2] + sphere.radius * z
            };
            markers.push_back(
                makeMarker(geom, options, sphere, body_index, point, qweight, true));
        }
    }
#endif
}

void
appendInteriorMarkers(amrex::Geometry const& geom,
                      DLMOptions const& options,
                      RigidSphere const& sphere,
                      std::size_t body_index,
                      std::vector<MarkerPoint>& markers)
{
    if (!options.use_interior_points) {
        return;
    }
    if (options.multiplier_layout == "pacific") {
        return;
    }

    const amrex::Real interior_radius =
        std::max(amrex::Real(0.0),
                 sphere.radius - options.interior_margin * minCellSize(geom));
    if (interior_radius <= 0.0) {
        return;
    }

    const amrex::Real radius_sq = interior_radius * interior_radius;
    const amrex::Real qweight = cellVolume(geom);
    for (amrex::BoxIterator bit(geom.Domain()); bit.ok(); ++bit) {
        const amrex::IntVect iv = bit();
        const auto point = cellCenter(geom, iv);
        amrex::Real dist_sq = 0.0;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            const amrex::Real d = point[dir] - sphere.center[dir];
            dist_sq += d * d;
        }
        if (dist_sq <= radius_sq) {
            markers.push_back(makeMarker(geom, options, sphere, body_index, point, qweight, false));
        }
    }
}

void
appendMarkers(amrex::Geometry const& geom,
              DLMOptions const& options,
              RigidSphere const& sphere,
              std::size_t body_index,
              std::vector<MarkerPoint>& markers)
{
    if (sphere.radius <= 0.0) {
        amrex::Abort("DLM sphere radius must be positive");
    }
    appendBoundaryMarkers(geom, options, sphere, body_index, markers);
    appendInteriorMarkers(geom, options, sphere, body_index, markers);
}

std::vector<MarkerPoint>
buildMarkers(amrex::Geometry const& geom, DLMOptions const& options)
{
    if (options.radius <= 0.0) {
        amrex::Abort("dlm.radius must be positive when dlm.enabled = true");
    }
    if (options.marker_spacing <= 0.0) {
        amrex::Abort("dlm.marker_spacing must be positive");
    }
    if (!supportedDlmInterpolation(options.interpolation)) {
        amrex::Abort("Supported marker DLM interpolation options are q1 and q2");
    }
    if (!supportedMultiplierLayout(options.multiplier_layout)) {
        amrex::Abort("dlm.multiplier_layout must be either 'legacy' or 'pacific'");
    }

    std::vector<MarkerPoint> markers;
    const RigidSphere sphere = optionSphere(options);
    appendMarkers(geom, options, sphere, 0, markers);
    if (markers.empty()) {
        amrex::Abort("DLM marker scheme generated no multiplier points");
    }
    return markers;
}

std::vector<MarkerPoint>
buildMarkers(amrex::Geometry const& geom,
             std::vector<RigidSphere> const& spheres,
             DLMOptions const& options)
{
    if (options.marker_spacing <= 0.0) {
        amrex::Abort("dlm.marker_spacing must be positive");
    }
    if (!supportedDlmInterpolation(options.interpolation)) {
        amrex::Abort("Supported marker DLM interpolation options are q1 and q2");
    }
    if (!supportedMultiplierLayout(options.multiplier_layout)) {
        amrex::Abort("dlm.multiplier_layout must be either 'legacy' or 'pacific'");
    }
    if (spheres.empty()) {
        amrex::Abort("DLM marker scheme requires at least one sphere");
    }

    std::vector<MarkerPoint> markers;
    for (std::size_t body_index = 0; body_index < spheres.size(); ++body_index) {
        appendMarkers(geom, options, spheres[body_index], body_index, markers);
    }
    if (markers.empty()) {
        amrex::Abort("DLM marker scheme generated no multiplier points");
    }
    return markers;
}

ConstraintPoint
makeConstraint(amrex::Geometry const& geom,
               DLMOptions const& options,
               MarkerPoint const& marker,
               int component)
{
    ConstraintPoint constraint;
    constraint.body_index = marker.body_index;
    constraint.component = component;
    constraint.position = marker.position;
    constraint.center_vector = marker.center_vector;
    constraint.quadrature_weight = marker.quadrature_weight;
    constraint.boundary = marker.boundary;
    if (options.interpolation == "q2" && marker.boundary) {
        constraint.stencil =
            q2StaggeredStencil(geom, marker.position, marker.center_vector, component);
    } else {
        constraint.stencil = q1StaggeredStencil(geom, marker.position, component);
    }
    if (constraint.stencil.size() == 1) {
        constraint.dof = constraint.stencil.front().cell;
    }
    return constraint;
}

ConstraintPoint
makeInteriorFaceConstraint(amrex::Geometry const& geom,
                           RigidSphere const& sphere,
                           std::size_t body_index,
                           amrex::IntVect const& face,
                           int component,
                           amrex::Real quadrature_weight)
{
    ConstraintPoint constraint;
    constraint.body_index = body_index;
    constraint.component = component;
    constraint.dof = face;
    constraint.position = faceCenter(geom, face, component);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        constraint.center_vector[dir] = constraint.position[dir] - sphere.center[dir];
    }
    constraint.quadrature_weight = quadrature_weight;
    constraint.boundary = false;
    constraint.stencil.push_back({face, amrex::Real(1.0)});
    return constraint;
}

void
appendBoundaryMarkerConstraints(amrex::Geometry const& geom,
                                DLMOptions const& options,
                                std::vector<MarkerPoint> const& markers,
                                std::vector<ConstraintPoint>& constraints)
{
    for (MarkerPoint const& marker : markers) {
        if (!marker.boundary) {
            continue;
        }
        for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
            constraints.push_back(makeConstraint(geom, options, marker, comp));
        }
    }
}

void
appendInteriorFaceConstraints(amrex::Geometry const& geom,
                              DLMOptions const& options,
                              RigidSphere const& sphere,
                              std::size_t body_index,
                              std::vector<MarkerPoint> const& markers,
                              std::set<std::tuple<int, int, int, int>>& occupied_faces,
                              std::vector<ConstraintPoint>& constraints)
{
    if (!options.use_interior_points) {
        return;
    }

    const bool pacific_layout = options.multiplier_layout == "pacific";
    const amrex::Real dx = minCellSize(geom);
    const amrex::Real interior_radius = pacific_layout
        ? sphere.radius
        : std::max(amrex::Real(0.0), sphere.radius - options.interior_margin * dx);
    if (interior_radius <= 0.0) {
        return;
    }

    const amrex::Real radius_sq = interior_radius * interior_radius;
    const amrex::Box& domain = geom.Domain();
    for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
        const amrex::Box face_domain = amrex::surroundingNodes(domain, comp);
        for (amrex::BoxIterator bit(face_domain); bit.ok(); ++bit) {
            const amrex::IntVect iv = bit();
            const auto point = faceCenter(geom, iv, comp);
            amrex::Real dist_sq = 0.0;
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                const amrex::Real d = point[dir] - sphere.center[dir];
                dist_sq += d * d;
            }
            if (dist_sq > radius_sq) {
                continue;
            }

            if (pacific_layout && dist_sq > (sphere.radius - dx) * (sphere.radius - dx)) {
                bool too_close_to_boundary = false;
                for (MarkerPoint const& marker : markers) {
                    if (!marker.boundary || marker.body_index != body_index) {
                        continue;
                    }
                    amrex::Real marker_distance_sq = 0.0;
                    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                        const amrex::Real d = point[dir] - marker.position[dir];
                        marker_distance_sq += d * d;
                    }
                    if (marker_distance_sq < dx * dx) {
                        too_close_to_boundary = true;
                        break;
                    }
                }
                if (too_close_to_boundary) {
                    continue;
                }
            }

            const auto key = std::make_tuple(comp, iv[0], iv[1], cellK(iv));
            if (!occupied_faces.insert(key).second) {
                continue;
            }
            constraints.push_back(makeInteriorFaceConstraint(
                geom,
                sphere,
                body_index,
                iv,
                comp,
                pacific_layout ? amrex::Real(1.0) : cellVolume(geom)));
        }
    }
}

std::vector<ConstraintPoint>
buildConstraints(amrex::Geometry const& geom, DLMOptions const& options)
{
    std::vector<MarkerPoint> markers = buildMarkers(geom, options);
    return buildConstraints(geom,
                            markers,
                            std::vector<RigidSphere>{optionSphere(options)},
                            options);
}

std::vector<ConstraintPoint>
buildConstraints(amrex::Geometry const& geom,
                 std::vector<MarkerPoint> const& markers,
                 std::vector<RigidSphere> const& spheres,
                 DLMOptions const& options)
{
    std::vector<ConstraintPoint> constraints;
    constraints.reserve(markers.size() * AMREX_SPACEDIM);
    appendBoundaryMarkerConstraints(geom, options, markers, constraints);

    std::set<std::tuple<int, int, int, int>> occupied_faces;
    if (options.multiplier_layout != "pacific") {
        for (ConstraintPoint const& constraint : constraints) {
            for (StencilEntry const& entry : constraint.stencil) {
                occupied_faces.insert(std::make_tuple(constraint.component,
                                                      entry.cell[0],
                                                      entry.cell[1],
                                                      cellK(entry.cell)));
            }
        }
    }
    for (std::size_t body_index = 0; body_index < spheres.size(); ++body_index) {
        appendInteriorFaceConstraints(geom,
                                      options,
                                      spheres[body_index],
                                      body_index,
                                      markers,
                                      occupied_faces,
                                      constraints);
    }
    if (constraints.empty()) {
        amrex::Abort("DLM marker scheme generated no constraints");
    }
    return constraints;
}

std::vector<ConstraintPoint>
buildConstraints(amrex::Geometry const& geom,
                 std::vector<RigidSphere> const& spheres,
                 DLMOptions const& options)
{
    const std::vector<MarkerPoint> markers = buildMarkers(geom, spheres, options);
    return buildConstraints(geom, markers, spheres, options);
}

void
writeConstraintCatalog(std::vector<ConstraintPoint> const& constraints,
                       std::string const& filename)
{
    if (filename.empty() || !amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

    std::ofstream output(filename, std::ios::out | std::ios::trunc);
    if (!output) {
        amrex::Abort("Unable to open DLM constraint catalog: " + filename);
    }
    output << "constraint,body,component,boundary,x,y,z,quadrature_weight,"
              "stencil_entry,i,j,k,stencil_weight\n";
    output << std::setprecision(17);
    for (std::size_t id = 0; id < constraints.size(); ++id) {
        ConstraintPoint const& constraint = constraints[id];
        for (std::size_t entry_id = 0; entry_id < constraint.stencil.size(); ++entry_id) {
            StencilEntry const& entry = constraint.stencil[entry_id];
            output << id << ','
                   << constraint.body_index << ','
                   << constraint.component << ','
                   << static_cast<int>(constraint.boundary) << ','
                   << constraint.position[0] << ','
                   << constraint.position[1] << ','
                   << constraint.position[2] << ','
                   << constraint.quadrature_weight << ','
                   << entry_id << ','
                   << entry.cell[0] << ','
                   << entry.cell[1] << ','
                   << cellK(entry.cell) << ','
                   << entry.weight << '\n';
        }
    }
}

amrex::Real
bodyVolume(DLMOptions const& options)
{
    return sphereVolume(options.radius);
}

amrex::Real
sphereVolume(amrex::Real radius)
{
#if (AMREX_SPACEDIM == 2)
    return pi * radius * radius;
#else
    return (amrex::Real(4.0) / amrex::Real(3.0)) * pi
         * radius * radius * radius;
#endif
}

} // namespace ELFF::Fluid::DLM::Detail
