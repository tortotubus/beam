#ifndef ELFF_FLUID_PRECISION_H_
#define ELFF_FLUID_PRECISION_H_

#include <AMReX_REAL.H>

#include <algorithm>
#include <limits>

namespace ELFF::Fluid {

/// @brief Clamp positive solver tolerances to 64 machine epsilons.
/// Nonpositive requests retain their disabling/absolute-tolerance semantics.
// Scale accuracy with the configured amrex::Real, not with the execution
// backend. This gives CPU and GPU builds one policy while avoiding requests
// below the resolution of float32.
inline constexpr amrex::Real
usableRelativeTolerance(amrex::Real requested) noexcept
{
    if (requested <= amrex::Real(0.0)) {
        return requested;
    }
    constexpr amrex::Real floor =
        amrex::Real(64.0) * std::numeric_limits<amrex::Real>::epsilon();
    return std::max(requested, floor);
}

/// @brief Clamp positive test tolerances to 128 machine epsilons.
/// Uses amrex::Real precision, independently of the execution backend.
inline constexpr amrex::Real
usableTestTolerance(amrex::Real requested) noexcept
{
    if (requested <= amrex::Real(0.0)) {
        return requested;
    }
    constexpr amrex::Real floor =
        amrex::Real(128.0) * std::numeric_limits<amrex::Real>::epsilon();
    return std::max(requested, floor);
}

} // namespace ELFF::Fluid

#endif
