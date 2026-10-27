#ifndef ELFF_FLUID_BOUNDARY_H_
#define ELFF_FLUID_BOUNDARY_H_

#include <elff/fluid/DLMTypes.H>
#include <elff/fluid/LevelData.H>
#include <elff/fluid/Options.H>

#include <AMReX_Array.H>
#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_MultiFab.H>
#include <AMReX_SPACE.H>
#include <AMReX_Vector.H>

#include <array>
#include <string>

namespace ELFF::Fluid::Boundary {

/** @brief Supported domain-face boundary conditions. */
enum class BoundaryKind
{
    Natural,
    Periodic,
    Inflow,
    Outflow,
    NoSlip,
    Slip
};

/// Low/high linear-operator boundary conditions indexed by direction.
using LinOpDomainBC =
    std::array<amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM>, 2>;

/** @brief Normalize a boundary name before parsing. */
[[nodiscard]] std::string canonicalBoundaryName(std::string name);
/** @brief Parse a supported boundary name; abort on unknown names. */
[[nodiscard]] BoundaryKind parseBoundaryKind(std::string const& name);
/** @brief Borrow the configured low or high face options in direction dir. */
[[nodiscard]] BoundaryFaceOptions const& faceOptions(BoundaryOptions const& options,
                                                     int dir,
                                                     bool lower);
/** @brief Resolve face configuration against geometry periodicity. */
[[nodiscard]] BoundaryKind effectiveBoundaryKind(amrex::Geometry const& geom,
                                                 SolverOptions const& options,
                                                 int dir,
                                                 bool lower);

/** @brief Build component-wise cell-velocity boundary records. */
[[nodiscard]] amrex::Vector<amrex::BCRec>
makeVelocityBCs(amrex::Geometry const& geom, SolverOptions const& options);

/** @brief Copy boundary records to device-accessible storage for GPU kernels. */
[[nodiscard]] amrex::Gpu::DeviceVector<amrex::BCRec>
makeDeviceBCs(amrex::Vector<amrex::BCRec> const& host_bcs);

/** @brief Map physical conditions to pressure-projection operator conditions. */
[[nodiscard]] LinOpDomainBC projectionBCs(amrex::Geometry const& geom,
                                          SolverOptions const& options);

/** @brief Map physical conditions to viscous velocity-operator conditions. */
[[nodiscard]] LinOpDomainBC velocityOperatorBCs(amrex::Geometry const& geom,
                                                SolverOptions const& options);

/** @brief Fill physical cell-velocity boundaries using the configured face values. */
void fillVelocityPhysicalBoundary(amrex::MultiFab& velocity,
                                  amrex::Geometry const& geom,
                                  SolverOptions const& options);

/** @brief Fill physical velocity boundaries on the supplied level. */
void fillVelocityPhysicalBoundary(LevelData& state, SolverOptions const& options);

/** @brief Enforce physical boundary conditions on staggered face velocities. */
void fillFaceVelocityPhysicalBoundary(DLM::FaceFields& face_velocity,
                                      amrex::Geometry const& geom,
                                      SolverOptions const& options);

} // namespace ELFF::Fluid::Boundary

#endif
