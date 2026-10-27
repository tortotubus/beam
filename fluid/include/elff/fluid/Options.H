#ifndef ELFF_FLUID_OPTIONS_H_
#define ELFF_FLUID_OPTIONS_H_

#include <AMReX_REAL.H>

#include <array>
#include <string>

namespace ELFF::Fluid {

/** @brief Marker/cell-mask constraints, rigid-body properties, and refinement controls.
 * Read from the dlm parameter namespace. Zero residual tolerances permit
 * iteration-count-driven operation; fail_on_nonconvergence controls failure policy. */
struct DLMOptions
{
    bool enabled = false;
    std::string scheme = "marker";
    std::string mode = "free";
    std::string interpolation = "q1";
    std::string fluid_response = "mass";
    std::string multiplier_layout = "pacific";
    std::string constraint_output;
    std::string history_file;
    amrex::Real radius = 0.15;
    amrex::Real marker_spacing = 1.0;
    amrex::Real interior_margin = 0.5;
    amrex::Real relaxation = 1.0;
    amrex::Real constraint_abs_tol = 0.0;
    amrex::Real constraint_rel_tol = 0.0;
    amrex::Real fluid_density = 1.0;
    amrex::Real solid_density = 2.0;
    int constraint_iterations = 2;
    int refinement_max_level = -1;
    bool explicit_added_mass = false;
    bool use_boundary_points = true;
    bool use_interior_points = true;
    bool enforce_interior_velocity = false;
    bool refine = true;
    bool refine_rigid_mask = true;
    bool refine_stencil_support = true;
    bool require_finest_stencil_support = true;
    bool apply_on_finest_level_only = true;
    bool fail_on_nonconvergence = false;
    bool profile_phases = false;
    amrex::Real refinement_buffer = 0.0;
    std::array<amrex::Real, 3> center{0.5, 0.5, 0.5};
    std::array<amrex::Real, 3> velocity{0.0, 0.0, 0.0};
    std::array<amrex::Real, 3> angular_velocity{0.0, 0.0, 0.0};
    std::array<amrex::Real, 3> gravity{0.0, 0.0, 0.0};
};

/** @brief Physical boundary type and prescribed velocity on one domain face. */
struct BoundaryFaceOptions
{
    std::string type = "natural";
    std::array<amrex::Real, 3> velocity{0.0, 0.0, 0.0};
};

/** @brief Low/high face configuration indexed by spatial direction. */
struct BoundaryOptions
{
    std::array<BoundaryFaceOptions, 3> lo{};
    std::array<BoundaryFaceOptions, 3> hi{};
};

/** @brief Runtime fluid, output, boundary, and coupling configuration.
 * The fluid hierarchy advances with fixed_dt on every level, shortened to
 * stop_time for the final step; there is no time subcycling. */
struct SolverOptions
{
    int max_step = 1;
    amrex::Real stop_time = -1.0;
    amrex::Real fixed_dt = 1.0e-3;
    amrex::Real projection_rel_tol = 1.0e-10;
    amrex::Real projection_abs_tol = 1.0e-12;
    amrex::Real viscosity_rel_tol = 1.0e-10;
    amrex::Real viscosity_abs_tol = 1.0e-12;
    amrex::Real kinematic_viscosity = 0.0;
    std::array<amrex::Real, 3> body_force{0.0, 0.0, 0.0};
    std::array<amrex::Real, 3> pressure_gradient{0.0, 0.0, 0.0};
    bool composite_projection = true;
    bool composite_viscosity = true;
    bool advective_reflux = true;
    bool viscous_reflux = true;
    bool mac_sync = true;
    bool sync_interpolation = true;
    bool sync_projection_residuals = true;
    std::string advection_type = "MOL";
    int plot_int = 0;
    std::string plot_file = "plt";
    bool plot_initial = false;
    int regrid_int = 0;
    BoundaryOptions boundary;
    DLMOptions dlm;
    int verbose = 1;
    bool profile_phases = false;

    /// Read AMReX ParmParse namespaces and apply precision-aware tolerances.
    static SolverOptions read();
    /// Print the effective runtime configuration on the I/O rank.
    void print() const;
};

} // namespace ELFF::Fluid

#endif
