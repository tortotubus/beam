#include <elff/fluid/Options.H>
#include <elff/fluid/Precision.H>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_SPACE.H>

#include <algorithm>
#include <array>

namespace ELFF::Fluid {

namespace {

bool
queryArray(amrex::ParmParse& pp,
           std::string const& name,
           std::array<amrex::Real, 3>& values,
           int count)
{
    amrex::Vector<amrex::Real> parsed(count);
    for (int i = 0; i < count; ++i) {
        parsed[i] = values[i];
    }

    if (pp.queryarr(name.c_str(), parsed, 0, count)) {
        for (int i = 0; i < count; ++i) {
            values[i] = parsed[i];
        }
        return true;
    }
    return false;
}

void
readBoundaryFace(amrex::ParmParse& pp,
                 std::string const& name,
                 BoundaryFaceOptions& face,
                 std::array<amrex::Real, 3> const& default_inflow_velocity,
                 bool has_default_inflow_velocity)
{
    pp.query(name.c_str(), face.type);
    const bool has_face_velocity =
        queryArray(pp, name + "_velocity", face.velocity, AMREX_SPACEDIM);
    if (!has_face_velocity && has_default_inflow_velocity && face.type == "inflow") {
        face.velocity = default_inflow_velocity;
    }
}

} // namespace

SolverOptions
SolverOptions::read()
{
    SolverOptions options;

    amrex::ParmParse pp;
    pp.query("max_step", options.max_step);
    pp.query("stop_time", options.stop_time);

    amrex::ParmParse pp_fluid("fluid");
    pp_fluid.query("max_step", options.max_step);
    pp_fluid.query("stop_time", options.stop_time);
    pp_fluid.query("fixed_dt", options.fixed_dt);
    pp_fluid.query("projection_rel_tol", options.projection_rel_tol);
    pp_fluid.query("projection_abs_tol", options.projection_abs_tol);
    pp_fluid.query("viscosity_rel_tol", options.viscosity_rel_tol);
    pp_fluid.query("viscosity_abs_tol", options.viscosity_abs_tol);
    options.projection_rel_tol = usableRelativeTolerance(options.projection_rel_tol);
    options.viscosity_rel_tol = usableRelativeTolerance(options.viscosity_rel_tol);
    pp_fluid.query("kinematic_viscosity", options.kinematic_viscosity);
    queryArray(pp_fluid, "body_force", options.body_force, AMREX_SPACEDIM);
    queryArray(pp_fluid, "pressure_gradient", options.pressure_gradient, AMREX_SPACEDIM);
    pp_fluid.query("composite_projection", options.composite_projection);
    pp_fluid.query("composite_viscosity", options.composite_viscosity);
    pp_fluid.query("advective_reflux", options.advective_reflux);
    pp_fluid.query("viscous_reflux", options.viscous_reflux);
    pp_fluid.query("mac_sync", options.mac_sync);
    pp_fluid.query("sync_interpolation", options.sync_interpolation);
    pp_fluid.query("sync_projection_residuals", options.sync_projection_residuals);
    pp_fluid.query("advection_type", options.advection_type);
    pp_fluid.query("verbose", options.verbose);
    pp_fluid.query("profile_phases", options.profile_phases);

    amrex::ParmParse pp_amr("amr");
    pp_amr.query("plot_int", options.plot_int);
    pp_amr.query("plot_file", options.plot_file);
    pp_amr.query("plot_initial", options.plot_initial);
    pp_amr.query("regrid_int", options.regrid_int);

    pp_fluid.query("plot_int", options.plot_int);
    pp_fluid.query("plot_file", options.plot_file);
    pp_fluid.query("plot_initial", options.plot_initial);
    pp_fluid.query("regrid_int", options.regrid_int);

    amrex::ParmParse pp_bc("bc");
    std::array<amrex::Real, 3> default_inflow_velocity{0.0, 0.0, 0.0};
    const bool has_default_inflow_velocity =
        queryArray(pp_bc, "inflow_velocity", default_inflow_velocity, AMREX_SPACEDIM);
    constexpr std::array<char const*, 3> dir_names{"x", "y", "z"};
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        readBoundaryFace(pp_bc,
                         std::string(dir_names[dir]) + "lo",
                         options.boundary.lo[dir],
                         default_inflow_velocity,
                         has_default_inflow_velocity);
        readBoundaryFace(pp_bc,
                         std::string(dir_names[dir]) + "hi",
                         options.boundary.hi[dir],
                         default_inflow_velocity,
                         has_default_inflow_velocity);
    }

    amrex::ParmParse pp_dlm("dlm");
    pp_dlm.query("enabled", options.dlm.enabled);
    pp_dlm.query("scheme", options.dlm.scheme);
    pp_dlm.query("mode", options.dlm.mode);
    pp_dlm.query("interpolation", options.dlm.interpolation);
    pp_dlm.query("fluid_response", options.dlm.fluid_response);
    pp_dlm.query("multiplier_layout", options.dlm.multiplier_layout);
    pp_dlm.query("constraint_output", options.dlm.constraint_output);
    pp_dlm.query("history_file", options.dlm.history_file);
    pp_dlm.query("radius", options.dlm.radius);
    pp_dlm.query("marker_spacing", options.dlm.marker_spacing);
    pp_dlm.query("interior_margin", options.dlm.interior_margin);
    pp_dlm.query("relaxation", options.dlm.relaxation);
    pp_dlm.query("constraint_abs_tol", options.dlm.constraint_abs_tol);
    pp_dlm.query("constraint_rel_tol", options.dlm.constraint_rel_tol);
    pp_dlm.query("fluid_density", options.dlm.fluid_density);
    pp_dlm.query("solid_density", options.dlm.solid_density);
    pp_dlm.query("constraint_iterations", options.dlm.constraint_iterations);
    pp_dlm.query("refinement_max_level", options.dlm.refinement_max_level);
    pp_dlm.query("explicit_added_mass", options.dlm.explicit_added_mass);
    pp_dlm.query("use_boundary_points", options.dlm.use_boundary_points);
    pp_dlm.query("use_interior_points", options.dlm.use_interior_points);
    pp_dlm.query("enforce_interior_velocity", options.dlm.enforce_interior_velocity);
    pp_dlm.query("refine", options.dlm.refine);
    pp_dlm.query("refine_rigid_mask", options.dlm.refine_rigid_mask);
    pp_dlm.query("refine_stencil_support", options.dlm.refine_stencil_support);
    pp_dlm.query("require_finest_stencil_support",
                 options.dlm.require_finest_stencil_support);
    pp_dlm.query("apply_on_finest_level_only",
                 options.dlm.apply_on_finest_level_only);
    pp_dlm.query("fail_on_nonconvergence", options.dlm.fail_on_nonconvergence);
    pp_dlm.query("refinement_buffer", options.dlm.refinement_buffer);
    queryArray(pp_dlm, "center", options.dlm.center, AMREX_SPACEDIM);
    queryArray(pp_dlm, "velocity", options.dlm.velocity, AMREX_SPACEDIM);
    queryArray(pp_dlm, "angular_velocity", options.dlm.angular_velocity, 3);
    queryArray(pp_dlm, "gravity", options.dlm.gravity, AMREX_SPACEDIM);

    return options;
}

void
SolverOptions::print() const
{
    if (verbose <= 0) {
        return;
    }

    amrex::Print() << "ELFF fluid solver options:\n"
                   << "  max_step = " << max_step << "\n"
                   << "  stop_time = " << stop_time << "\n"
                   << "  fixed_dt = " << fixed_dt << "\n"
                   << "  projection_rel_tol = " << projection_rel_tol << "\n"
                   << "  projection_abs_tol = " << projection_abs_tol << "\n"
                   << "  viscosity_rel_tol = " << viscosity_rel_tol << "\n"
                   << "  viscosity_abs_tol = " << viscosity_abs_tol << "\n"
                   << "  kinematic_viscosity = " << kinematic_viscosity << "\n"
                   << "  body_force = "
                   << body_force[0] << " " << body_force[1]
#if (AMREX_SPACEDIM == 3)
                   << " " << body_force[2]
#endif
                   << "\n"
                   << "  pressure_gradient = "
                   << pressure_gradient[0] << " " << pressure_gradient[1]
#if (AMREX_SPACEDIM == 3)
                   << " " << pressure_gradient[2]
#endif
                   << "\n"
                   << "  composite_projection = " << composite_projection << "\n"
                   << "  composite_viscosity = " << composite_viscosity << "\n"
                   << "  advective_reflux = " << advective_reflux << "\n"
                   << "  viscous_reflux = " << viscous_reflux << "\n"
                   << "  mac_sync = " << mac_sync << "\n"
                   << "  sync_interpolation = " << sync_interpolation << "\n"
                   << "  sync_projection_residuals = " << sync_projection_residuals << "\n"
                   << "  advection_type = " << advection_type << "\n"
                   << "  profile_phases = " << profile_phases << "\n"
                   << "  plot_int = " << plot_int << "\n"
                   << "  plot_file = " << plot_file << "\n"
                   << "  plot_initial = " << plot_initial << "\n"
                   << "  regrid_int = " << regrid_int << "\n"
                   << "  bc.xlo = " << boundary.lo[0].type << " "
                   << boundary.lo[0].velocity[0] << " "
                   << boundary.lo[0].velocity[1] << " "
                   << boundary.lo[0].velocity[2] << "\n"
                   << "  bc.xhi = " << boundary.hi[0].type << " "
                   << boundary.hi[0].velocity[0] << " "
                   << boundary.hi[0].velocity[1] << " "
                   << boundary.hi[0].velocity[2] << "\n"
                   << "  bc.ylo = " << boundary.lo[1].type << " "
                   << boundary.lo[1].velocity[0] << " "
                   << boundary.lo[1].velocity[1] << " "
                   << boundary.lo[1].velocity[2] << "\n"
                   << "  bc.yhi = " << boundary.hi[1].type << " "
                   << boundary.hi[1].velocity[0] << " "
                   << boundary.hi[1].velocity[1] << " "
                   << boundary.hi[1].velocity[2] << "\n"
#if (AMREX_SPACEDIM == 3)
                   << "  bc.zlo = " << boundary.lo[2].type << " "
                   << boundary.lo[2].velocity[0] << " "
                   << boundary.lo[2].velocity[1] << " "
                   << boundary.lo[2].velocity[2] << "\n"
                   << "  bc.zhi = " << boundary.hi[2].type << " "
                   << boundary.hi[2].velocity[0] << " "
                   << boundary.hi[2].velocity[1] << " "
                   << boundary.hi[2].velocity[2] << "\n"
#endif
                   << "  dlm.enabled = " << dlm.enabled << "\n"
                   << "  dlm.scheme = " << dlm.scheme << "\n"
                   << "  dlm.mode = " << dlm.mode << "\n"
                   << "  dlm.interpolation = " << dlm.interpolation << "\n"
                   << "  dlm.fluid_response = " << dlm.fluid_response << "\n"
                   << "  dlm.multiplier_layout = " << dlm.multiplier_layout << "\n"
                   << "  dlm.constraint_output = " << dlm.constraint_output << "\n"
                   << "  dlm.history_file = " << dlm.history_file << "\n"
                   << "  dlm.radius = " << dlm.radius << "\n"
                   << "  dlm.marker_spacing = " << dlm.marker_spacing << "\n"
                   << "  dlm.interior_margin = " << dlm.interior_margin << "\n"
                   << "  dlm.relaxation = " << dlm.relaxation << "\n"
                   << "  dlm.constraint_abs_tol = " << dlm.constraint_abs_tol << "\n"
                   << "  dlm.constraint_rel_tol = " << dlm.constraint_rel_tol << "\n"
                   << "  dlm.fluid_density = " << dlm.fluid_density << "\n"
                   << "  dlm.solid_density = " << dlm.solid_density << "\n"
                   << "  dlm.constraint_iterations = " << dlm.constraint_iterations << "\n"
                   << "  dlm.refinement_max_level = " << dlm.refinement_max_level << "\n"
                   << "  dlm.explicit_added_mass = " << dlm.explicit_added_mass << "\n"
                   << "  dlm.use_boundary_points = " << dlm.use_boundary_points << "\n"
                   << "  dlm.use_interior_points = " << dlm.use_interior_points << "\n"
                   << "  dlm.enforce_interior_velocity = "
                   << dlm.enforce_interior_velocity << "\n"
                   << "  dlm.refine = " << dlm.refine << "\n"
                   << "  dlm.refine_rigid_mask = " << dlm.refine_rigid_mask << "\n"
                   << "  dlm.refine_stencil_support = " << dlm.refine_stencil_support << "\n"
                   << "  dlm.require_finest_stencil_support = "
                   << dlm.require_finest_stencil_support << "\n"
                   << "  dlm.apply_on_finest_level_only = "
                   << dlm.apply_on_finest_level_only << "\n"
                   << "  dlm.fail_on_nonconvergence = "
                   << dlm.fail_on_nonconvergence << "\n"
                   << "  dlm.refinement_buffer = " << dlm.refinement_buffer << "\n"
                   << "  dlm.center = "
                   << dlm.center[0] << " " << dlm.center[1]
#if (AMREX_SPACEDIM == 3)
                   << " " << dlm.center[2]
#endif
                   << "\n"
                   << "  dlm.gravity = "
                   << dlm.gravity[0] << " " << dlm.gravity[1]
#if (AMREX_SPACEDIM == 3)
                   << " " << dlm.gravity[2]
#endif
                   << "\n";
}

} // namespace ELFF::Fluid
