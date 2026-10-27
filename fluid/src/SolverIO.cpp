#include <elff/fluid/Solver.H>

#include <AMReX_GpuLaunch.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <string>

namespace ELFF::Fluid {

namespace {

amrex::Vector<std::string>
plotVariableNames()
{
    amrex::Vector<std::string> names;
    names.reserve(6 * AMREX_SPACEDIM + 6);

    names.push_back("x_velocity");
    names.push_back("y_velocity");
#if (AMREX_SPACEDIM == 3)
    names.push_back("z_velocity");
#endif
    names.push_back("density");
    names.push_back("pressure");
    names.push_back("gradp_x");
    names.push_back("gradp_y");
#if (AMREX_SPACEDIM == 3)
    names.push_back("gradp_z");
#endif
    names.push_back("divu");
    names.push_back("mac_velocity_x");
    names.push_back("mac_velocity_y");
#if (AMREX_SPACEDIM == 3)
    names.push_back("mac_velocity_z");
#endif
    names.push_back("kinetic_energy");
    names.push_back("rigid_mask");
    names.push_back("rigid_velocity_x");
    names.push_back("rigid_velocity_y");
#if (AMREX_SPACEDIM == 3)
    names.push_back("rigid_velocity_z");
#endif
    names.push_back("dlm_force_x");
    names.push_back("dlm_force_y");
#if (AMREX_SPACEDIM == 3)
    names.push_back("dlm_force_z");
#endif
    names.push_back("projection_phi");
    names.push_back("grad_phi_x");
    names.push_back("grad_phi_y");
#if (AMREX_SPACEDIM == 3)
    names.push_back("grad_phi_z");
#endif

    return names;
}

void
fillKineticEnergy(amrex::MultiFab& plot_mf, int rho_comp, int ke_comp)
{
    for (amrex::MFIter mfi(plot_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto const plot = plot_mf.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            amrex::Real velocity_squared = 0.0;
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                velocity_squared += plot(i,j,k,dir) * plot(i,j,k,dir);
            }
            plot(i,j,k,ke_comp) = amrex::Real(0.5) * plot(i,j,k,rho_comp) * velocity_squared;
        });
    }
}

amrex::MultiFab
makePlotMultiFab(LevelData const& state)
{
    constexpr int nplot = 6 * AMREX_SPACEDIM + 6;
    constexpr int velocity_comp = 0;
    constexpr int density_comp = velocity_comp + AMREX_SPACEDIM;
    constexpr int pressure_comp = density_comp + 1;
    constexpr int gradp_comp = pressure_comp + 1;
    constexpr int divu_comp = gradp_comp + AMREX_SPACEDIM;
    constexpr int mac_velocity_comp = divu_comp + 1;
    constexpr int kinetic_energy_comp = mac_velocity_comp + AMREX_SPACEDIM;
    constexpr int rigid_mask_comp = kinetic_energy_comp + 1;
    constexpr int rigid_velocity_comp = rigid_mask_comp + 1;
    constexpr int dlm_force_comp = rigid_velocity_comp + AMREX_SPACEDIM;
    constexpr int projection_phi_comp = dlm_force_comp + AMREX_SPACEDIM;
    constexpr int grad_phi_comp = projection_phi_comp + 1;

    amrex::MultiFab plot_mf(state.velocity().boxArray(),
                            state.velocity().DistributionMap(),
                            nplot,
                            0,
                            amrex::MFInfo(),
                            state.velocity().Factory());
    plot_mf.setVal(0.0);

    amrex::MultiFab::Copy(plot_mf, state.velocity(), 0, velocity_comp,
                          AMREX_SPACEDIM, 0);
    amrex::MultiFab::Copy(plot_mf, state.density(), 0, density_comp, 1, 0);
    amrex::average_node_to_cellcenter(plot_mf, pressure_comp, state.pressure(), 0, 1);
    amrex::MultiFab::Copy(plot_mf, state.gradp(), 0, gradp_comp, AMREX_SPACEDIM, 0);
    amrex::MultiFab::Copy(plot_mf, state.divu(), 0, divu_comp, 1, 0);
    amrex::average_face_to_cellcenter(plot_mf, mac_velocity_comp,
                                      amrex::GetArrOfConstPtrs(state.faceVelocity()));
    fillKineticEnergy(plot_mf, density_comp, kinetic_energy_comp);
    amrex::MultiFab::Copy(plot_mf, state.rigidMask(), 0, rigid_mask_comp, 1, 0);
    amrex::MultiFab::Copy(plot_mf, state.rigidVelocity(), 0, rigid_velocity_comp,
                          AMREX_SPACEDIM, 0);
    amrex::MultiFab::Copy(plot_mf, state.dlmForce(), 0, dlm_force_comp,
                          AMREX_SPACEDIM, 0);
    amrex::MultiFab::Copy(plot_mf, plot_mf, pressure_comp, projection_phi_comp, 1, 0);
    amrex::MultiFab::Copy(plot_mf, plot_mf, gradp_comp, grad_phi_comp,
                          AMREX_SPACEDIM, 0);

    const amrex::Real dt = state.newTime() - state.oldTime();
    const amrex::Real inv_dt = dt > 0.0 ? amrex::Real(1.0) / dt : amrex::Real(0.0);
    for (amrex::MFIter mfi(plot_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto const plot = plot_mf.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const amrex::Real density = plot(i,j,k,density_comp);
            const amrex::Real pressure_scale = density * inv_dt;
            plot(i,j,k,pressure_comp) *= pressure_scale;
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                plot(i,j,k,gradp_comp + dir) *= pressure_scale;
                plot(i,j,k,dlm_force_comp + dir) *= density;
            }
        });
    }

    return plot_mf;
}

} // namespace

std::string
Solver::plotFileName() const
{
    return amrex::Concatenate(options_.plot_file, step_, 5);
}

void
Solver::writePlotFile() const
{
    const int nlevels = finestLevel() + 1;
    AMREX_ALWAYS_ASSERT(nlevels > 0);

    amrex::Vector<amrex::MultiFab> plot_mfs(nlevels);
    amrex::Vector<amrex::Geometry> geoms(nlevels);
    amrex::Vector<int> level_steps(nlevels, step_);
    amrex::Vector<amrex::IntVect> ref_ratios;
    ref_ratios.reserve(std::max(0, nlevels - 1));

    for (int lev = 0; lev < nlevels; ++lev) {
        AMREX_ALWAYS_ASSERT(lev < static_cast<int>(levels_.size()));
        AMREX_ALWAYS_ASSERT(levels_[lev] != nullptr);

        plot_mfs[lev] = makePlotMultiFab(*levels_[lev]);
        geoms[lev] = Geom(lev);
        if (lev < nlevels - 1) {
            ref_ratios.push_back(refRatio(lev));
        }
    }

    const auto filename = plotFileName();
    amrex::Print() << "Writing plotfile " << filename << "\n";

    amrex::WriteMultiLevelPlotfile(filename,
                                   nlevels,
                                   amrex::GetVecOfConstPtrs(plot_mfs),
                                   plotVariableNames(),
                                   geoms,
                                   time_,
                                   level_steps,
                                   ref_ratios,
                                   "HyperCLaw-V1.1");
}

void
Solver::writeDlmHistory(amrex::Real output_time) const
{
    if (options_.dlm.history_file.empty() ||
        !amrex::ParallelDescriptor::IOProcessor())
    {
        return;
    }

    const DLMOptions dlm_options = dlmOptionsForCurrentStep();
    const int source_lev = dlm_options.apply_on_finest_level_only ? finestLevel() : 0;
    if (source_lev < 0 || source_lev >= static_cast<int>(dlm_states_.size())) {
        return;
    }
    DLM::State const& dlm_state = dlm_states_[source_lev];
    if (dlm_state.bodies.empty()) {
        return;
    }

    const bool write_header = step_ == 0;
    std::ofstream output(
        options_.dlm.history_file,
        std::ios::out | (write_header ? std::ios::trunc : std::ios::app));
    if (!output) {
        amrex::Abort("Unable to open DLM history file: " + options_.dlm.history_file);
    }
    if (write_header) {
        output << "time,step,id,position_x,position_y,position_z,"
                  "velocity_x,velocity_y,velocity_z,omega_x,omega_y,omega_z,"
                  "multiplier_force_x,multiplier_force_y,multiplier_force_z,"
                  "fictitious_force_x,fictitious_force_y,fictitious_force_z,"
                  "hydro_force_x,hydro_force_y,hydro_force_z,"
                  "multiplier_torque_x,multiplier_torque_y,multiplier_torque_z,"
                  "fictitious_torque_x,fictitious_torque_y,fictitious_torque_z,"
                  "hydro_torque_x,hydro_torque_y,hydro_torque_z,iterations,"
                  "algebraic_initial_l2,algebraic_final_l2,realized_l2,"
                  "realized_relative_l2,realized_linf,converged,termination\n";
    }
    output << std::setprecision(17);
    for (DLM::BodyState const& body : dlm_state.bodies) {
        output << output_time << ',' << step_ + 1 << ',' << body.id;
        for (amrex::Real value : body.center) {
            output << ',' << value;
        }
        for (amrex::Real value : body.velocity) {
            output << ',' << value;
        }
        for (amrex::Real value : body.angular_velocity) {
            output << ',' << value;
        }
        for (amrex::Real value : body.multiplier_force) {
            output << ',' << value;
        }
        for (amrex::Real value : body.fictitious_force) {
            output << ',' << value;
        }
        for (amrex::Real value : body.hydrodynamic_force) {
            output << ',' << value;
        }
        for (amrex::Real value : body.multiplier_torque) {
            output << ',' << value;
        }
        for (amrex::Real value : body.fictitious_torque) {
            output << ',' << value;
        }
        for (amrex::Real value : body.hydrodynamic_torque) {
            output << ',' << value;
        }
        output << ',' << dlm_state.iterations_performed
               << ',' << dlm_state.initial_residual_norm
               << ',' << dlm_state.final_residual_norm
               << ',' << dlm_state.realized_residual_norm
               << ',' << dlm_state.realized_relative_residual
               << ',' << dlm_state.realized_residual_linf
               << ',' << static_cast<int>(dlm_state.converged)
               << ',' << dlm_state.termination_reason << '\n';
    }
}

} // namespace ELFF::Fluid
