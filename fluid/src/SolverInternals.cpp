#include "SolverInternals.H"

#include <elff/fluid/Boundary.H>
#include <elff/fluid/StaggeredHelmholtzOp.H>

#include <AMReX_BC_TYPES.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFabUtil.H>

#include <hydro_MacProjector.H>
#include <hydro_NodalProjector.H>
#include <hydro_mol.H>

#include <algorithm>

namespace ELFF::Fluid::SolverDetail {

amrex::Vector<amrex::BCRec>
makeExtrapolationBCs(amrex::Geometry const& geom, int ncomp)
{
    amrex::Vector<amrex::BCRec> bcs(ncomp);
    for (int comp = 0; comp < ncomp; ++comp) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            const int bc_type = geom.isPeriodic(dir)
                ? amrex::BCType::int_dir
                : amrex::BCType::foextrap;
            bcs[comp].setLo(dir, bc_type);
            bcs[comp].setHi(dir, bc_type);
        }
    }
    return bcs;
}

DLM::RigidSphere
optionSphere(DLMOptions const& options)
{
    DLM::RigidSphere sphere;
    sphere.id = 1;
    sphere.geometry_id = 0;
    sphere.particle = true;
    sphere.prescribed = (options.mode == "prescribed");
    sphere.radius = options.radius;
    sphere.density = options.solid_density;
    sphere.center = options.center;
    sphere.velocity = options.velocity;
    sphere.angular_velocity = options.angular_velocity;
    return sphere;
}

void
projectFaceVelocity(DLM::FaceFields& face_velocity,
                    DLM::FaceFields const& face_beta,
                    amrex::Geometry const& geom,
                    SolverOptions const& options)
{
    BL_PROFILE("ELFF::Fluid::projectFaceVelocity");
    Boundary::fillFaceVelocityPhysicalBoundary(face_velocity, geom, options);

    amrex::LPInfo lp_info;
    Hydro::MacProjector mac_projector({amrex::GetArrOfPtrs(face_velocity)},
                                      {amrex::GetArrOfConstPtrs(face_beta)},
                                      {geom},
                                      lp_info);
    auto [lo, hi] = Boundary::projectionBCs(geom, options);
    mac_projector.setDomainBC(lo, hi);
    mac_projector.setVerbose(std::max(0, options.verbose - 1));
    mac_projector.project(options.projection_rel_tol, options.projection_abs_tol);
    Boundary::fillFaceVelocityPhysicalBoundary(face_velocity, geom, options);
}

void
projectFaceVelocity(LevelData& state, SolverOptions const& options)
{
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        state.faceBeta()[dir].setVal(1.0);
    }

    projectFaceVelocity(state.faceVelocity(), state.faceBeta(), state.geom(), options);

    amrex::computeDivergence(state.divu(), amrex::GetArrOfConstPtrs(state.faceVelocity()),
                             state.geom());
}

void
extrapolateFaceVelocityFromCell(LevelData& state, SolverOptions const& options)
{
    Boundary::fillVelocityPhysicalBoundary(state, options);

    auto velocity_bcs = Boundary::makeVelocityBCs(state.geom(), options);
    auto device_velocity_bcs = Boundary::makeDeviceBCs(velocity_bcs);
    MOL::ExtrapVelToFaces(state.velocity(),
                          AMREX_D_DECL(state.faceVelocity()[0],
                                       state.faceVelocity()[1],
                                       state.faceVelocity()[2]),
                          state.geom(),
                          velocity_bcs,
                          device_velocity_bcs.dataPtr());
}

void
refreshFaceVelocityFromCell(LevelData& state, SolverOptions const& options)
{
    extrapolateFaceVelocityFromCell(state, options);
    projectFaceVelocity(state, options);
}

void
projectVelocityField(amrex::MultiFab& velocity,
                     amrex::Geometry const& geom,
                     SolverOptions const& options,
                     amrex::MultiFab* pressure,
                     amrex::MultiFab* gradp,
                     amrex::MultiFab* sync_resid_crse,
                     amrex::MultiFab* sync_resid_fine,
                     amrex::IntVect const& sync_ref_ratio,
                     amrex::BoxArray const& sync_fine_grids)
{
    BL_PROFILE("ELFF::Fluid::projectVelocityField");
    Boundary::fillVelocityPhysicalBoundary(velocity, geom, options);

    amrex::LPInfo lp_info;
    Hydro::NodalProjector nodal_projector({&velocity},
                                          amrex::Real(1.0),
                                          {geom},
                                          lp_info);
    auto [lo, hi] = Boundary::projectionBCs(geom, options);
    nodal_projector.setDomainBC(lo, hi);
    nodal_projector.setVerbose(std::max(0, options.verbose - 1));
    if (sync_resid_fine != nullptr) {
        nodal_projector.setSyncResidualFine(sync_resid_fine);
    }
    if (sync_resid_crse != nullptr) {
        nodal_projector.setSyncResidualCrse(sync_resid_crse,
                                           sync_ref_ratio,
                                           sync_fine_grids);
    }
    nodal_projector.project(options.projection_rel_tol, options.projection_abs_tol);

    auto phi = nodal_projector.getPhi();
    if (pressure != nullptr && !phi.empty()) {
        pressure->ParallelCopy(*phi[0], 0, 0, 1, 0, 0, geom.periodicity());
    }

    auto grad_phi = nodal_projector.getGradPhi();
    if (gradp != nullptr && !grad_phi.empty()) {
        amrex::MultiFab::Copy(*gradp, *grad_phi[0], 0, 0, AMREX_SPACEDIM, 0);
    }

    Boundary::fillVelocityPhysicalBoundary(velocity, geom, options);
}

void
defineFaceFieldsLike(LevelData const& state, DLM::FaceFields& faces, int ngrow)
{
    for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
        faces[comp].define(state.faceVelocity()[comp].boxArray(),
                           state.faceVelocity()[comp].DistributionMap(),
                           1,
                           ngrow,
                           amrex::MFInfo(),
                           state.faceVelocity()[comp].Factory());
    }
}

void
defineStaggeredResponseWorkspace(LevelData const& state,
                                 StaggeredResponseWorkspace& workspace)
{
    defineFaceFieldsLike(state, workspace.scaled_rhs, 1);
    defineFaceFieldsLike(state, workspace.solution, 1);
    defineFaceFieldsLike(state, workspace.beta, 0);

    for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
        workspace.beta[comp].setVal(1.0);
    }
}

void
copyFaceFields(DLM::FaceFields& dst, DLM::FaceFields const& src)
{
    for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
        dst[comp].setVal(0.0, 0, 1, dst[comp].nGrowVect());
        amrex::MultiFab::Copy(dst[comp], src[comp], 0, 0, 1, 0);
    }
}

namespace {

void
scaleFaceFields(DLM::FaceFields& fields, amrex::Real scale)
{
    for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
        fields[comp].mult(scale, 0, 1, 0);
    }
}

void
solveFaceHelmholtz(LevelData const& state,
                   SolverOptions const& options,
                   DLM::FaceFields const& rhs_delta,
                   StaggeredResponseWorkspace& workspace,
                   amrex::Real dt)
{
    BL_PROFILE("ELFF::Fluid::solveFaceHelmholtz");
    copyFaceFields(workspace.scaled_rhs, rhs_delta);
    const amrex::Real fluid_density =
        std::max(options.dlm.fluid_density, amrex::Real(1.0e-30));
    scaleFaceFields(workspace.scaled_rhs, dt / fluid_density);
    copyFaceFields(workspace.solution, workspace.scaled_rhs);

    const amrex::Real nu = options.kinematic_viscosity;
    const amrex::Real helmholtz_scale = dt * nu;
    if (helmholtz_scale <= 0.0) {
        return;
    }

    auto [lo, hi] = Boundary::velocityOperatorBCs(state.geom(), options);
    for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
        amrex::LPInfo lp_info;
        StaggeredHelmholtzOp op({state.geom()},
                                {state.velocity().boxArray()},
                                {state.velocity().DistributionMap()},
                                lp_info,
                                comp,
                                amrex::Real(1.0),
                                helmholtz_scale);
        op.setDomainBC(lo, hi);
        op.setLevelBC(0, nullptr);

        amrex::MLMG mlmg(op);
        mlmg.setVerbose(std::max(0, options.verbose - 1));
        mlmg.setFinalFillBC(true);
        mlmg.solve({&workspace.solution[comp]},
                   {&workspace.scaled_rhs[comp]},
                   options.viscosity_rel_tol,
                   options.viscosity_abs_tol);
    }
}

} // namespace

void
solveStaggeredVelocityResponse(LevelData const& state,
                               SolverOptions const& options,
                               DLM::FaceFields const& rhs_delta,
                               DLM::FaceFields& velocity_delta,
                               StaggeredResponseWorkspace& workspace,
                               amrex::Real dt)
{
    BL_PROFILE("ELFF::Fluid::solveStaggeredVelocityResponse");
    if (options.dlm.fluid_response == "mass") {
        copyFaceFields(velocity_delta, rhs_delta);
        const amrex::Real fluid_density =
            std::max(options.dlm.fluid_density, amrex::Real(1.0e-30));
        scaleFaceFields(velocity_delta, dt / fluid_density);
        return;
    }
    if (options.dlm.fluid_response != "projected_helmholtz") {
        amrex::Abort(
            "dlm.fluid_response must be either 'mass' or 'projected_helmholtz'");
    }
    solveFaceHelmholtz(state, options, rhs_delta, workspace, dt);
    projectFaceVelocity(workspace.solution, workspace.beta, state.geom(), options);
    copyFaceFields(velocity_delta, workspace.solution);
}

} // namespace ELFF::Fluid::SolverDetail
