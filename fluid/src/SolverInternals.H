#ifndef ELFF_FLUID_SOLVER_INTERNALS_H_
#define ELFF_FLUID_SOLVER_INTERNALS_H_

#include <elff/fluid/DLMTypes.H>
#include <elff/fluid/LevelData.H>
#include <elff/fluid/Options.H>

#include <AMReX_BCRec.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

namespace ELFF::Fluid::SolverDetail {

struct StaggeredResponseWorkspace
{
    DLM::FaceFields scaled_rhs;
    DLM::FaceFields solution;
    DLM::FaceFields beta;
};

amrex::Vector<amrex::BCRec>
makeExtrapolationBCs(amrex::Geometry const& geom, int ncomp);

DLM::RigidSphere optionSphere(DLMOptions const& options);

void projectFaceVelocity(DLM::FaceFields& face_velocity,
                         DLM::FaceFields const& face_beta,
                         amrex::Geometry const& geom,
                         SolverOptions const& options);

void projectFaceVelocity(LevelData& state, SolverOptions const& options);

void extrapolateFaceVelocityFromCell(LevelData& state, SolverOptions const& options);

void refreshFaceVelocityFromCell(LevelData& state, SolverOptions const& options);

void projectVelocityField(amrex::MultiFab& velocity,
                          amrex::Geometry const& geom,
                          SolverOptions const& options,
                          amrex::MultiFab* pressure,
                          amrex::MultiFab* gradp,
                          amrex::MultiFab* sync_resid_crse = nullptr,
                          amrex::MultiFab* sync_resid_fine = nullptr,
                          amrex::IntVect const& sync_ref_ratio = amrex::IntVect(1),
                          amrex::BoxArray const& sync_fine_grids = amrex::BoxArray());

void defineFaceFieldsLike(LevelData const& state, DLM::FaceFields& faces, int ngrow);

void defineStaggeredResponseWorkspace(LevelData const& state,
                                      StaggeredResponseWorkspace& workspace);

void copyFaceFields(DLM::FaceFields& dst, DLM::FaceFields const& src);

void solveStaggeredVelocityResponse(LevelData const& state,
                                    SolverOptions const& options,
                                    DLM::FaceFields const& rhs_delta,
                                    DLM::FaceFields& velocity_delta,
                                    StaggeredResponseWorkspace& workspace,
                                    amrex::Real dt);

} // namespace ELFF::Fluid::SolverDetail

#endif
