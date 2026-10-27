#include <elff/fluid/StaggeredHelmholtzOp.H>

#include <AMReX_GpuLaunch.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFabUtil.H>

#include <algorithm>

namespace ELFF::Fluid {

namespace {

struct FaceStencilGeometry
{
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> inv_dx_sq{};
    amrex::GpuArray<int, AMREX_SPACEDIM> domain_lo{};
    amrex::GpuArray<int, AMREX_SPACEDIM> face_hi{};
    amrex::GpuArray<int, AMREX_SPACEDIM> periodic{};
};

FaceStencilGeometry
makeFaceStencilGeometry(amrex::Geometry const& geom, int face_component)
{
    FaceStencilGeometry stencil_geom;
    const auto* dx = geom.CellSize();
    const amrex::Box& domain = geom.Domain();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        stencil_geom.inv_dx_sq[dir] = amrex::Real(1.0) / (dx[dir] * dx[dir]);
        stencil_geom.domain_lo[dir] = domain.smallEnd(dir);
        stencil_geom.face_hi[dir] = domain.bigEnd(dir)
                                  + ((dir == face_component) ? 1 : 0);
        stencil_geom.periodic[dir] = geom.isPeriodic(dir) ? 1 : 0;
    }
    return stencil_geom;
}

amrex::Dim3
makeDim3(amrex::IntVect const& v)
{
    amrex::Dim3 result{v[0], 1, 1};
#if (AMREX_SPACEDIM >= 2)
    result.y = v[1];
#endif
#if (AMREX_SPACEDIM == 3)
    result.z = v[2];
#endif
    return result;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
helmholtzDiagonal(int i,
                  int j,
                  int k,
                  FaceStencilGeometry const& geom,
                  amrex::Real alpha,
                  amrex::Real beta) noexcept
{
    amrex::ignore_unused(k);
    amrex::GpuArray<int, AMREX_SPACEDIM> idx{AMREX_D_DECL(i, j, k)};
    amrex::Real diag = alpha;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const amrex::Real coeff = beta * geom.inv_dx_sq[dir];
        const bool has_lower = (geom.periodic[dir] != 0)
                            || (idx[dir] > geom.domain_lo[dir]);
        const bool has_upper = (geom.periodic[dir] != 0)
                            || (idx[dir] < geom.face_hi[dir]);
        if (has_lower) {
            diag += coeff;
        }
        if (has_upper) {
            diag += coeff;
        }
    }
    return diag;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
helmholtzNeighborSum(int i,
                     int j,
                     int k,
                     amrex::Array4<amrex::Real const> const& in,
                     FaceStencilGeometry const& geom,
                     amrex::Real beta) noexcept
{
    amrex::GpuArray<int, AMREX_SPACEDIM> idx{AMREX_D_DECL(i, j, k)};
    amrex::Real sum = 0.0;

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const amrex::Real coeff = beta * geom.inv_dx_sq[dir];
        const bool has_lower = (geom.periodic[dir] != 0)
                            || (idx[dir] > geom.domain_lo[dir]);
        const bool has_upper = (geom.periodic[dir] != 0)
                            || (idx[dir] < geom.face_hi[dir]);

        if (has_lower) {
            int ni = i;
            int nj = j;
            int nk = k;
            if (dir == 0) {
                --ni;
            } else if (dir == 1) {
                --nj;
#if (AMREX_SPACEDIM == 3)
            } else {
                --nk;
#endif
            }
            sum += coeff * in(ni, nj, nk);
        }

        if (has_upper) {
            int ni = i;
            int nj = j;
            int nk = k;
            if (dir == 0) {
                ++ni;
            } else if (dir == 1) {
                ++nj;
#if (AMREX_SPACEDIM == 3)
            } else {
                ++nk;
#endif
            }
            sum += coeff * in(ni, nj, nk);
        }
    }

    return sum;
}

} // namespace

StaggeredHelmholtzOp::StaggeredHelmholtzOp(
    amrex::Vector<amrex::Geometry> const& geom,
    amrex::Vector<amrex::BoxArray> const& grids,
    amrex::Vector<amrex::DistributionMapping> const& dmap,
    amrex::LPInfo const& info,
    int face_component,
    amrex::Real alpha,
    amrex::Real beta)
    : face_component_(face_component),
      alpha_(alpha),
      beta_(beta)
{
    AMREX_ALWAYS_ASSERT(face_component_ >= 0);
    AMREX_ALWAYS_ASSERT(face_component_ < AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(!grids.empty());
    AMREX_ALWAYS_ASSERT(grids[0].ixType().cellCentered());
    m_ixtype = amrex::IntVect::TheDimensionVector(face_component_);
    amrex::Vector<amrex::FabFactory<amrex::FArrayBox> const*> factories;
    amrex::MLLinOp::define(geom, grids, dmap, info, factories, false);
}

std::string
StaggeredHelmholtzOp::name() const
{
    return "StaggeredHelmholtzOp";
}

int
StaggeredHelmholtzOp::getNComp() const
{
    return 1;
}

int
StaggeredHelmholtzOp::getNGrow(int, int) const
{
    return 1;
}

bool
StaggeredHelmholtzOp::isSingular(int) const
{
    return false;
}

bool
StaggeredHelmholtzOp::isBottomSingular() const
{
    return false;
}

void
StaggeredHelmholtzOp::setLevelBC(int amrlev,
                                 amrex::MultiFab const* levelbcdata,
                                 amrex::MultiFab const* robinbc_a,
                                 amrex::MultiFab const* robinbc_b,
                                 amrex::MultiFab const* robinbc_f)
{
    amrex::ignore_unused(amrlev, levelbcdata, robinbc_a, robinbc_b, robinbc_f);
}

void
StaggeredHelmholtzOp::restriction(int amrlev,
                                  int cmglev,
                                  amrex::MultiFab& crse,
                                  amrex::MultiFab& fine) const
{
    const amrex::IntVect ratio =
        (amrlev > 0) ? amrex::IntVect(2) : mg_coarsen_ratio_vec[cmglev - 1];
    amrex::average_down_faces(fine, crse, ratio, 0);
}

void
StaggeredHelmholtzOp::interpolation(int amrlev,
                                    int fmglev,
                                    amrex::MultiFab& fine,
                                    amrex::MultiFab const& crse) const
{
    const amrex::IntVect ratio =
        (amrlev > 0) ? amrex::IntVect(2) : mg_coarsen_ratio_vec[fmglev];
    const amrex::Dim3 ratio3 = makeDim3(ratio);

    for (amrex::MFIter mfi(fine, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto const fine_arr = fine.array(mfi);
        auto const crse_arr = crse.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const int ci = amrex::coarsen(i, ratio3.x);
            const int cj = amrex::coarsen(j, ratio3.y);
            const int ck = amrex::coarsen(k, ratio3.z);
            fine_arr(i,j,k) += crse_arr(ci,cj,ck);
        });
    }
}

void
StaggeredHelmholtzOp::interpAssign(int amrlev,
                                   int fmglev,
                                   amrex::MultiFab& fine,
                                   amrex::MultiFab& crse) const
{
    const amrex::Geometry& crse_geom = Geom(amrlev, fmglev + 1);
    const amrex::IntVect ratio =
        (amrlev > 0) ? amrex::IntVect(2) : mg_coarsen_ratio_vec[fmglev];
    const amrex::IntVect ng = crse.nGrowVect();

    amrex::MultiFab crse_on_fine_layout;
    amrex::MultiFab const* crse_for_interp = nullptr;
    if (amrex::isMFIterSafe(crse, fine)) {
        crse.FillBoundary(crse_geom.periodicity());
        crse_for_interp = &crse;
    } else {
        amrex::BoxArray cba = fine.boxArray();
        cba.coarsen(ratio);
        crse_on_fine_layout.define(cba, fine.DistributionMap(), 1, ng);
        crse_on_fine_layout.setVal(0.0);
        crse_on_fine_layout.ParallelCopy(crse,
                                         0,
                                         0,
                                         1,
                                         amrex::IntVect(0),
                                         ng,
                                         crse_geom.periodicity());
        crse_for_interp = &crse_on_fine_layout;
    }

    const amrex::Dim3 ratio3 = makeDim3(ratio);
    for (amrex::MFIter mfi(fine, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto const fine_arr = fine.array(mfi);
        auto const crse_arr = crse_for_interp->const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const int ci = amrex::coarsen(i, ratio3.x);
            const int cj = amrex::coarsen(j, ratio3.y);
            const int ck = amrex::coarsen(k, ratio3.z);
            fine_arr(i,j,k) = crse_arr(ci,cj,ck);
        });
    }
}

void
StaggeredHelmholtzOp::applyBC(int amrlev,
                              int mglev,
                              amrex::MultiFab& mf,
                              bool skip_fillboundary) const
{
    if (!skip_fillboundary) {
        mf.FillBoundary(Geom(amrlev, mglev).periodicity());
    }
}

void
StaggeredHelmholtzOp::apply(int amrlev,
                            int mglev,
                            amrex::MultiFab& out,
                            amrex::MultiFab& in,
                            BCMode bc_mode,
                            StateMode s_mode,
                            amrex::MLMGBndry const* bndry) const
{
    amrex::ignore_unused(bc_mode, s_mode, bndry);
    applyBC(amrlev, mglev, in);
    Fapply(amrlev, mglev, out, in);
}

void
StaggeredHelmholtzOp::smooth(int amrlev,
                             int mglev,
                             amrex::MultiFab& sol,
                             amrex::MultiFab const& rhs,
                             bool skip_fillboundary,
                             int niter) const
{
    for (int iter = 0; iter < niter; ++iter) {
        for (int redblack = 0; redblack < 2; ++redblack) {
            applyBC(amrlev, mglev, sol, skip_fillboundary);
            Fsmooth(amrlev, mglev, sol, rhs, redblack);
            skip_fillboundary = false;
        }
    }
}

void
StaggeredHelmholtzOp::Fapply(int amrlev,
                             int mglev,
                             amrex::MultiFab& out,
                             amrex::MultiFab const& in) const
{
    const FaceStencilGeometry stencil_geom =
        makeFaceStencilGeometry(Geom(amrlev, mglev), face_component_);
    const amrex::Real alpha = alpha_;
    const amrex::Real beta = beta_;

    for (amrex::MFIter mfi(out, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto const in_arr = in.const_array(mfi);
        auto const out_arr = out.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const amrex::Real diag =
                helmholtzDiagonal(i, j, k, stencil_geom, alpha, beta);
            const amrex::Real neighbor_sum =
                helmholtzNeighborSum(i, j, k, in_arr, stencil_geom, beta);
            out_arr(i,j,k) = diag * in_arr(i,j,k) - neighbor_sum;
        });
    }
}

void
StaggeredHelmholtzOp::Fsmooth(int amrlev,
                              int mglev,
                              amrex::MultiFab& sol,
                              amrex::MultiFab const& rhs,
                              int redblack) const
{
    const FaceStencilGeometry stencil_geom =
        makeFaceStencilGeometry(Geom(amrlev, mglev), face_component_);
    const amrex::Real alpha = alpha_;
    const amrex::Real beta = beta_;

    for (amrex::MFIter mfi(sol, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto const rhs_arr = rhs.const_array(mfi);
        auto const sol_arr = sol.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (((i + j + k) & 1) != redblack) {
                return;
            }
            const amrex::Real diag =
                helmholtzDiagonal(i, j, k, stencil_geom, alpha, beta);
            const amrex::Real neighbor_sum =
                helmholtzNeighborSum(i, j, k, sol_arr, stencil_geom, beta);
            sol_arr(i,j,k) = (rhs_arr(i,j,k) + neighbor_sum) / diag;
        });
    }
}

void
StaggeredHelmholtzOp::normalize(int amrlev, int mglev, amrex::MultiFab& mf) const
{
    const FaceStencilGeometry stencil_geom =
        makeFaceStencilGeometry(Geom(amrlev, mglev), face_component_);
    const amrex::Real alpha = alpha_;
    const amrex::Real beta = beta_;

    for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto const mf_arr = mf.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            mf_arr(i,j,k) /= helmholtzDiagonal(i, j, k, stencil_geom, alpha, beta);
        });
    }
}

void
StaggeredHelmholtzOp::solutionResidual(int amrlev,
                                       amrex::MultiFab& resid,
                                       amrex::MultiFab& x,
                                       amrex::MultiFab const& b,
                                       amrex::MultiFab const* crse_bcdata)
{
    amrex::ignore_unused(crse_bcdata);
    apply(amrlev, 0, resid, x, BCMode::Inhomogeneous, StateMode::Solution);
    amrex::MultiFab::Xpay(resid, amrex::Real(-1.0), b, 0, 0, 1, amrex::IntVect(0));
}

void
StaggeredHelmholtzOp::correctionResidual(int amrlev,
                                         int mglev,
                                         amrex::MultiFab& resid,
                                         amrex::MultiFab& x,
                                         amrex::MultiFab const& b,
                                         BCMode bc_mode,
                                         amrex::MultiFab const* crse_bcdata)
{
    amrex::ignore_unused(bc_mode, crse_bcdata);
    apply(amrlev, mglev, resid, x, BCMode::Homogeneous, StateMode::Correction);
    amrex::MultiFab::Xpay(resid, amrex::Real(-1.0), b, 0, 0, 1, amrex::IntVect(0));
}

void
StaggeredHelmholtzOp::prepareForSolve()
{
}

amrex::Real
StaggeredHelmholtzOp::xdoty(int amrlev,
                            int mglev,
                            amrex::MultiFab const& x,
                            amrex::MultiFab const& y,
                            bool local) const
{
    amrex::ignore_unused(amrlev, mglev);
    amrex::Real result =
        amrex::Dot(x, 0, y, 0, 1, amrex::IntVect(0), true);
    if (!local) {
        amrex::ParallelAllReduce::Sum(result, amrex::ParallelContext::CommunicatorSub());
    }
    return result;
}

amrex::Real
StaggeredHelmholtzOp::normInf(int amrlev,
                              amrex::MultiFab const& mf,
                              bool local) const
{
    amrex::ignore_unused(amrlev);
    return mf.norminf(0, 1, amrex::IntVect(0), local);
}

void
StaggeredHelmholtzOp::averageDownAndSync(amrex::Vector<amrex::MultiFab>& sol) const
{
    const int finest_amr_lev = NAMRLevels() - 1;
    sol[finest_amr_lev].FillBoundary(Geom(finest_amr_lev, 0).periodicity());

    for (int fine_lev = finest_amr_lev; fine_lev > 0; --fine_lev) {
        amrex::average_down_faces(sol[fine_lev],
                                  sol[fine_lev - 1],
                                  AMRRefRatioVect(fine_lev - 1),
                                  0);
        sol[fine_lev - 1].FillBoundary(Geom(fine_lev - 1, 0).periodicity());
    }
}

void
StaggeredHelmholtzOp::avgDownResMG(int clev,
                                   amrex::MultiFab& cres,
                                   amrex::MultiFab const& fres) const
{
    amrex::average_down_faces(fres, cres, mg_coarsen_ratio_vec[clev - 1], 0);
}

} // namespace ELFF::Fluid
