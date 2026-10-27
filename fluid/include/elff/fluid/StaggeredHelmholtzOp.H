#ifndef ELFF_FLUID_STAGGERED_HELMHOLTZ_OP_H_
#define ELFF_FLUID_STAGGERED_HELMHOLTZ_OP_H_

#include <AMReX_MLLinOp.H>
#include <AMReX_REAL.H>

namespace ELFF::Fluid {

/** @brief Scalar face-centered operator alpha * u - beta * Laplacian(u).
 * Supplies AMReX multigrid application, smoothing, transfer, and residual
 * operations for one velocity component. Its use in the DLM response is
 * level-local; it does not imply a composite AMR DLM solve. */
class StaggeredHelmholtzOp final : public amrex::MLLinOp {
public:
  /// Configure geometry, layout, face direction, and constant coefficients.
  StaggeredHelmholtzOp(amrex::Vector<amrex::Geometry> const &geom,
                       amrex::Vector<amrex::BoxArray> const &grids,
                       amrex::Vector<amrex::DistributionMapping> const &dmap,
                       amrex::LPInfo const &info, int face_component,
                       amrex::Real alpha, amrex::Real beta);

  [[nodiscard]] std::string name() const override;
  [[nodiscard]] int getNComp() const override;
  [[nodiscard]] int getNGrow(int a_lev = 0, int mg_lev = 0) const override;
  [[nodiscard]] bool isSingular(int amrlev) const override;
  [[nodiscard]] bool isBottomSingular() const override;

  void setLevelBC(int amrlev, amrex::MultiFab const *levelbcdata,
                  amrex::MultiFab const *robinbc_a = nullptr,
                  amrex::MultiFab const *robinbc_b = nullptr,
                  amrex::MultiFab const *robinbc_f = nullptr) override;
  void restriction(int amrlev, int cmglev, amrex::MultiFab &crse,
                   amrex::MultiFab &fine) const override;
  void interpolation(int amrlev, int fmglev, amrex::MultiFab &fine,
                     amrex::MultiFab const &crse) const override;
  void interpAssign(int amrlev, int fmglev, amrex::MultiFab &fine,
                    amrex::MultiFab &crse) const override;
  void apply(int amrlev, int mglev, amrex::MultiFab &out, amrex::MultiFab &in,
             BCMode bc_mode, StateMode s_mode,
             amrex::MLMGBndry const *bndry = nullptr) const override;
  void smooth(int amrlev, int mglev, amrex::MultiFab &sol,
              amrex::MultiFab const &rhs, bool skip_fillboundary,
              int niter) const override;
  void Fapply(int amrlev, int mglev, amrex::MultiFab &out,
              amrex::MultiFab const &in) const;
  void Fsmooth(int amrlev, int mglev, amrex::MultiFab &sol,
               amrex::MultiFab const &rhs, int redblack) const;
  void normalize(int amrlev, int mglev, amrex::MultiFab &mf) const override;
  void solutionResidual(int amrlev, amrex::MultiFab &resid, amrex::MultiFab &x,
                        amrex::MultiFab const &b,
                        amrex::MultiFab const *crse_bcdata = nullptr) override;
  void
  correctionResidual(int amrlev, int mglev, amrex::MultiFab &resid,
                     amrex::MultiFab &x, amrex::MultiFab const &b,
                     BCMode bc_mode,
                     amrex::MultiFab const *crse_bcdata = nullptr) override;
  void prepareForSolve() override;
  [[nodiscard]] amrex::Real xdoty(int amrlev, int mglev,
                                  amrex::MultiFab const &x,
                                  amrex::MultiFab const &y,
                                  bool local) const override;
  [[nodiscard]] amrex::Real normInf(int amrlev, amrex::MultiFab const &mf,
                                    bool local) const override;
  void averageDownAndSync(amrex::Vector<amrex::MultiFab> &sol) const override;
  void avgDownResMG(int clev, amrex::MultiFab &cres,
                    amrex::MultiFab const &fres) const override;

private:
  void applyBC(int amrlev, int mglev, amrex::MultiFab &mf,
               bool skip_fillboundary = false) const;

  int face_component_ = 0;
  amrex::Real alpha_ = 1.0;
  amrex::Real beta_ = 0.0;
};

} // namespace ELFF::Fluid

#endif
