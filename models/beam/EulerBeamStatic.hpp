#pragma once

#include <iostream>

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>

#include "EulerBeam.hpp"
#include "Shapes.hpp"
#include "config/config.hpp"

namespace beam {

/**
 * @brief A class to solve the static Euler–Bernoulli beam equation; this
 * model is valid for only small deflections.
 *
 * The strong form of an Euler–Bernoulli beam of length \f(L\f), bending
 * stiffness \f(EI\f), subjected to a distributed transverse load \f(q(x)\f),
 * seeks the transverse displacement
 * \f(w(x)\f) satisfying
 * \f[
 *   EI \,\frac{d^4 w}{d x^4}(x) \;=\; q(x),
 *   \quad x \in (0, L),
 * \f]
 * subject to general boundary conditions of the form
 * \f{align*}{
 *   w(a)        &= w_a,          & w'(a)        &= \theta_a,\\
 *   w(b)        &= w_b,          & w'(b)        &= \theta_b,
 * \f}
 * where \f(a=0\f) and \f(b=L\f) in the canonical case.
 *
 * To derive the weak form, let \f(v(x)\f) be a test function in
 * \f(H^2(a,b)\f) satisfying the essential (Dirichlet) boundary conditions
 * \f(v(a)=0,\,v(b)=0\f) (or the homogeneous counterparts of \f(w(a),w(b)\f)).
 * Multiply the strong form by \f(v(x)\f) and integrate over \f([a,b]\f):
 * \f{align*}{
 *   \int_a^b EI\,w''''(x)\,v(x)\,\mathrm{d}x
 *   &= \int_a^b q(x)\,v(x)\,\mathrm{d}x.
 * \f}
 * Integrate by parts twice to shift derivatives onto \f(v\f), invoking
 * the natural (Neumann) boundary conditions for bending moment and shear:
 * \f{align*}{
 *   \int_a^b EI\,w''(x)\,v''(x)\,\mathrm{d}x
 *   &= \int_a^b q(x)\,v(x)\,\mathrm{d}x
 *     + \Bigl[\,EI\,w'''(x)\,v(x) - EI\,w''(x)\,v'(x)\Bigr]_a^b.
 * \f}
 * Enforcing the specified natural BCs (e.g.\ prescribed bending moment
 * \f(EI\,w''\f) or shear force \f(EI\,w'''\f) at the ends) eliminates the
 * boundary term, yielding the variational problem:
 * \f[
 *   \text{Find } w \in V \subset H^2(a,b)
 *   \text{ such that }
 *   \int_a^b EI\,w''\,v''\,\mathrm{d}x
 *   = \int_a^b q\,v\,\mathrm{d}x
 *   \quad\forall\,v\in V_0.
 * \f]
 *
 * In a finite‐element discretization, choose a subspace \f(V_h =
 * \mathrm{span}\{N_i\})\f) of \f(V\f), with \f(w_h(x)=\sum_j w_j N_j(x)\f). The
 * Galerkin condition
 * \f(a(w_h,N_i)=\ell(N_i)\;\forall i)\f) leads to the linear system
 * \f[
 *   K\,\mathbf w = \mathbf f,
 *   \quad
 *   K_{ij} = \int_a^b EI\,N_j''(x)\,N_i''(x)\,\mathrm{d}x,
 *   \quad
 *   f_i    = \int_a^b q(x)\,N_i(x)\,\mathrm{d}x.
 * \f]
 *
 * Element‐level contributions are computed via local integrals in a reference
 * coordinate \f(\xi\in[-1,1]\f) and assembled into the global stiffness matrix
 * and load vector.
 */

class EulerBeamStatic : public EulerBeam
{
public:
  EulerBeamStatic(real_t length, real_t EI, size_t nodes, EulerBeamBCs bcs)
    : EulerBeam(length, EI, nodes, bcs)
    , elements(nodes - 1)
    , dof(2 * nodes)
    , K_global(Eigen::MatrixXd::Zero(dof, dof))
    , F_global(Eigen::VectorXd::Zero(dof))
    , M_global(Eigen::MatrixXd::Zero(dof, dof))
    , u(Eigen::VectorXd::Zero(dof))
    , residual(Eigen::VectorXd::Zero(dof))
  {
  }

  ~EulerBeamStatic() = default;

  void apply_initial_condition(EulerBeamMesh& bmesh)
  {
    size_t nodes = mesh.get_nodes();
    BEAM_ASSERT(
      nodes == bmesh.get_nodes(),
      "Provided mesh must have same node count as the previous mesh.\n");
    auto centerline = bmesh.get_centerline();
    auto slopes = bmesh.get_slope();
    for (size_t ni = 0; ni < nodes; ni++) {
      u(2 * ni + 0) = centerline[ni][1];
      u(2 * ni + 1) = slopes[ni][1];
    }
  }

  void solve()
  {
    assemble_load({0.,0.,0.});
    assemble_stiffness();
    assemble_residual();
    apply_boundary_conditions();

    Eigen::ConjugateGradient<
      Eigen::MatrixXd,                      // or SparseMatrix<double>
      Eigen::Lower | Eigen::Upper,          // tell it K is symmetric
      Eigen::DiagonalPreconditioner<double> // Jacobi preconditioner
      >
      cg;

    cg.compute(K_global);
    if (cg.info() != Eigen::Success)
      BEAM_ABORT("CG decomposition failed");

    u = cg.solve(F_global);
    if (cg.info() != Eigen::Success)
      BEAM_ABORT("CG did not converge");

    std::cout << "CG iters: " << cg.iterations()
              << ", final error est.: " << cg.error() << "\n";

    update_mesh();
  }

  void solve(std::array<real_t, 3> load)
  {
    assemble_load(load);
    assemble_stiffness();
    assemble_residual();
    apply_boundary_conditions();

    Eigen::ConjugateGradient<
      Eigen::MatrixXd,                      // or SparseMatrix<double>
      Eigen::Lower | Eigen::Upper,          // tell it K is symmetric
      Eigen::DiagonalPreconditioner<double> // Jacobi preconditioner
      >
      cg;

    cg.compute(K_global);
    if (cg.info() != Eigen::Success)
      BEAM_ABORT("CG decomposition failed");

    u = cg.solve(F_global);
    if (cg.info() != Eigen::Success)
      BEAM_ABORT("CG did not converge");

    std::cout << "CG iters: " << cg.iterations()
              << ", final error est.: " << cg.error() << "\n";

    update_mesh();
  }

protected:
  size_t elements, dof;
  Eigen::MatrixXd K_global, M_global;
  Eigen::VectorXd F_global, u, residual;

  EulerBeamStatic(real_t length,
                  real_t EI,
                  real_t mu,
                  size_t nodes,
                  EulerBeamBCs bcs)
    : EulerBeam(length, EI, nodes, bcs)
    , elements(nodes - 1)
    , dof(2 * nodes)
    , K_global(Eigen::MatrixXd::Zero(dof, dof))
    , F_global(Eigen::VectorXd::Zero(dof))
    , M_global(Eigen::MatrixXd::Zero(dof, dof))
    , u(Eigen::VectorXd::Zero(dof))
    , residual(Eigen::VectorXd::Zero(dof))
  {
  }

  void assemble_residual() { residual = F_global - K_global * u; }

  void assemble_load(std::array<real_t, 3> load)
  {
    real_t h = mesh.get_ds();
    // zero out the global load
    F_global.setZero();

    // 3-point Gauss–Legendre
    constexpr std::array<real_t, 3> q_points = { real_t(0.1127016654),
                                                 real_t(0.5),
                                                 real_t(0.8872983346) };
    constexpr std::array<real_t, 3> q_weights = { real_t(0.2777777778),
                                                  real_t(0.4444444444),
                                                  real_t(0.2777777778) };

    // local element load (4 DOFs per beam element)
    Eigen::Matrix<real_t, 4, 1> F_e = Eigen::Matrix<real_t, 4, 1>::Zero();

    // integrate q(x) * H(xi) over the element
    for (int i = 0; i < 3; ++i) {
      auto H_arr = CubicHermite<real_t>::values(q_points[i], h);
      // wrap it in a 4×1 Eigen vector (no copy)
      Eigen::Map<const Eigen::Matrix<real_t, 4, 1>> H_map(H_arr.data());

      // then assemble
      F_e += q_weights[i] * (h * load[1]) * H_map;
    }

    // assemble into the global vector
    for (size_t e = 0; e < elements; ++e) {
      const int idx[4] = { static_cast<int>(2 * (e + 0) + 0),
                           static_cast<int>(2 * (e + 0) + 1),
                           static_cast<int>(2 * (e + 1) + 0),
                           static_cast<int>(2 * (e + 1) + 1) };
      for (int a = 0; a < 4; ++a)
        F_global(idx[a]) += F_e(a);
    }
  }

  void assemble_stiffness()
  {
    real_t h = mesh.get_ds();

    // zero out the global stiffness
    K_global.setZero();

    // 3-point Gauss–Legendre on [0,1]
    constexpr std::array<real_t, 3> q_points = { real_t(0.1127016654),
                                                 real_t(0.5),
                                                 real_t(0.8872983346) };
    constexpr std::array<real_t, 3> q_weights = { real_t(0.2777777778),
                                                  real_t(0.4444444444),
                                                  real_t(0.2777777778) };

    // element stiffness (4×4) accumulator
    Eigen::Matrix<real_t, 4, 4> K_e = Eigen::Matrix<real_t, 4, 4>::Zero();

    // integrate EI/h³ ∫ (d²H/dxi² ⊗ d²H/dxi²) dxi
    for (int q = 0; q < 3; ++q) {
      auto d2H_arr = CubicHermite<real_t>::second_derivs(q_points[q], h);
      Eigen::Map<const Eigen::Matrix<real_t, 4, 1>> d2H_map(d2H_arr.data());

      // then accumulate
      const real_t factor = (EI)*q_weights[q] * (h);
      // outer product yields a 4×4 matrix
      K_e.noalias() += factor * (d2H_map * d2H_map.transpose());
    }

    // scatter into global K
    for (size_t e = 0; e < elements; ++e) {
      const int idx[4] = {
        int(2 * e), int(2 * e + 1), int(2 * (e + 1)), int(2 * (e + 1) + 1)
      };
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
          K_global(idx[i], idx[j]) += K_e(i, j);
    }
  }

  void assemble_mass_matrix()
  {
    const real_t h = mesh.get_ds();
    M_global.setZero();

    // 3-point Gauss–Legendre on [0,1]
    constexpr std::array<real_t, 3> q_points = { real_t(0.1127016654),
                                                 real_t(0.5),
                                                 real_t(0.8872983346) };
    constexpr std::array<real_t, 3> q_weights = { real_t(0.2777777778),
                                                  real_t(0.4444444444),
                                                  real_t(0.2777777778) };

    // Local element matrix (4×4 for cubic Hermite)
    for (size_t e = 0; e < elements; ++e) {
      Eigen::Matrix<real_t, 4, 4> M_e = Eigen::Matrix<real_t, 4, 4>::Zero();

      // Quadrature
      for (int qp = 0; qp < 3; ++qp) {
        auto H_arr = CubicHermite<real_t>::values(q_points[qp], h);
        Eigen::Map<const Eigen::Matrix<real_t, 4, 1>> H_map(H_arr.data());

        // accumulate: ∫ N_i N_j dx ≈ ∑ h * w_q * H * H^T
        M_e.noalias() += (h * q_weights[qp]) * (H_map * H_map.transpose());
      }

      // scatter into global M_global
      const int idx[4] = {
        int(2 * e), int(2 * e + 1), int(2 * (e + 1)), int(2 * (e + 1) + 1)
      };
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
          M_global(idx[i], idx[j]) += M_e(i, j);
    }
  }

  void apply_boundary_conditions()
  {
    size_t nodes = mesh.get_nodes();

    size_t ndof_y = 2 * nodes;
    size_t offset_y = 0;

    for (size_t bi = 0; bi < 2; ++bi) {
      EulerBeamBCEnd bcend = boundary_conditions.end[bi];
      size_t ni = 0;
      std::vector<size_t> idx(2);

      switch (bcend) {
        case left:
          ni = 0;
          break;
        case right:
          ni = nodes - 1;
          break;
      }

      EulerBeamBCType bctype = boundary_conditions.type[bi];
      EulerBeamBCVals bcvals = boundary_conditions.vals[bi];
      std::vector<double> vals(2);

      switch (bctype) {
        case free_bc:
          idx = {};
          vals = {};
          break;
        case simple_bc:
          idx = { offset_y + 2 * ni + 0 };
          vals = { bcvals.position[1] };
          break;
        case clamped_bc:
          idx = { offset_y + 2 * ni + 0, offset_y + 2 * ni + 1 };
          vals = { bcvals.position[1], bcvals.slope[1] };
          break;
        case point_force_bc:
          idx = { offset_y + 2 * ni + 0 };
          vals = { bcend == left ? -bcvals.force[1] : bcvals.force[1] };
          break;
        case point_torque_bc:
          idx = { offset_y + 2 * ni + 1 };
          vals = { bcend == left ? bcvals.torque[1] : -bcvals.torque[1] };
          break;
      }

      switch (bctype) {
        case point_force_bc:
          for (size_t i = 0; i < vals.size(); ++i) {
            F_global(idx[i]) += vals[i];
          }
          break;
        case point_torque_bc:
          for (size_t i = 0; i < vals.size(); ++i) {
            F_global(idx[i]) += vals[i];
          }
          break;
        default:
          for (size_t i = 0; i < vals.size(); ++i) {
            K_global.row(idx[i]).setZero();
            K_global.col(idx[i]).setZero();
            K_global(idx[i], idx[i]) = 1.;
            F_global(idx[i]) = vals[i];
          }
          break;
      }
    }
  }

  void update_mesh()
  {
    size_t nodes = mesh.get_nodes();
    std::vector<std::array<real_t, 3>>& centerline = mesh.get_centerline();
    std::vector<real_t>& s = mesh.get_curvilinear_axis();
    for (size_t i = 0; i < nodes; ++i) {
      centerline[i][0] = s[i];
      centerline[i][1] = u(2 * i);
      centerline[i][2] = 0.;
    }
  }

};

}
