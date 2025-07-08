#pragma once

// #include <beam/LinAlg/Matrix.hpp>
// #include <beam/LinAlg/Vector.hpp>

#include "EulerBeam.hpp"
#include "Shapes.hpp"

#include <cstdio> // for popen, pclose, fprintf
#include <iostream>
#include <stdlib.h>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/AutoDiff>

namespace beam {

class EulerBeamStaticInextensibleADDM : public EulerBeam
{
public:
  EulerBeamStaticInextensibleADDM(real_t length,
                                  real_t EI,
                                  real_t load,
                                  real_t area,
                                  size_t nodes,
                                  beam::EulerBeamBCs bcs,
                                  real_t r_penalty)
    : EulerBeam(length, EI, area, nodes, bcs)
    , elements(nodes - 1)
    , dimension(2)
    , dof((2 * nodes * dimension))
    , jacobian(Eigen::MatrixXd::Zero(dof, dof))
    , residual(Eigen::VectorXd::Zero(dof))
    , lambda(Eigen::VectorXd::Zero(nodes))
    , u(Eigen::VectorXd::Zero(dof))
    , r_penalty(r_penalty)
  {
    set_load({ 0., load, 0. });
  };

  ~EulerBeamStaticInextensibleADDM() {};

  const Eigen::VectorXd& get_solution() const { return u; }
  const Eigen::VectorXd& get_residual() const { return residual; }
  const Eigen::VectorXd& get_lambda() const { return lambda; }
  const Eigen::MatrixXd& get_jacobian() const { return jacobian; }

  void set_solution(Eigen::VectorXd s) { this->u = s; }
  void set_residual(Eigen::VectorXd r) { this->residual = r; }
  void set_lambda(Eigen::VectorXd l) { this->lambda = l; }
  void set_jacobian(Eigen::MatrixXd j) { this->jacobian = j; }

  void solve()
  {
    real_t S_norm = 0;
    size_t max_iter_inner = 3;
    size_t max_iter_outer = 50;
    real_t tol_inner = 1e-6;
    real_t tol_outer = 1e-6;

    Eigen::LDLT<Eigen::MatrixXd> ldlt;

    Eigen::ConjugateGradient<
      Eigen::MatrixXd,                      // or SparseMatrix<double>
      Eigen::Lower | Eigen::Upper,          // tell it K is symmetric
      Eigen::DiagonalPreconditioner<real_t> // Jacobi preconditioner
      >
      cg;

    apply_initial_condition();

    for (size_t it_o = 0; it_o < max_iter_outer; it_o++) {

      std::cout << "iter " << it_o << "  ‖S‖ = " << S_norm << "\n";

      if (S_norm < tol_outer && it_o > 0) {
        std::cout << "Converged in " << it_o << " iters.\n";
        break;
      }

      for (size_t it_i = 0; it_i < max_iter_inner; it_i++) {

        assemble_system();
        apply_boundary_conditions();

        double res_norm = residual.norm();
        std::cout << "iter " << it_i << "  ‖R‖ = " << res_norm << "\n";
        if (res_norm < tol_inner) {
          std::cout << "Converged in " << it_i << " iters.\n";
          break;
        }

        Eigen::VectorXd delta_u;

        cg.compute(jacobian);
        if (cg.info() != Eigen::Success) {
          throw std::runtime_error("IC preconditioner failed");
        }
        delta_u = cg.solve(-residual);

        // ldlt.compute(jacobian);
        // delta_u = jacobian.colPivHouseholderQr().solve(-residual);

        // 4) update
        u += delta_u;
      }

      S_norm = update_lambda();
    }
    update_mesh();
  }

  void assemble_residual() { residual = assemble_residual_template<real_t>(u); }

  real_t update_lambda()
  {
    size_t nodes = mesh.get_nodes();
    size_t ndof_x = 2 * nodes;
    size_t ndof_y = 2 * nodes;
    size_t ndof_l = nodes;

    size_t offset_x = 0;
    size_t offset_y = ndof_x;

    Eigen::VectorXd S = Eigen::VectorXd::Zero(ndof_l);

    for (size_t i = 0; i < nodes; i++) {
      real_t xp = u[2 * i + offset_x + 1];
      real_t yp = u[2 * i + offset_y + 1];

      S[i] = xp * xp + yp * yp - 1.0;

      lambda[i] += S[i] * r_penalty;
    }

    return S.norm();
  }

  void update_mesh()
  {
    size_t nodes = this->mesh.get_nodes();

    size_t ndof_x = 2 * nodes;
    size_t ndof_y = 2 * nodes;
    size_t ndof_l = nodes;

    size_t offset_x = 0;
    size_t offset_y = ndof_x;
    size_t offset_l = ndof_x + ndof_y;

    std::vector<std::array<real_t, 3>>& centerline = mesh.get_centerline();
    std::vector<real_t>& s = mesh.get_curvilinear_axis();

    for (size_t i = 0; i < nodes; ++i) {
      centerline[i][0] = u(offset_x + 2 * i);
      centerline[i][1] = u(offset_y + 2 * i);
      centerline[i][2] = 0.;
    }
  }

  void assemble_system()
  {
    using AD = Eigen::AutoDiffScalar<Eigen::VectorXd>;
    using ADVec = Eigen::Matrix<AD, Eigen::Dynamic, 1>;

    ADVec x_ad(dof);
    for (int i = 0; i < int(dof); ++i) {
      Eigen::VectorXd seed = Eigen::VectorXd::Zero(dof);
      seed(i) = 1.0;
      x_ad(i) = AD(u(i), seed);
    }

    ADVec R_ad = assemble_residual_template<AD>(x_ad);

    residual.resize(dof);
    jacobian.resize(dof, dof);
    for (int i = 0; i < int(dof); ++i) {
      residual(i) = R_ad(i).value();
      jacobian.row(i) = R_ad(i).derivatives().transpose();
    }
  }

  void apply_initial_condition()
  {
    real_t h = mesh.get_ds();
    size_t nodes = mesh.get_nodes();

    size_t ndof_x = 2 * nodes;
    size_t ndof_y = 2 * nodes;

    size_t offset_x = 0;
    size_t offset_y = ndof_x;

    for (size_t i = 0; i < nodes; i++) {
      u(offset_x + 2 * i + 0) = h * i;
      u(offset_x + 2 * i + 1) = 1.;
      u(offset_y + 2 * i + 0) = 0.;
      u(offset_y + 2 * i + 1) = 0.;
    }
  }

  /**
   * @brief Apply boundary conditions to the residual and jacobian
   *
   * Modifies the residual and jacobian according to boundary conditions:
   * - free_bc: No constraints
   * - simple_bc: Position constraints only
   * - clamped_bc: Position and slope constraints
   */

  void apply_boundary_conditions()
  {
    size_t nodes = mesh.get_nodes();
    real_t h = mesh.get_ds();

    size_t ndof_x = 2 * nodes;
    size_t ndof_y = 2 * nodes;
    size_t ndof_l = nodes;
    size_t offset_x = 0;
    size_t offset_y = ndof_x;
    size_t offset_l = ndof_x + ndof_y;

    for (size_t bi = 0; bi < 2; ++bi) {
      EulerBeamBCEnd bcend = boundary_conditions.end[bi];
      size_t ni = 0;
      std::vector<size_t> idx(4);

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
      std::vector<real_t> vals(4);

      switch (bctype) {
        case free_bc:
          idx = {};
          vals = {};
          break;
        case simple_bc:
          idx = { offset_x + 2 * ni + 0, offset_y + 2 * ni + 0 };
          vals = { bcvals.position[0], bcvals.position[1] };
          break;
        case clamped_bc:
          idx = { offset_x + 2 * ni + 0,
                  offset_x + 2 * ni + 1,
                  offset_y + 2 * ni + 0,
                  offset_y + 2 * ni + 1 };
          vals = { bcvals.position[0],
                   bcvals.slope[0],
                   bcvals.position[1],
                   bcvals.slope[1] };
          break;
      }

      for (size_t i = 0; i < vals.size(); ++i) {
        for (size_t j = 0; j < dof; ++j) {
          jacobian(idx[i], j) = 0.;
          jacobian(j, idx[i]) = 0.;
        }
        residual(idx[i]) = u(idx[i]) - vals[i];
        jacobian(idx[i], idx[i]) = 1.;
      }
    }
  }

  /**
   * @brief Assemble the residual vector for the nonlinear system
   *
   * @tparam T Scalar type (real_t or autodiff)
   * @param u Current solution vector
   * @param lambda Lagrangian multiplier
   * @return Assembled residual vector containing:
   *         - Bending energy terms
   *         - External load contributions
   *         - Inextensibility constraints
   *
   * Uses Gauss quadrature with 3 points for numerical integration.
   */
  template<typename T>
  Eigen::Matrix<T, Eigen::Dynamic, 1> assemble_residual_template(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& u)
  {
    real_t h = mesh.get_ds();
    size_t nodes = mesh.get_nodes();

    real_t xi_q[] = { 0.1127016654, 0.5, 0.8872983346 };
    real_t w_q[] = { 0.2777777778, 0.4444444444, 0.2777777778 };

    size_t ndof_x = 2 * nodes;
    size_t ndof_y = 2 * nodes;

    size_t offset_x = 0;
    size_t offset_y = ndof_x;

    Eigen::Matrix<T, Eigen::Dynamic, 1> residual =
      Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(dof);

    for (size_t e = 0; e < elements; ++e) {
      std::vector<size_t> elem_nodes = { e, e + 1 };
      std::vector<size_t> idx_x = { offset_x + 2 * elem_nodes[0],
                                    offset_x + 2 * elem_nodes[0] + 1,
                                    offset_x + 2 * elem_nodes[1],
                                    offset_x + 2 * elem_nodes[1] + 1 };
      std::vector<size_t> idx_y = { offset_y + 2 * elem_nodes[0],
                                    offset_y + 2 * elem_nodes[0] + 1,
                                    offset_y + 2 * elem_nodes[1],
                                    offset_y + 2 * elem_nodes[1] + 1 };
      std::vector<size_t> idx_l = { elem_nodes[0], elem_nodes[1] };

      std::array<T, 4> ux = {
        u[idx_x[0]], u[idx_x[1]], u[idx_x[2]], u[idx_x[3]]
      };
      std::array<T, 4> uy = {
        u[idx_y[0]], u[idx_y[1]], u[idx_y[2]], u[idx_y[3]]
      };
      std::array<T, 2> ul = { lambda[idx_l[0]], lambda[idx_l[1]] };

      std::vector<T> R_loc_x(4, 0);
      std::vector<T> R_loc_y(4, 0);
      // std::vector<T> R_loc_l(2, 0);

      for (size_t qi = 0; qi < 3; ++qi) {
        real_t xi = xi_q[qi];
        real_t w = w_q[qi];

        auto H = CubicHermite<real_t>::values(xi, h);
        auto dH = CubicHermite<real_t>::derivs(xi, h);
        auto ddH = CubicHermite<real_t>::second_derivs(xi, h);
        auto M = LinearShape<real_t>::values(xi);

        T x = 0, xp = 0, xpp = 0;
        T y = 0, yp = 0, ypp = 0;
        for (size_t i = 0; i < 4; i++) {
          x += H[i] * ux[i];
          xp += dH[i] * ux[i];
          xpp += ddH[i] * ux[i];
          y += H[i] * uy[i];
          yp += dH[i] * uy[i];
          ypp += ddH[i] * uy[i];
        }

        T l = 0;
        for (size_t i = 0; i < 2; i++) {
          l += M[i] * ul[i];
        }

        T S = xp * xp + yp * yp - 1.0;

        for (size_t a = 0; a < 4; ++a) {
          // Bending Energy Contributions
          R_loc_x[a] += EI * xpp * ddH[a] * w * h;
          R_loc_y[a] += EI * ypp * ddH[a] * w * h;
          // External load contribution in y
          R_loc_x[a] -= uniform_load[0] * H[a] * w * h;
          R_loc_y[a] -= uniform_load[1] * H[a] * w * h;
          // Constraint contributions
          R_loc_x[a] += 2 * (l + r_penalty * S) * xp * dH[a] * w * h;
          R_loc_y[a] += 2 * (l + r_penalty * S) * yp * dH[a] * w * h;
        }

        // for (size_t a = 0; a < 2; ++a) {
        //   // Variation w.r.t. lambda
        //   R_loc_l[a] += S * M[a] * w * h;
        // }
      }

      // Scatter local residual contribution to global residual
      for (size_t i = 0; i < 4; ++i) {
        residual[idx_x[i]] += R_loc_x[i];
        residual[idx_y[i]] += R_loc_y[i];
      }

      // for (size_t i = 0; i < 2; ++i) {
      //   residual[idx_l[i]] += R_loc_l[i];
      // }
    }
    return residual;
  }

  size_t dimension;
  size_t elements;
  size_t dof;
  real_t r_penalty;

  Eigen::VectorXd residual, lambda;
  Eigen::MatrixXd jacobian;
  Eigen::VectorXd u;
};

} // namespace beam