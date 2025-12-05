#pragma once

// #include <beam/LinAlg/Matrix.hpp>
// #include <beam/LinAlg/Vector.hpp>

#include "EulerBeam.hpp"
#include "Shapes.hpp"

#include <cmath>
#include <cstdio> // for popen, pclose, fprintf
#include <iostream>
#include <stdlib.h>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky> // for SimplicialLLT
#include <unsupported/Eigen/AutoDiff>

namespace beam {

class EulerBeamStaticInextensibleMoM : public EulerBeam
{
public:
  EulerBeamStaticInextensibleMoM(real_t length,
                                 real_t EI,
                                 size_t nodes,
                                 beam::EulerBeamBCs bcs,
                                 real_t r_penalty)
    : EulerBeam(length, EI, nodes, bcs)
    , elements(nodes - 1)
    , nodes(nodes)
    , ds(mesh.get_ds())
    , ndof_x(2 * nodes)
    , ndof_y(2 * nodes)
    , ndof_z(2 * nodes)
    , ndof_l(1 * nodes)
    , ndof(ndof_x + ndof_y + ndof_z)
    , offset_x(0)
    , offset_y(ndof_x)
    , offset_z(ndof_x + ndof_y)
    , offset_l(ndof_x + ndof_y + ndof_z)
    , r_penalty(r_penalty)
    , max_iter_outer(1000)
    , max_iter_inner(1000)
    , tol_outer(1e-5)
    , tol_inner(1e-5)
    , jacobian(Eigen::MatrixXd::Zero(ndof, ndof))
    , residual(Eigen::VectorXd::Zero(ndof))
    , u(Eigen::VectorXd::Zero(ndof))
    , lambda(Eigen::VectorXd::Zero(ndof_l))
  {
    apply_initial_condition(this->mesh);
  };

  ~EulerBeamStaticInextensibleMoM() {};

  /**
   *
   */
  virtual void solve() override { solve({ 0., 0., 0. }); }

  /**
   *
   */

  void solve(std::array<real_t, 3> load) override
  {
    real_t S_norm = 0;

    Eigen::LDLT<Eigen::MatrixXd> solver;
    // Eigen::LLT<Eigen::MatrixXd> solver;
    // Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solver;
    // Eigen::BiCGSTAB<Eigen::MatrixXd> solver;
    // Eigen::ConjugateGradient<Eigen::MatrixXd, Eigen::Upper | Eigen::Lower> solver; solver.setTolerance(tol_inner);

    for (size_t iter_outer = 0; iter_outer < max_iter_outer; iter_outer++) {
      assemble_system(load);
      apply_boundary_conditions();

      real_t res_norm = residual.norm();

      if (res_norm < tol_outer) {
        std::cout << iter_outer << ": ||r|| = " << res_norm;
        std::cout << "\t ||S|| = " << S_norm << std::endl;
        break;
      } else if (iter_outer == max_iter_outer - 1) {
        BEAM_ABORT(
          "EulerBeamStaticInexntensibleMoM::solve() did not converge.\n");
      } else {
        std::cout << iter_outer << ": ||r|| = " << res_norm;
        std::cout << "\t ||S|| = " << S_norm << std::endl;
      }

      solver.compute(jacobian);

      if (solver.info() != Eigen::Success) {
        BEAM_ABORT("EulerBeamStaticInextensibleMoM::solve(): "
                   "Preconditioner failed.\n");
      }

      Eigen::VectorXd delta_u = solver.solve(-residual);
      u += delta_u;

      S_norm = update_lambda();
    }

    update_mesh();
  }

  void solve(std::vector<std::array<real_t,3>> load)
  {
    real_t S_norm = 0;

    Eigen::LDLT<Eigen::MatrixXd> solver;
    //Eigen::LLT<Eigen::MatrixXd> solver;
    // Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solver;
    // Eigen::BiCGSTAB<Eigen::MatrixXd> solver;
    //Eigen::ConjugateGradient<Eigen::MatrixXd, Eigen::Upper | Eigen::Lower> solver; solver.setTolerance(tol_inner);

    for (size_t iter_outer = 0; iter_outer < max_iter_outer; iter_outer++) {
      assemble_system(load);
      apply_boundary_conditions();

      real_t res_norm = residual.norm();

      if (res_norm < tol_outer) {
        std::cout << iter_outer << ": ||r|| = " << res_norm;
        std::cout << "\t ||S|| = " << S_norm << std::endl;
        break;
      } else if (iter_outer == max_iter_outer - 1) {
        BEAM_ABORT(
          "EulerBeamStaticInexntensibleMoM::solve() did not converge.\n");
      } else {
        std::cout << iter_outer << ": ||r|| = " << res_norm;
        std::cout << "\t ||S|| = " << S_norm << std::endl;
      }

      solver.compute(jacobian);

      if (solver.info() != Eigen::Success) {
        BEAM_ABORT("EulerBeamStaticInextensibleMoM::solve(): "
                   "Preconditioner failed.\n");
      }

      Eigen::VectorXd delta_u = solver.solve(-residual);
      u += delta_u;

      S_norm = update_lambda();
    }

    update_mesh();
  }

  virtual void apply_initial_condition()
  {
    for (size_t i = 0; i < nodes; i++) {
      u(offset_x + 2 * i + 0) = ds * i;
      u(offset_x + 2 * i + 1) = 1.;
      u(offset_y + 2 * i + 0) = 0.;
      u(offset_y + 2 * i + 1) = 0.;
      u(offset_z + 2 * i + 0) = 0.;
      u(offset_z + 2 * i + 1) = 0.;
      lambda(i) = 0.;
    }
    update_mesh();
  }

  virtual void apply_initial_condition(EulerBeamMesh& bmesh) override
  {
    if (bmesh.get_nodes() != mesh.get_nodes()) {
      BEAM_ABORT("");
    }

    auto centerline = bmesh.get_centerline();
    auto slope = bmesh.get_slope();

    for (size_t i = 0; i < nodes; i++) {
      u(offset_x + 2 * i + 0) = centerline[i][0];
      u(offset_y + 2 * i + 0) = centerline[i][1];
      u(offset_z + 2 * i + 0) = centerline[i][2];
      u(offset_x + 2 * i + 1) = slope[i][0];
      u(offset_y + 2 * i + 1) = slope[i][1];
      u(offset_z + 2 * i + 1) = slope[i][2];

      lambda(i) = 0.;
    }

    update_mesh();
  }

protected:
  size_t dimension;
  size_t elements, nodes;
  real_t ds;
  size_t ndof_x, ndof_y, ndof_z, ndof_l;
  size_t offset_x, offset_y, offset_z, offset_l;
  size_t ndof;
  real_t r_penalty;
  size_t max_iter_outer, max_iter_inner;
  real_t tol_outer, tol_inner;

  Eigen::VectorXd residual;
  Eigen::MatrixXd jacobian, mass;
  Eigen::VectorXd u, lambda;

  EulerBeamStaticInextensibleMoM(real_t length,
                                 real_t EI,
                                 real_t mu,
                                 size_t nodes,
                                 beam::EulerBeamBCs bcs,
                                 real_t r_penalty)
    : EulerBeam(length, EI, mu, nodes, bcs)
    , elements(nodes - 1)
    , nodes(nodes)
    , ds(mesh.get_ds())
    , ndof_x(2 * nodes)
    , ndof_y(2 * nodes)
    , ndof_z(2 * nodes)
    , ndof_l(1 * nodes)
    , ndof(ndof_x + ndof_y + ndof_z)
    , offset_x(0)
    , offset_y(ndof_x)
    , offset_z(ndof_x + ndof_y)
    , r_penalty(r_penalty)
    , max_iter_outer(1000)
    , max_iter_inner(1000)
    , tol_outer(1e-5)
    , tol_inner(1e-5)
    , jacobian(Eigen::MatrixXd::Zero(ndof, ndof))
    , residual(Eigen::VectorXd::Zero(ndof))
    , u(Eigen::VectorXd::Zero(ndof))
    , lambda(Eigen::VectorXd::Zero(ndof_l))
  {
    apply_initial_condition(this->mesh);
  };

  /**
   *
   */
  real_t update_lambda(real_t omega = 1.0)
  {

    static constexpr std::array<real_t, 3> xi_q = { 0.1127016654,
                                                    0.5,
                                                    0.8872983346 };
    static constexpr std::array<real_t, 3> w_q = { 0.2777777778,
                                                   0.4444444444,
                                                   0.2777777778 };
    Eigen::Matrix<real_t, Eigen::Dynamic, 1> lambda_n =
      Eigen::Matrix<real_t, Eigen::Dynamic, 1>::Zero(ndof_l);

    // #pragma omp for
    for (size_t e = 0; e < elements; ++e) {

      size_t elem_nodes[] = { e, e + 1 };
      size_t idx_x[] = { offset_x + 2 * elem_nodes[0],
                         offset_x + 2 * elem_nodes[0] + 1,
                         offset_x + 2 * elem_nodes[1],
                         offset_x + 2 * elem_nodes[1] + 1 };
      size_t idx_y[] = { offset_y + 2 * elem_nodes[0],
                         offset_y + 2 * elem_nodes[0] + 1,
                         offset_y + 2 * elem_nodes[1],
                         offset_y + 2 * elem_nodes[1] + 1 };
      size_t idx_z[] = { offset_z + 2 * elem_nodes[0],
                         offset_z + 2 * elem_nodes[0] + 1,
                         offset_z + 2 * elem_nodes[1],
                         offset_z + 2 * elem_nodes[1] + 1 };
      size_t idx_l[] = { elem_nodes[0], elem_nodes[1] };

      real_t ux[] = {
        u[idx_x[0]], u[idx_x[1]], u[idx_x[2]], u[idx_x[3]]
      };
      real_t uy[] = {
        u[idx_y[0]], u[idx_y[1]], u[idx_y[2]], u[idx_y[3]]
      };
      real_t uz[] = {
        u[idx_z[0]], u[idx_z[1]], u[idx_z[2]], u[idx_z[3]]
      };
      real_t ul[] = { lambda[idx_l[0]], lambda[idx_l[1]] };

      real_t R_loc_l[2] = {0.,0.};

      for (size_t qi = 0; qi < 3; ++qi) {
        real_t xi = xi_q[qi];
        real_t w = w_q[qi];

        auto dH = CubicHermite<real_t>::derivs(xi, ds);
        auto M = LinearShape<real_t>::values(xi);

        real_t xp = 0;
        real_t yp = 0;
        real_t zp = 0;

        for (size_t i = 0; i < 4; i++) {
          xp += dH[i] * ux[i];
          yp += dH[i] * uy[i];
          zp += dH[i] * uz[i];
        }

        real_t l = 0;
        for (size_t i = 0; i < 2; i++) {
          l += M[i] * ul[i];
        }

        real_t S = xp * xp + yp * yp + zp * zp - 1.0;

        for (size_t a = 0; a < 2; ++a) {
          // Variation w.r.t. lambda
          R_loc_l[a] += S * M[a] * w * ds;
        }
      }

      for (size_t i = 0; i < 2; ++i) {
        lambda_n[idx_l[i]] += R_loc_l[i];
      }
    }

    lambda = lambda_n;

    return lambda_n.norm();
  }

  /**
   *
   */
  void assemble_residual(std::array<real_t, 3> load)
  {
    residual = assemble_residual_template<real_t>(u, load);
  }


  /**
   *
   */
  void assemble_residual(std::vector<std::array<real_t, 3>> load)
  {
    residual = assemble_residual_template<real_t>(u, load);
  }

  /**
   *
   */
  virtual void assemble_system(std::array<real_t, 3> load)
  {
    using AD = Eigen::AutoDiffScalar<Eigen::VectorXd>;
    using ADVec = Eigen::Matrix<AD, Eigen::Dynamic, 1>;

    ADVec x_ad(ndof);
    for (int i = 0; i < int(ndof); ++i) {
      Eigen::VectorXd seed = Eigen::VectorXd::Zero(ndof);
      seed(i) = 1.0;
      x_ad(i) = AD(u(i), seed);
    }

    ADVec R_ad = assemble_residual_template<AD>(x_ad, load);

    residual.resize(ndof);
    jacobian.resize(ndof, ndof);
    for (int i = 0; i < int(ndof); ++i) {
      residual(i) = R_ad(i).value();
      jacobian.row(i) = R_ad(i).derivatives().transpose();
    }
  }

  /**
   *
   */
  virtual void assemble_system(std::vector<std::array<real_t, 3>> load)
  {
    using AD = Eigen::AutoDiffScalar<Eigen::VectorXd>;
    using ADVec = Eigen::Matrix<AD, Eigen::Dynamic, 1>;

    ADVec x_ad(ndof);
    for (int i = 0; i < int(ndof); ++i) {
      Eigen::VectorXd seed = Eigen::VectorXd::Zero(ndof);
      seed(i) = 1.0;
      x_ad(i) = AD(u(i), seed);
    }

    ADVec R_ad = assemble_residual_template<AD>(x_ad, load);

    residual.resize(ndof);
    jacobian.resize(ndof, ndof);
    for (int i = 0; i < int(ndof); ++i) {
      residual(i) = R_ad(i).value();
      jacobian.row(i) = R_ad(i).derivatives().transpose();
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
  void apply_boundary_conditions() {
    apply_boundary_conditions(jacobian, residual);
  }

  /**
   * @brief Apply boundary conditions to the residual and jacobian
   *
   * Modifies the residual and jacobian according to boundary conditions:
   * - free_bc: No constraints
   * - simple_bc: Position constraints only
   * - clamped_bc: Position and slope constraints
   */
  void apply_boundary_conditions(Eigen::MatrixXd &A, Eigen::VectorXd &R)
  {
    for (size_t bi = 0; bi < 2; ++bi) {
      EulerBeamBCEnd bcend = boundary_conditions.end[bi];
      size_t ni = 0;
      std::vector<size_t> idx(6);

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
      std::vector<real_t> vals(6);

      switch (bctype) {
        case free_bc:
          idx = {};
          vals = {};
          break;
        case simple_bc:
          idx = { offset_x + 2 * ni + 0,
                  offset_y + 2 * ni + 0,
                  offset_z + 2 * ni + 0 };
          vals = { bcvals.position[0], bcvals.position[1], bcvals.position[2] };
          break;
        case clamped_bc:
          idx = { offset_x + 2 * ni + 0, offset_x + 2 * ni + 1,
                  offset_y + 2 * ni + 0, offset_y + 2 * ni + 1,
                  offset_z + 2 * ni + 0, offset_z + 2 * ni + 1 };
          vals = { bcvals.position[0], bcvals.slope[0],    bcvals.position[1],
                   bcvals.slope[1],    bcvals.position[2], bcvals.slope[2] };
          break;
        case point_force_bc:
          idx = {
            offset_x + 2 * ni + 0,
            offset_y + 2 * ni + 0,
            offset_z + 2 * ni + 0,
          };
          switch (bcend) {
            case left:
              vals = { bcvals.force[0], bcvals.force[1], bcvals.force[2] };
              break;
            case right:
              vals = { bcvals.force[0], bcvals.force[1], bcvals.force[2] };
              break;
          };
          break;
        case point_torque_bc:
          idx = {
            offset_x + 2 * ni + 1,
            offset_y + 2 * ni + 1,
            offset_z + 2 * ni + 1,
          };
          switch (bcend) {
            case left:
              vals = { bcvals.torque[0], bcvals.torque[1], bcvals.torque[2] };
              break;
            case right:
              vals = { bcvals.torque[0], bcvals.torque[1], bcvals.torque[2] };
              break;
          };
          break;
        default:
          break;
      }

      switch (bctype) {
        case point_force_bc:
          for (size_t i = 0; i < vals.size(); ++i)
            residual(idx[i]) -= vals[i];
          break;
        case point_torque_bc:
          for (size_t i = 0; i < vals.size(); ++i)
            residual(idx[i]) -= vals[i];
          break;
        default:
          for (size_t i = 0; i < vals.size(); ++i) {
            jacobian.row(idx[i]).setZero();
            jacobian.col(idx[i]).setZero();
            residual(idx[i]) = u(idx[i]) - vals[i];
            jacobian(idx[i], idx[i]) = 1.;
          }
          break;
      }
    }
  }

  /**
   * @brief Update the mesh object with the current solution data
   */
  void update_mesh()
  {

    std::vector<std::array<real_t, 3>>& centerline = mesh.get_centerline();
    std::vector<std::array<real_t, 3>>& slope = mesh.get_slope();
    std::vector<real_t>& s = mesh.get_curvilinear_axis();

    for (size_t i = 0; i < nodes; ++i) {
      centerline[i][0] = u(offset_x + 2 * i + 0);
      centerline[i][1] = u(offset_y + 2 * i + 0);
      centerline[i][2] = u(offset_z + 2 * i + 0);
      slope[i][0] = u(offset_x + 2 * i + 1);
      slope[i][1] = u(offset_y + 2 * i + 1);
      slope[i][2] = u(offset_z + 2 * i + 1);
    }
  }

  /**
   * @brief Assemble the residual vector for the nonlinear system
   *
   * @tparam T Scalar type (real_t or autodiff)
   * @param u Current solution vector 
   * @param load Current load on the beam
   * @return Assembled residual vector containing:
   *         - Bending energy terms
   *         - External load contributions
   *         - Inextensibility constraints
   *
   * Uses Gauss quadrature with 3 points for numerical integration.
   */
  template<typename T>
  Eigen::Matrix<T, Eigen::Dynamic, 1> assemble_residual_template(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& u,
    const std::array<real_t, 3> load) const
  {
    static constexpr real_t xi_q[] = { 0.1127016654, 0.5, 0.8872983346 };
    static constexpr real_t w_q[] = { 0.2777777778,
                                      0.4444444444,
                                      0.2777777778 };

    Eigen::Matrix<T, Eigen::Dynamic, 1> residual =
      Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(ndof);

    real_t Hq[3][4], dHq[3][4], ddHq[3][4], Mq[3][2];
    for (int qi = 0; qi < 3; ++qi) {
      auto H = CubicHermite<real_t>::values(xi_q[qi], ds);
      auto dH = CubicHermite<real_t>::derivs(xi_q[qi], ds);
      auto ddH = CubicHermite<real_t>::second_derivs(xi_q[qi], ds);
      auto M = LinearShape<real_t>::values(xi_q[qi]);
      for (int i = 0; i < 4; ++i) {
        Hq[qi][i] = H[i];
        dHq[qi][i] = dH[i];
        ddHq[qi][i] = ddH[i];
      }
      Mq[qi][0] = M[0];
      Mq[qi][1] = M[1];
    }

    // #pragma omp parallel for
    // #pragma unroll
    for (size_t e = 0; e < elements; ++e) {
      const size_t elem_nodes[] = { e, e + 1 };

      const size_t idx_x[] = { offset_x + 2 * elem_nodes[0] + 0,
                               offset_x + 2 * elem_nodes[0] + 1,
                               offset_x + 2 * elem_nodes[1] + 0,
                               offset_x + 2 * elem_nodes[1] + 1 };
      const size_t idx_y[] = { offset_y + 2 * elem_nodes[0] + 0,
                               offset_y + 2 * elem_nodes[0] + 1,
                               offset_y + 2 * elem_nodes[1] + 0,
                               offset_y + 2 * elem_nodes[1] + 1 };
      const size_t idx_z[] = { offset_z + 2 * elem_nodes[0] + 0,
                               offset_z + 2 * elem_nodes[0] + 1,
                               offset_z + 2 * elem_nodes[1] + 0,
                               offset_z + 2 * elem_nodes[1] + 1 };
      const size_t idx_l[] = { elem_nodes[0], elem_nodes[1] };

      T ux[] = { u[idx_x[0]], u[idx_x[1]], u[idx_x[2]], u[idx_x[3]] };
      T uy[] = { u[idx_y[0]], u[idx_y[1]], u[idx_y[2]], u[idx_y[3]] };
      T uz[] = { u[idx_z[0]], u[idx_z[1]], u[idx_z[2]], u[idx_z[3]] };
      T ul[] = { lambda[idx_l[0]], lambda[idx_l[1]] };

      T R_loc_x[4] = { T(0), T(0), T(0), T(0) };
      T R_loc_y[4] = { T(0), T(0), T(0), T(0) };
      T R_loc_z[4] = { T(0), T(0), T(0), T(0) };

      for (size_t qi = 0; qi < 3; ++qi) {
        T xp = 0, xpp = 0;
        T yp = 0, ypp = 0;
        T zp = 0, zpp = 0;

        for (size_t i = 0; i < 4; i++) {
          xp += dHq[qi][i] * ux[i];
          yp += dHq[qi][i] * uy[i];
          zp += dHq[qi][i] * uz[i];
          xpp += ddHq[qi][i] * ux[i];
          ypp += ddHq[qi][i] * uy[i];
          zpp += ddHq[qi][i] * uz[i];
        }

        T l = 0;
        for (size_t i = 0; i < 2; i++) {
          l += Mq[qi][i] * ul[i];
        }

        T S = xp * xp + yp * yp + zp * zp - 1.0;

        for (size_t a = 0; a < 4; ++a) {
          // Bending Energy Contributions
          const T coeff_bending = EI * ddHq[qi][a] * w_q[qi] * ds; 
          R_loc_x[a] += xpp * coeff_bending;
          R_loc_y[a] += ypp * coeff_bending;
          R_loc_z[a] += zpp * coeff_bending;
          // External load contribution
          R_loc_x[a] -= load[0] * Hq[qi][a] * w_q[qi] * ds;
          R_loc_y[a] -= load[1] * Hq[qi][a] * w_q[qi] * ds;
          R_loc_z[a] -= load[2] * Hq[qi][a] * w_q[qi] * ds;
          // Constraint contributions
          const T coeff_constraint = 2 * (l + r_penalty * S) * dHq[qi][a] * w_q[qi] * ds;
          R_loc_x[a] += xp * coeff_constraint;
          R_loc_y[a] += yp * coeff_constraint;
          R_loc_z[a] += zp * coeff_constraint;
        }
      }

      // Scatter local residual contribution to global residual
      for (size_t i = 0; i < 4; ++i) {
        residual[idx_x[i]] += R_loc_x[i];
        residual[idx_y[i]] += R_loc_y[i];
        residual[idx_z[i]] += R_loc_z[i];
      }
    }
    return residual;
  }


  /**
   * @brief Assemble the residual vector for the nonlinear system
   *
   * @tparam T Scalar type (real_t or autodiff)
   * @param u Current solution vector 
   * @param load Current load on the beam
   * @return Assembled residual vector containing:
   *         - Bending energy terms
   *         - External load contributions
   *         - Inextensibility constraints
   *
   * Uses Gauss quadrature with 3 points for numerical integration.
   */
  template<typename T>
  Eigen::Matrix<T, Eigen::Dynamic, 1> assemble_residual_template(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& u,
    const std::vector<std::array<real_t, 3>> load) const
  {

    BEAM_ASSERT(load.size() == nodes, "Nodes does not match load vector size.\n");

    static constexpr real_t xi_q[] = { 0.1127016654, 0.5, 0.8872983346 };
    static constexpr real_t w_q[] = { 0.2777777778,
                                      0.4444444444,
                                      0.2777777778 };

    Eigen::Matrix<T, Eigen::Dynamic, 1> residual =
      Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(ndof);

    real_t Hq[3][4], dHq[3][4], ddHq[3][4], Mq[3][2];
    for (int qi = 0; qi < 3; ++qi) {
      auto H = CubicHermite<real_t>::values(xi_q[qi], ds);
      auto dH = CubicHermite<real_t>::derivs(xi_q[qi], ds);
      auto ddH = CubicHermite<real_t>::second_derivs(xi_q[qi], ds);
      auto M = LinearShape<real_t>::values(xi_q[qi]);
      for (int i = 0; i < 4; ++i) {
        Hq[qi][i] = H[i];
        dHq[qi][i] = dH[i];
        ddHq[qi][i] = ddH[i];
      }
      Mq[qi][0] = M[0];
      Mq[qi][1] = M[1];
    }

    for (size_t e = 0; e < elements; ++e) {
      const size_t elem_nodes[] = { e, e + 1 };

      const size_t idx_x[] = { offset_x + 2 * elem_nodes[0] + 0,
                               offset_x + 2 * elem_nodes[0] + 1,
                               offset_x + 2 * elem_nodes[1] + 0,
                               offset_x + 2 * elem_nodes[1] + 1 };
      const size_t idx_y[] = { offset_y + 2 * elem_nodes[0] + 0,
                               offset_y + 2 * elem_nodes[0] + 1,
                               offset_y + 2 * elem_nodes[1] + 0,
                               offset_y + 2 * elem_nodes[1] + 1 };
      const size_t idx_z[] = { offset_z + 2 * elem_nodes[0] + 0,
                               offset_z + 2 * elem_nodes[0] + 1,
                               offset_z + 2 * elem_nodes[1] + 0,
                               offset_z + 2 * elem_nodes[1] + 1 };
      
      const size_t idx_load[] = { elem_nodes[0], elem_nodes[1] };
      
      const size_t idx_lam[] = { elem_nodes[0], elem_nodes[1] };

      T ux[] = { u[idx_x[0]], u[idx_x[1]], u[idx_x[2]], u[idx_x[3]] };
      T uy[] = { u[idx_y[0]], u[idx_y[1]], u[idx_y[2]], u[idx_y[3]] };
      T uz[] = { u[idx_z[0]], u[idx_z[1]], u[idx_z[2]], u[idx_z[3]] };

      real_t fx [] = { load[idx_load[0]][0], load[idx_load[1]][0] }; 
      real_t fy [] = { load[idx_load[0]][1], load[idx_load[1]][1] }; 
      real_t fz [] = { load[idx_load[0]][2], load[idx_load[1]][2] }; 

      T ul[] = { lambda[idx_lam[0]], lambda[idx_lam[1]] };

      T R_loc_x[4] = { T(0), T(0), T(0), T(0) };
      T R_loc_y[4] = { T(0), T(0), T(0), T(0) };
      T R_loc_z[4] = { T(0), T(0), T(0), T(0) };

      for (size_t qi = 0; qi < 3; ++qi) {
        T xp = 0, xpp = 0;
        T yp = 0, ypp = 0;
        T zp = 0, zpp = 0;


        for (size_t i = 0; i < 4; i++) {
          xp += dHq[qi][i] * ux[i];
          yp += dHq[qi][i] * uy[i];
          zp += dHq[qi][i] * uz[i];
          xpp += ddHq[qi][i] * ux[i];
          ypp += ddHq[qi][i] * uy[i];
          zpp += ddHq[qi][i] * uz[i];
        }

        T l = 0;
        real_t fxp = 0;
        real_t fyp = 0;
        real_t fzp = 0;

        for (size_t i = 0; i < 2; i++) {
          l   += Mq[qi][i] * ul[i];
          fxp += Mq[qi][i] * fx[i];
          fyp += Mq[qi][i] * fy[i];
          fzp += Mq[qi][i] * fz[i];
        }

        T S = xp * xp + yp * yp + zp * zp - 1.0;

        for (size_t a = 0; a < 4; ++a) {
          // Bending Energy Contributions
          const T coeff_bending = EI * ddHq[qi][a] * w_q[qi] * ds; 
          R_loc_x[a] += xpp * coeff_bending;
          R_loc_y[a] += ypp * coeff_bending;
          R_loc_z[a] += zpp * coeff_bending;
 
          // External load contribution
          const T coeff_load = Hq[qi][a] * w_q[qi] * ds;
          R_loc_x[a] -= fxp * Hq[qi][a] * w_q[qi] * ds;
          R_loc_y[a] -= fyp * Hq[qi][a] * w_q[qi] * ds;
          R_loc_z[a] -= fzp * Hq[qi][a] * w_q[qi] * ds;
 
          // Constraint contributions
          const T coeff_constraint = 2 * (l + r_penalty * S) * dHq[qi][a] * w_q[qi] * ds;
          R_loc_x[a] += xp * coeff_constraint;
          R_loc_y[a] += yp * coeff_constraint;
          R_loc_z[a] += zp * coeff_constraint;
        }
      }

      // Scatter local residual contribution to global residual
      for (size_t i = 0; i < 4; ++i) {
        residual[idx_x[i]] += R_loc_x[i];
        residual[idx_y[i]] += R_loc_y[i];
        residual[idx_z[i]] += R_loc_z[i];
      }
    }
    return residual;
  }
};

} // namespace beam