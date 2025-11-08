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
    , dof((2 * nodes * 3) + nodes)
    , r_penalty(r_penalty)
    , jacobian(Eigen::MatrixXd::Zero(dof, dof))
    , residual(Eigen::VectorXd::Zero(dof))
    , u(Eigen::VectorXd::Zero(dof))
  {
    apply_initial_condition(this->mesh);
  };

  ~EulerBeamStaticInextensibleMoM() {};


  /**
   * 
   */
  virtual void solve() override {
    solve({0.,0.,0.});
  }

  /**
   * 
   */
  virtual void solve(std::array<real_t, 3> load) override
  {
    real_t tol = 1e-6;
    size_t max_iter = 100;

    Eigen::LLT<Eigen::MatrixXd> cholesky;
    Eigen::LDLT<Eigen::MatrixXd> ldlt;

    for (size_t it = 0; it < max_iter; it++) {
      assemble_system(load);
      apply_boundary_conditions();

      real_t res_norm = residual.norm();
      std::cout << "iter " << it << "  ‖R‖ = " << res_norm << "\n";
      if (res_norm < tol) {
        std::cout << "Converged in " << it << " iters.\n";
        break;
      }

      Eigen::VectorXd delta_u;
      ldlt.compute(jacobian);
      if (1 == 1) { // ldlt.info() == Eigen::Success) {
        delta_u = jacobian.colPivHouseholderQr().solve(-residual);
      } else {
        delta_u = ldlt.solve(-residual);
      }

      // 4) update
      u += delta_u;

      update_mesh();
    }
  }

  virtual void apply_initial_condition()
  {
    size_t nodes = mesh.get_nodes();
    real_t ds = mesh.get_ds();

    size_t ndof_x = 2 * nodes;
    size_t ndof_y = 2 * nodes;
    size_t ndof_z = 2 * nodes;
    size_t ndof_l = nodes;

    size_t offset_x = 0;
    size_t offset_y = ndof_x;
    size_t offset_z = ndof_x + ndof_y;
    size_t offset_l = ndof_x + ndof_y + ndof_z;

    for (size_t i = 0; i < nodes; i++) {
      u(offset_x + 2 * i + 0) = ds * i;
      u(offset_x + 2 * i + 1) = 1.;
      u(offset_y + 2 * i + 0) = 0.;
      u(offset_y + 2 * i + 1) = 0.;
      u(offset_z + 2 * i + 0) = 0.;
      u(offset_z + 2 * i + 1) = 0.;
      u(offset_l + i) = 0.;
    }
  }

  virtual void apply_initial_condition(EulerBeamMesh& bmesh) override
  {

    size_t nodes = bmesh.get_nodes();

    if (bmesh.get_nodes() != mesh.get_nodes()) {
      BEAM_ABORT("");
    }

    auto centerline = bmesh.get_centerline();
    auto slope = bmesh.get_slope();

    size_t ndof_x = 2 * nodes;
    size_t ndof_y = 2 * nodes;
    size_t ndof_z = 2 * nodes;
    size_t ndof_l = nodes;

    size_t offset_x = 0;
    size_t offset_y = ndof_x;
    size_t offset_z = ndof_x + ndof_y;
    size_t offset_l = ndof_x + ndof_y + ndof_z;

    for (size_t i = 0; i < nodes; i++) {
      u(offset_x + 2 * i + 0) = centerline[i][0];
      u(offset_x + 2 * i + 1) = slope[i][0];
      u(offset_y + 2 * i + 0) = centerline[i][1];
      u(offset_y + 2 * i + 1) = slope[i][1];
      u(offset_z + 2 * i + 0) = centerline[i][2];
      u(offset_z + 2 * i + 1) = slope[i][2];
      u(offset_l + i) = 0.;
    }
  }

protected:
  size_t dimension;
  size_t elements;
  size_t dof;
  real_t r_penalty;

  Eigen::VectorXd residual;
  Eigen::MatrixXd jacobian;
  Eigen::VectorXd u;

  EulerBeamStaticInextensibleMoM(real_t length,
                                 real_t EI,
                                 real_t mu,
                                 size_t nodes,
                                 beam::EulerBeamBCs bcs,
                                 real_t r_penalty)
    : EulerBeam(length, EI, mu, nodes, bcs)
    , elements(nodes - 1)
    , dof((2 * nodes * 3) + nodes)
    , r_penalty(r_penalty)
    , jacobian(Eigen::MatrixXd::Zero(dof, dof))
    , residual(Eigen::VectorXd::Zero(dof))
    , u(Eigen::VectorXd::Zero(dof))
  {
    apply_initial_condition(this->mesh);
  };

  /**
   * 
   */
  void assemble_residual(std::array<real_t, 3> load) { residual = assemble_residual_impl<real_t>(u, load); }

  /**
   * 
   */
  virtual void assemble_system(std::array<real_t, 3> load)
  {
    using AD = Eigen::AutoDiffScalar<Eigen::VectorXd>;
    using ADVec = Eigen::Matrix<AD, Eigen::Dynamic, 1>;

    ADVec x_ad(dof);
    for (int i = 0; i < int(dof); ++i) {
      Eigen::VectorXd seed = Eigen::VectorXd::Zero(dof);
      seed(i) = 1.0;
      x_ad(i) = AD(u(i), seed);
    }

    ADVec R_ad = assemble_residual_impl<AD>(x_ad, load);

    residual.resize(dof);
    jacobian.resize(dof, dof);
    for (int i = 0; i < int(dof); ++i) {
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
  void apply_boundary_conditions()
  {
    size_t nodes = mesh.get_nodes();
    real_t h = mesh.get_ds();

    size_t ndof_x = 2 * nodes;
    size_t ndof_y = 2 * nodes;
    size_t ndof_z = 2 * nodes;
    size_t ndof_l = nodes;

    size_t offset_x = 0;
    size_t offset_y = ndof_x;
    size_t offset_z = ndof_x + ndof_y;
    size_t offset_l = ndof_x + ndof_y + ndof_z;

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
    size_t nodes = this->mesh.get_nodes();

    size_t ndof_x = 2 * nodes;
    size_t ndof_y = 2 * nodes;
    size_t ndof_z = 2 * nodes;
    size_t ndof_l = nodes;

    size_t offset_x = 0;
    size_t offset_y = ndof_x;
    size_t offset_z = ndof_x + ndof_y;
    size_t offset_l = ndof_x + ndof_y + ndof_z;

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
   * @return Assembled residual vector containing:
   *         - Bending energy terms
   *         - External load contributions
   *         - Inextensibility constraints
   *
   * Uses Gauss quadrature with 3 points for numerical integration.
   */
  template<typename T>
  Eigen::Matrix<T, Eigen::Dynamic, 1> assemble_residual_impl(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& u, const std::array<real_t,3> load) const
  {
    size_t nodes = mesh.get_nodes();
    real_t h = mesh.get_ds();

    real_t xi_q[] = { 0.1127016654, 0.5, 0.8872983346 };
    real_t w_q[] = { 0.2777777778, 0.4444444444, 0.2777777778 };

    size_t ndof_x = 2 * nodes;
    size_t ndof_y = 2 * nodes;
    size_t ndof_z = 2 * nodes;
    size_t ndof_l = nodes;

    size_t offset_x = 0;
    size_t offset_y = ndof_x;
    size_t offset_z = ndof_x + ndof_y;
    size_t offset_l = ndof_x + ndof_y + ndof_z;

    Eigen::Matrix<T, Eigen::Dynamic, 1> residual =
      Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(dof);

    for (size_t e = 0; e < elements; ++e) {
      std::vector<size_t> elem_nodes = { e, e + 1 };
      std::vector<size_t> idx_x = { offset_x + 2 * elem_nodes[0] + 0,
                                    offset_x + 2 * elem_nodes[0] + 1,
                                    offset_x + 2 * elem_nodes[1] + 0,
                                    offset_x + 2 * elem_nodes[1] + 1 };
      std::vector<size_t> idx_y = { offset_y + 2 * elem_nodes[0] + 0,
                                    offset_y + 2 * elem_nodes[0] + 1,
                                    offset_y + 2 * elem_nodes[1] + 0,
                                    offset_y + 2 * elem_nodes[1] + 1 };
      std::vector<size_t> idx_z = { offset_z + 2 * elem_nodes[0] + 0,
                                    offset_z + 2 * elem_nodes[0] + 1,
                                    offset_z + 2 * elem_nodes[1] + 0,
                                    offset_z + 2 * elem_nodes[1] + 1 };
      std::vector<size_t> idx_l = { offset_l + 1 * elem_nodes[0] + 0,
                                    offset_l + 1 * elem_nodes[1] + 0 };

      std::array<T, 4> ux = {
        u[idx_x[0]], u[idx_x[1]], u[idx_x[2]], u[idx_x[3]]
      };
      std::array<T, 4> uy = {
        u[idx_y[0]], u[idx_y[1]], u[idx_y[2]], u[idx_y[3]]
      };
      std::array<T, 4> uz = {
        u[idx_z[0]], u[idx_z[1]], u[idx_z[2]], u[idx_z[3]]
      };
      std::array<T, 2> ul = { u[idx_l[0]], u[idx_l[1]] };

      std::vector<T> R_loc_x(4, 0);
      std::vector<T> R_loc_y(4, 0);
      std::vector<T> R_loc_z(4, 0);
      std::vector<T> R_loc_l(2, 0);

      for (size_t qi = 0; qi < 3; ++qi) {
        real_t xi = xi_q[qi];
        real_t w = w_q[qi];

        auto H = CubicHermite<real_t>::values(xi, h);
        auto dH = CubicHermite<real_t>::derivs(xi, h);
        auto ddH = CubicHermite<real_t>::second_derivs(xi, h);
        auto M = LinearShape<real_t>::values(xi);

        T x = 0, xp = 0, xpp = 0;
        T y = 0, yp = 0, ypp = 0;
        T z = 0, zp = 0, zpp = 0;
        for (size_t i = 0; i < 4; i++) {
          x += H[i] * ux[i];
          xp += dH[i] * ux[i];
          xpp += ddH[i] * ux[i];
          y += H[i] * uy[i];
          yp += dH[i] * uy[i];
          ypp += ddH[i] * uy[i];
          z += H[i] * uz[i];
          zp += dH[i] * uz[i];
          zpp += ddH[i] * uz[i];
        }

        // T l = beam::dot<T,2>(M, ul);
        T l = 0;
        for (size_t i = 0; i < 2; i++) {
          l += M[i] * ul[i];
        }

        T S = xp * xp + yp * yp + zp * zp - 1.0;

        for (size_t a = 0; a < 4; ++a) {
          // Bending Energy Contributions
          R_loc_x[a] += EI * xpp * ddH[a] * w * h;
          R_loc_y[a] += EI * ypp * ddH[a] * w * h;
          R_loc_z[a] += EI * zpp * ddH[a] * w * h;
          // External load contribution
          R_loc_x[a] -= load[0] * H[a] * w * h;
          R_loc_y[a] -= load[1] * H[a] * w * h;
          R_loc_z[a] -= load[2] * H[a] * w * h;
          // Constraint contributions
          R_loc_x[a] += 2 * (l + r_penalty * S) * xp * dH[a] * w * h;
          R_loc_y[a] += 2 * (l + r_penalty * S) * yp * dH[a] * w * h;
          R_loc_z[a] += 2 * (l + r_penalty * S) * zp * dH[a] * w * h;
        }

        for (size_t a = 0; a < 2; ++a) {
          // Variation w.r.t. lambda
          R_loc_l[a] += S * M[a] * w * h;
        }
      }

      // Scatter local residual contribution to global residual
      for (size_t i = 0; i < 4; ++i) {
        residual[idx_x[i]] += R_loc_x[i];
        residual[idx_y[i]] += R_loc_y[i];
        residual[idx_z[i]] += R_loc_z[i];
      }

      for (size_t i = 0; i < 2; ++i) {
        residual[idx_l[i]] += R_loc_l[i];
      }
    }
    return residual;
  }
};

} // namespace beam