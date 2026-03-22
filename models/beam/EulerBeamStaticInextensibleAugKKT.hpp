#pragma once

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

using namespace Eigen;

namespace ELFF {
namespace Models {

class EulerBeamStaticInextensibleAugKKT : public EulerBeam
{
public:
  EulerBeamStaticInextensibleAugKKT(real_t length,
                                    real_t EI,
                                    size_t nodes,
                                    EulerBeam::EulerBeamBCs bcs,
                                    real_t r_penalty);

  ~EulerBeamStaticInextensibleAugKKT();

  /**
   *
   */
  virtual void solve() override;

  /**
   *
   */
  virtual void solve(std::array<real_t, 3> load) override;

  virtual void apply_initial_condition();

  virtual void apply_initial_condition(EulerBeamMesh& bmesh) override;

protected:
  size_t dimension;
  size_t elements, nodes;
  real_t ds;
  size_t ndof_x, ndof_y, ndof_z, ndof_l;
  size_t offset_x, offset_y, offset_z, offset_l;
  size_t ndof;
  real_t r_penalty;
  size_t max_iter;
  real_t tol;

  VectorXd residual;
  MatrixXd jacobian;
  VectorXd u;

  EulerBeamStaticInextensibleAugKKT(real_t length,
                                    real_t EI,
                                    real_t mu,
                                    size_t nodes,
                                    EulerBeam::EulerBeamBCs bcs,
                                    real_t r_penalty);

  /**
   *
   */
  void assemble_residual(std::array<real_t, 3> load);

  /**
   *
   */
  virtual void assemble_system(std::array<real_t, 3> load);

  /**
   * @brief Apply boundary conditions to the residual and jacobian
   *
   * Modifies the residual and jacobian according to boundary conditions:
   * - free_bc: No constraints
   * - simple_bc: Position constraints only
   * - clamped_bc: Position and slope constraints
   */
  void apply_boundary_conditions();

  /**
   * @brief Update the mesh object with the current solution data
   */
  void update_mesh();

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
  Matrix<T, Dynamic, 1> assemble_residual_template(
    const Matrix<T, Dynamic, 1>& u,
    const std::array<real_t, 3> load) const
  {
    real_t xi_q[] = { 0.1127016654, 0.5, 0.8872983346 };
    real_t w_q[] = { 0.2777777778, 0.4444444444, 0.2777777778 };

    Matrix<T, Dynamic, 1> residual =
      Matrix<T, Dynamic, 1>::Zero(ndof);

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

        auto H = CubicHermite<real_t>::values(xi, ds);
        auto dH = CubicHermite<real_t>::derivs(xi, ds);
        auto ddH = CubicHermite<real_t>::second_derivs(xi, ds);
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

        // T l = ELFF::dot<T,2>(M, ul);
        T l = 0;
        for (size_t i = 0; i < 2; i++) {
          l += M[i] * ul[i];
        }

        T S = xp * xp + yp * yp + zp * zp - 1.0;

        for (size_t a = 0; a < 4; ++a) {
          // Bending Energy Contributions
          R_loc_x[a] += EI * xpp * ddH[a] * w * ds;
          R_loc_y[a] += EI * ypp * ddH[a] * w * ds;
          R_loc_z[a] += EI * zpp * ddH[a] * w * ds;
          // External load contribution
          R_loc_x[a] -= load[0] * H[a] * w * ds;
          R_loc_y[a] -= load[1] * H[a] * w * ds;
          R_loc_z[a] -= load[2] * H[a] * w * ds;
          // Constraint contributions
          R_loc_x[a] += 2 * (l + r_penalty * S) * xp * dH[a] * w * ds;
          R_loc_y[a] += 2 * (l + r_penalty * S) * yp * dH[a] * w * ds;
          R_loc_z[a] += 2 * (l + r_penalty * S) * zp * dH[a] * w * ds;
        }

        for (size_t a = 0; a < 2; ++a) {
          // Variation w.r.t. lambda
          R_loc_l[a] += S * M[a] * w * ds;
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

} // namespace Models
} // namespace ELFF
