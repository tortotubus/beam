#pragma once

#include "elff/models/rope/Rope.hpp"
#include "elff/fem/Shapes.hpp"

#include <array>
#include <vector>

#include <Eigen/Dense>

using namespace Eigen;

namespace ELFF {
namespace Models {

/**
 * @brief Static inextensible Euler beam solved with a method-of-multipliers
 * formulation.
 */
class RopeStatic : public Rope
{
public:
  /**
   * @brief Constructs a static inextensible beam model.
   *
   * @param length Beam length
   * @param nodes Number of discretization nodes
   * @param bcs Boundary conditions at the beam ends
   * @param r_penalty Penalty parameter used in the inextensibility constraint
   */
  RopeStatic(real_t length, size_t nodes, Rope::RopeBCs bcs, real_t r_penalty);

  /**
   * @brief Destroys the beam model and releases any owned resources.
   */
  ~RopeStatic();

  /**
   * @brief Solves the static beam problem with no external load.
   */
  virtual void solve() override;

  /**
   * @brief Solves the static beam problem under uniform loading.
   *
   * @param load Uniform distributed load vector
   */
  void solve(std::array<real_t, 3> load) override;

  /**
   * @brief Solves the static beam problem under nodal loading.
   *
   * @param load Load vector specified at the mesh nodes
   */
  void solve(std::vector<std::array<real_t, 3>> load);

  /**
   * @brief Applies the model's default initial condition.
   */
  virtual void apply_initial_condition();

  /**
   * @brief Initializes the solver state from a supplied beam mesh.
   *
   * @param bmesh Beam mesh containing the initial geometry
   */
  virtual void apply_initial_condition(RopeMesh& bmesh) override;

protected:
  /**
   * @brief Spatial dimension of the beam formulation.
   */
  size_t dimension;
  /**
   * @brief Element and node counts used by the discretization.
   */
  size_t elements, nodes;
  /**
   * @brief Uniform spacing between nodes along the beam centerline.
   */
  real_t ds;
  /**
   * @brief Degrees of freedom in each Cartesian block and in the constraint
   * block.
   */
  size_t ndof_x, ndof_y, ndof_z, ndof_l;
  /**
   * @brief Offsets into the global state vector for each degree-of-freedom
   * block.
   */
  size_t offset_x, offset_y, offset_z, offset_l;
  /**
   * @brief Total number of algebraic degrees of freedom.
   */
  size_t ndof;
  /**
   * @brief Penalty parameter used to enforce inextensibility.
   */
  real_t r_penalty;
  /**
   * @brief Outer and inner nonlinear-iteration limits.
   */
  size_t max_iter_outer, max_iter_inner;
  /**
   * @brief Outer and inner nonlinear solver tolerances.
   */
  real_t tol_outer, tol_inner;

  /**
   * @brief Current residual vector.
   */
  VectorXd residual;
  /**
   * @brief Current Jacobian and mass-matrix storage.
   */
  MatrixXd jacobian, mass;
  /**
   * @brief Global state vector and Lagrange multiplier vector.
   */
  VectorXd u, lambda;

  /**
   * @brief Constructs the shared base state for dynamic derived classes.
   *
   * @param length Beam length 
   * @param mu Mass per unit length
   * @param nodes Number of discretization nodes
   * @param bcs Boundary conditions at the beam ends
   * @param r_penalty Penalty parameter used in the inextensibility constraint
   */
  RopeStatic(real_t length, 
             real_t mu,
             size_t nodes,
             EulerBeam::EulerBeamBCs bcs,
             real_t r_penalty);

  /**
   * @brief Updates the Lagrange multiplier iterate.
   *
   * @param omega Relaxation parameter for the multiplier update
   * @return Norm of the multiplier correction
   */
  real_t update_lambda(real_t omega = 1.0);

  /**
   * @brief Assembles the residual for uniform loading.
   *
   * @param load Uniform distributed load vector
   */
  void assemble_residual(std::array<real_t, 3> load);

  /**
   * @brief Assembles the residual for nodal loading.
   *
   * @param load Load vector specified at the mesh nodes
   */
  void assemble_residual(std::vector<std::array<real_t, 3>> load);

  /**
   * @brief Assembles the nonlinear system for uniform loading.
   *
   * @param load Uniform distributed load vector
   */
  virtual void assemble_system(std::array<real_t, 3> load);

  /**
   * @brief Assembles the nonlinear system for nodal loading.
   *
   * @param load Load vector specified at the mesh nodes
   */
  virtual void assemble_system(std::vector<std::array<real_t, 3>> load);

  /**
   * @brief Applies the configured boundary conditions to the current system.
   */
  void apply_boundary_conditions();

  /**
   * @brief Applies the configured boundary conditions to a supplied matrix and
   * residual pair.
   *
   * @param A System matrix to modify
   * @param R Residual vector to modify
   */
  void apply_boundary_conditions(MatrixXd& A, VectorXd& R);

  /**
   * @brief Updates the beam mesh from the current solution state.
   */
  void update_mesh();

  template<typename T>
  /**
   * @brief Assembles the static residual for uniform loading.
   *
   * @tparam T Scalar type used for residual assembly
   * @param u State vector at which to evaluate the residual
   * @param load Uniform distributed load vector
   * @return Residual vector for the current state
   */
  Matrix<T, Dynamic, 1> assemble_residual_template(
    const Matrix<T, Dynamic, 1>& u,
    std::array<real_t, 3> load) const
  {
    static constexpr real_t xi_q[] = { 0.1127016654, 0.5, 0.8872983346 };
    static constexpr real_t w_q[] = { 0.2777777778,
                                      0.4444444444,
                                      0.2777777778 };

    Matrix<T, Dynamic, 1> residual = Matrix<T, Dynamic, 1>::Zero(ndof);

    real_t Hq[3][4], dHq[3][4], Mq[3][2];
    for (int qi = 0; qi < 3; ++qi) {
      auto H = ELFF::FEM::CubicHermite<real_t>::values(xi_q[qi], ds);
      auto dH = ELFF::FEM::CubicHermite<real_t>::derivs(xi_q[qi], ds);
      auto M = ELFF::FEM::LinearShape<real_t>::values(xi_q[qi]);
      for (int i = 0; i < 4; ++i) {
        Hq[qi][i] = H[i];
        dHq[qi][i] = dH[i];
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
      const size_t idx_l[] = { elem_nodes[0], elem_nodes[1] };

      T ux[] = { u[idx_x[0]], u[idx_x[1]], u[idx_x[2]], u[idx_x[3]] };
      T uy[] = { u[idx_y[0]], u[idx_y[1]], u[idx_y[2]], u[idx_y[3]] };
      T uz[] = { u[idx_z[0]], u[idx_z[1]], u[idx_z[2]], u[idx_z[3]] };
      T ul[] = { lambda[idx_l[0]], lambda[idx_l[1]] };

      T R_loc_x[4] = { T(0), T(0), T(0), T(0) };
      T R_loc_y[4] = { T(0), T(0), T(0), T(0) };
      T R_loc_z[4] = { T(0), T(0), T(0), T(0) };

      for (size_t qi = 0; qi < 3; ++qi) {
        T xp = 0;
        T yp = 0;
        T zp = 0;

        for (size_t i = 0; i < 4; i++) {
          xp += dHq[qi][i] * ux[i];
          yp += dHq[qi][i] * uy[i];
          zp += dHq[qi][i] * uz[i];
        }

        T l = 0;
        for (size_t i = 0; i < 2; i++) {
          l += Mq[qi][i] * ul[i];
        }

        T S = xp * xp + yp * yp + zp * zp - 1.0;

        for (size_t a = 0; a < 4; ++a) {
          R_loc_x[a] -= load[0] * Hq[qi][a] * w_q[qi] * ds;
          R_loc_y[a] -= load[1] * Hq[qi][a] * w_q[qi] * ds;
          R_loc_z[a] -= load[2] * Hq[qi][a] * w_q[qi] * ds;

          const T coeff_constraint =
            2 * (l + r_penalty * S) * dHq[qi][a] * w_q[qi] * ds;
          R_loc_x[a] += xp * coeff_constraint;
          R_loc_y[a] += yp * coeff_constraint;
          R_loc_z[a] += zp * coeff_constraint;
        }
      }

      for (size_t i = 0; i < 4; ++i) {
        residual[idx_x[i]] += R_loc_x[i];
        residual[idx_y[i]] += R_loc_y[i];
        residual[idx_z[i]] += R_loc_z[i];
      }
    }

    return residual;
  }

  template<typename T>
  /**
   * @brief Assembles the static residual for nodal loading.
   *
   * @tparam T Scalar type used for residual assembly
   * @param u State vector at which to evaluate the residual
   * @param load Load vector specified at the mesh nodes
   * @return Residual vector for the current state
   */
  Matrix<T, Dynamic, 1> assemble_residual_template(
    const Matrix<T, Dynamic, 1>& u,
    const std::vector<std::array<real_t, 3>> load) const
  {
    ELFF_ASSERT(load.size() == nodes,
                "Nodes does not match load vector size.\n");

    static constexpr real_t xi_q[] = { 0.1127016654, 0.5, 0.8872983346 };
    static constexpr real_t w_q[] = { 0.2777777778,
                                      0.4444444444,
                                      0.2777777778 };

    Matrix<T, Dynamic, 1> residual = Matrix<T, Dynamic, 1>::Zero(ndof);

    real_t Hq[3][4];
    real_t dHq[3][4];
    real_t Mq[3][2];
    for (int qi = 0; qi < 3; ++qi) {
      auto H = ELFF::FEM::CubicHermite<real_t>::values(xi_q[qi], ds);
      auto dH = ELFF::FEM::CubicHermite<real_t>::derivs(xi_q[qi], ds);
      auto M = ELFF::FEM::LinearShape<real_t>::values(xi_q[qi]);
      for (int i = 0; i < 4; ++i) {
        Hq[qi][i] = H[i];
        dHq[qi][i] = dH[i];
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

      real_t fx[] = { load[idx_load[0]][0], load[idx_load[1]][0] };
      real_t fy[] = { load[idx_load[0]][1], load[idx_load[1]][1] };
      real_t fz[] = { load[idx_load[0]][2], load[idx_load[1]][2] };

      T ul[] = { lambda[idx_lam[0]], lambda[idx_lam[1]] };

      T R_loc_x[4] = { T(0), T(0), T(0), T(0) };
      T R_loc_y[4] = { T(0), T(0), T(0), T(0) };
      T R_loc_z[4] = { T(0), T(0), T(0), T(0) };

      for (size_t qi = 0; qi < 3; ++qi) {
        T xp = 0;
        T yp = 0;
        T zp = 0;

        for (size_t i = 0; i < 4; i++) {
          xp += dHq[qi][i] * ux[i];
          yp += dHq[qi][i] * uy[i];
          zp += dHq[qi][i] * uz[i];
        }

        T l = 0;
        real_t fxp = 0;
        real_t fyp = 0;
        real_t fzp = 0;

        for (size_t i = 0; i < 2; i++) {
          l += Mq[qi][i] * ul[i];
          fxp += Mq[qi][i] * fx[i];
          fyp += Mq[qi][i] * fy[i];
          fzp += Mq[qi][i] * fz[i];
        }

        T S = xp * xp + yp * yp + zp * zp - 1.0;

        for (size_t a = 0; a < 4; ++a) {
          R_loc_x[a] -= fxp * Hq[qi][a] * w_q[qi] * ds;
          R_loc_y[a] -= fyp * Hq[qi][a] * w_q[qi] * ds;u
          R_loc_z[a] -= fzp * Hq[qi][a] * w_q[qi] * ds;

          const T coeff_constraint =
            2 * (l + r_penalty * S) * dHq[qi][a] * w_q[qi] * ds;
          R_loc_x[a] += xp * coeff_constraint;
          R_loc_y[a] += yp * coeff_constraint;
          R_loc_z[a] += zp * coeff_constraint;
        }
      }

      for (size_t i = 0; i < 4; ++i) {
        residual[idx_x[i]] += R_loc_x[i];
        residual[idx_y[i]] += R_loc_y[i];
        residual[idx_z[i]] += R_loc_z[i];
      }
    }

    return residual;
  }
};

} // namespace Models
} // namespace ELFF
