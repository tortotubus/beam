#pragma once

#include "EulerBeamStaticInextensibleMoM.hpp"

namespace beam {
class EulerBeamDynamicInextensibleMoM : public EulerBeamStaticInextensibleMoM
{
public:
  EulerBeamDynamicInextensibleMoM(real_t length,
                                  real_t EI,
                                  real_t mu,
                                  size_t nodes,
                                  beam::EulerBeamBCs bcs,
                                  real_t r_penalty)
    : EulerBeamStaticInextensibleMoM(length, EI, mu, nodes, bcs, r_penalty)
    , u_prev(Eigen::VectorXd::Zero(dof))
    , u_prev_prev(Eigen::VectorXd::Zero(dof))
    , v_prev(Eigen::VectorXd::Zero(dof))
    , a_prev(Eigen::VectorXd::Zero(dof)) {};

  virtual void solve(real_t dt, std::array<real_t, 3> load) override
  {
    solve_be(dt, load);
  }

  void solve_newmark(real_t dt, std::array<real_t, 3> load)
  {
    real_t tol = 1e-6;
    size_t max_iter = 100;

    for (size_t it = 0; it < max_iter; it++) {
      assemble_be_system(dt, load);
      apply_boundary_conditions();
      real_t res_norm = residual.norm();
      if (res_norm < tol) {
        break;
      }
      Eigen::VectorXd delta_u = jacobian.colPivHouseholderQr().solve(-residual);
      u += delta_u;
    }

    update_mesh();
    u_prev_prev = u_prev; // carry old n-1
    u_prev = u;           // store new n
  }

  void solve_be(real_t dt, std::array<real_t, 3> load)
  {
    real_t tol = 1e-6;
    size_t max_iter = 100;

    for (size_t it = 0; it < max_iter; it++) {
      assemble_be_system(dt, load);
      apply_boundary_conditions();
      real_t res_norm = residual.norm();
      if (res_norm < tol) {
        break;
      }
      Eigen::VectorXd delta_u = jacobian.colPivHouseholderQr().solve(-residual);
      u += delta_u;
    }

    update_mesh();
    u_prev_prev = u_prev; // carry old n-1
    u_prev = u;           // store new n
  }

  void apply_initial_condition() override
  {
    EulerBeamStaticInextensibleMoM::apply_initial_condition();
    u_prev = u;
    u_prev_prev = u;
  }

  void apply_initial_condition(EulerBeamMesh& bmesh) override
  {
    EulerBeamStaticInextensibleMoM::apply_initial_condition(bmesh);
  }

protected:
  Eigen::VectorXd v_prev;
  Eigen::VectorXd a_prev;
  Eigen::VectorXd u_prev, u_prev_prev;

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
  Eigen::Matrix<T, Eigen::Dynamic, 1> assemble_residual_newmark(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& u,
    real_t dt,
    std::array<real_t, 3> load,
    real_t beta,
    real_t gamma) const
  {
    Eigen::Matrix<T, Eigen::Dynamic, 1> res =
      EulerBeamStaticInextensibleMoM::assemble_residual_impl<T>(u, load);

    size_t nodes = mesh.get_nodes();

    size_t ndof_x = 2 * nodes;
    size_t ndof_y = 2 * nodes;
    size_t ndof_z = 2 * nodes;
    size_t ndof_l = nodes;

    size_t offset_x = 0;
    size_t offset_y = ndof_x;
    size_t offset_z = ndof_x + ndof_y;
    size_t offset_l = ndof_x + ndof_y + ndof_z;

    const T inv = T(1.0) / (T(beta) * T(dt * dt));
    const T kappa = (T(1.0) - T(2.0 * beta)) / (T(2.0 * beta));

    auto newmark_a = [&](size_t i) -> T {
      return inv * (u(i) - T(u_prev(i))) -
             (T(dt) / T(beta * dt * dt)) * T(v_prev(i)) - kappa * T(a_prev(i));
    };

    auto newmark_v = [&](const T& a_new, size_t i) -> T {
      return T(v_prev(i)) +
             T(dt) * ((T(1.0) - T(gamma)) * T(a_prev(i)) + T(gamma) + a_new);
    };

    for (size_t ni = 0; ni < nodes; ++ni) {
      size_t ix = offset_x + 2 * ni;
      size_t iy = offset_y + 2 * ni;
      size_t iz = offset_z + 2 * ni;

      T ax = newmark_a(ix);
      T ay = newmark_a(iy);
      T az = newmark_a(iz);

      res(ix) += T(mu) * ax;
      res(iy) += T(mu) * ay;
      res(iz) += T(mu) * az;
    }

    return res;
  }

  void assemble_be_system(real_t dt, std::array<real_t, 3> load)
  {
    std::cout << "Dynamic::assemble_be_system()\n";
    using AD = Eigen::AutoDiffScalar<Eigen::VectorXd>;
    using ADVec = Eigen::Matrix<AD, Eigen::Dynamic, 1>;

    ADVec x_ad(dof);

    for (int i = 0; i < int(dof); ++i) {
      Eigen::VectorXd seed = Eigen::VectorXd::Zero(dof);
      seed(i) = 1.0;
      x_ad(i) = AD(u(i), seed);
    }

    ADVec R_ad = assemble_residual_be<AD>(x_ad, dt, load);

    residual.resize(dof);
    jacobian.resize(dof, dof);
    for (int i = 0; i < int(dof); ++i) {
      residual(i) = R_ad(i).value();
      jacobian.row(i) = R_ad(i).derivatives().transpose();
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
  Eigen::Matrix<T, Eigen::Dynamic, 1> assemble_residual_be(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& u,
    real_t dt,
    std::array<real_t, 3> load) const
  {
    Eigen::Matrix<T, Eigen::Dynamic, 1> res =
      EulerBeamStaticInextensibleMoM::assemble_residual_impl<T>(u, load);

    size_t nodes = mesh.get_nodes();

    size_t ndof_x = 2 * nodes;
    size_t ndof_y = 2 * nodes;
    size_t ndof_z = 2 * nodes;
    size_t ndof_l = nodes;

    size_t offset_x = 0;
    size_t offset_y = ndof_x;
    size_t offset_z = ndof_x + ndof_y;
    size_t offset_l = ndof_x + ndof_y + ndof_z;

    for (size_t n = 0; n < nodes; n++) {
      res(offset_x + 2 * n) +=
        (mu / (dt * dt)) * (u[offset_x + 2 * n] - 2 * u_prev[offset_x + 2 * n] +
                            u_prev_prev[offset_x + 2 * n]);
      res(offset_y + 2 * n) +=
        (mu / (dt * dt)) * (u[offset_y + 2 * n] - 2 * u_prev[offset_y + 2 * n] +
                            u_prev_prev[offset_y + 2 * n]);
      res(offset_z + 2 * n) +=
        (mu / (dt * dt)) * (u[offset_z + 2 * n] - 2 * u_prev[offset_z + 2 * n] +
                            u_prev_prev[offset_z + 2 * n]);
    }

    return res;
  }
};
}