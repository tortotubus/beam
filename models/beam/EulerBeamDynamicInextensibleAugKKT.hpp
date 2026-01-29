#pragma once

#include "EulerBeamStaticInextensibleAugKKT.hpp"

namespace ELFF {
class EulerBeamDynamicInextensibleAugKKT
  : public EulerBeamStaticInextensibleAugKKT
{
public:
  EulerBeamDynamicInextensibleAugKKT(real_t length,
                                     real_t EI,
                                     real_t mu,
                                     size_t nodes,
                                     EulerBeam::EulerBeamBCs bcs,
                                     real_t r_penalty)
    : EulerBeamStaticInextensibleAugKKT(length, EI, mu, nodes, bcs, r_penalty)
    , u_prev(Eigen::VectorXd::Zero(ndof))
    , u_prev_prev(Eigen::VectorXd::Zero(ndof))
    , v_prev(Eigen::VectorXd::Zero(ndof))
    , a_prev(Eigen::VectorXd::Zero(ndof)) {};

  virtual void solve(real_t dt, std::array<real_t, 3> load) override
  {
    // solve_be(dt, load);
    // const real_t alpha = -0.1;
    const real_t alpha = 0.0;
    const real_t gamma = 0.5 - alpha;
    const real_t beta = 0.25 * (1 - alpha) * (1 - alpha);
    solve_newmark(dt, load, beta, gamma);
  }

  void solve_newmark(real_t dt,
                     std::array<real_t, 3> load,
                     real_t beta,
                     real_t gamma)
  {
    if (time_iter == 0) {
      Eigen::VectorXd R0 =
        EulerBeamStaticInextensibleAugKKT ::assemble_residual_template<real_t>(
          u, load);

      for (size_t n = 0; n < nodes; ++n) {
        size_t ix = offset_x + 2 * n;
        size_t iy = offset_y + 2 * n;
        size_t iz = offset_z + 2 * n;
        a_prev(ix) = (-R0(ix)) / (mu*ds);
        a_prev(iy) = (-R0(iy)) / (mu*ds);
        a_prev(iz) = (-R0(iz)) / (mu*ds);
      }
    }

    // real_t tol = 1e-9;
    // size_t max_iter = 500;

    u_prev = u;

    for (size_t it = 0; it < max_iter; it++) {
      assemble_system_newmark(dt, load, beta, gamma);
      apply_boundary_conditions();
      real_t res_norm = residual.norm();
      if (res_norm < tol) {
        break;
      }
      Eigen::VectorXd delta_u = jacobian.colPivHouseholderQr().solve(-residual);
      u += delta_u;
    }

    update_mesh();

    size_t nodes = mesh.get_nodes();

    for (size_t ni = 1; ni < nodes; ++ni) {
      size_t ix = offset_x + 2 * ni;
      size_t iy = offset_y + 2 * ni;
      size_t iz = offset_z + 2 * ni;

      auto upd = [&](size_t i) {
        real_t a_new = (u(i) - u_prev(i) - dt * v_prev(i)) / (beta * dt * dt) -
                       ((1.0 - 2.0 * beta) / (2.0 * beta)) * a_prev(i);
        real_t v_new =
          v_prev(i) + dt * ((1.0 - gamma) * a_prev(i) + gamma * a_new);
        a_prev(i) = a_new;
        v_prev(i) = v_new;
      };

      upd(ix);
      upd(iy);
      upd(iz);
    }

    u_prev_prev = u_prev; // carry old n-1
    u_prev = u;           // store new n

    time_iter++;
    t += dt;
  }

  void solve_be(real_t dt, std::array<real_t, 3> load)
  {
    real_t tol = 1e-6;
    size_t max_iter = 100;

    for (size_t it = 0; it < max_iter; it++) {
      assemble_system_be(dt, load);
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

    time_iter++;
    t += dt;
  }

  void apply_initial_condition() override
  {
    EulerBeamStaticInextensibleAugKKT::apply_initial_condition();
    u_prev = u;
    u_prev_prev = u;
  }

  void apply_initial_condition(EulerBeamMesh& bmesh) override
  {
    EulerBeamStaticInextensibleAugKKT::apply_initial_condition(bmesh);
    u_prev = u;
    u_prev_prev = u;
  }

protected:
  Eigen::VectorXd v_prev;
  Eigen::VectorXd a_prev;
  Eigen::VectorXd u_prev, u_prev_prev;
  std::array<real_t, 3> load_prev;

  /**
   *
   */
  void assemble_system_newmark(real_t dt,
                               std::array<real_t, 3> load,
                               real_t beta,
                               real_t gamma)
  {
    using AD = Eigen::AutoDiffScalar<Eigen::VectorXd>;
    using ADVec = Eigen::Matrix<AD, Eigen::Dynamic, 1>;

    ADVec x_ad(ndof);

    for (int i = 0; i < int(ndof); ++i) {
      Eigen::VectorXd seed = Eigen::VectorXd::Zero(ndof);
      seed(i) = 1.0;
      x_ad(i) = AD(u(i), seed);
    }

    ADVec R_ad = assemble_residual_newmark<AD>(x_ad, dt, load, beta, gamma);

    residual.resize(ndof);
    jacobian.resize(ndof, ndof);
    for (int i = 0; i < int(ndof); ++i) {
      residual(i) = R_ad(i).value();
      jacobian.row(i) = R_ad(i).derivatives().transpose();
    }
  }

  template<typename T>
  Eigen::Matrix<T, Eigen::Dynamic, 1> assemble_residual_newmark(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& u,
    real_t dt,
    std::array<real_t, 3> load,
    real_t beta,
    real_t /*gamma*/) const
  {
    Eigen::Matrix<T, Eigen::Dynamic, 1> res =
      EulerBeamStaticInextensibleAugKKT::assemble_residual_template<T>(u, load);

    // Guards to avoid NaNs/Infs:
    if (!(dt > 0.0))
      throw std::runtime_error("Newmark: dt must be > 0");
    if (!(beta > 0.0))
      throw std::runtime_error("Newmark: beta must be > 0");

    // Keep these as doubles (scalars), not T:
    const double inv = 1.0 / (beta * dt * dt); // = 1/(β Δt²)
    const double inv_bt = 1.0 / (beta * dt);   // = 1/(β Δt)
    const double kappa = (1.0 - 2.0 * beta) / (2.0 * beta);

    auto newmark_a = [&](Eigen::Index i) -> T {
      // u(i) is AD (T); u_prev/v_prev/a_prev are doubles
      return inv * (u(i) - u_prev(i)) - inv_bt * v_prev(i) - kappa * a_prev(i);
    };

    for (size_t n = 1; n < nodes; ++n) {
      const Eigen::Index ix = static_cast<Eigen::Index>(offset_x + 2 * n);
      const Eigen::Index iy = static_cast<Eigen::Index>(offset_y + 2 * n);
      const Eigen::Index iz = static_cast<Eigen::Index>(offset_z + 2 * n);

      const T ax = newmark_a(ix);
      const T ay = newmark_a(iy);
      const T az = newmark_a(iz);

      res(ix) += mu * ds * ax; // mu is double; AD ⊗ scalar is fine
      res(iy) += mu * ds * ay;
      res(iz) += mu * ds * az;
    }

    return res;
  }

  /**
   *
   */
  void assemble_system_be(real_t dt, std::array<real_t, 3> load)
  {
    std::cout << "Dynamic::assemble_be_system()\n";
    using AD = Eigen::AutoDiffScalar<Eigen::VectorXd>;
    using ADVec = Eigen::Matrix<AD, Eigen::Dynamic, 1>;

    ADVec x_ad(ndof);

    for (int i = 0; i < int(ndof); ++i) {
      Eigen::VectorXd seed = Eigen::VectorXd::Zero(ndof);
      seed(i) = 1.0;
      x_ad(i) = AD(u(i), seed);
    }

    ADVec R_ad = assemble_residual_be<AD>(x_ad, dt, load);

    residual.resize(ndof);
    jacobian.resize(ndof, ndof);
    for (int i = 0; i < int(ndof); ++i) {
      residual(i) = R_ad(i).value();
      jacobian.row(i) = R_ad(i).derivatives().transpose();
    }
  }

  /**
   * @brief Assemble the residual vector for the nonlinear system
   *
   * @tparam T Scalar type (real_t or autodiff)
   * @param u Current solution vector
   * @param dt Size of time step
   * @param load Current load on the beam
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
      EulerBeamStaticInextensibleAugKKT::assemble_residual_template<T>(u, load);

    for (size_t n = 0; n < nodes; n++) {
      res(offset_x + 2 * n) +=
        ((mu*ds) / (dt * dt)) * (u[offset_x + 2 * n] - 2 * u_prev[offset_x + 2 * n] +
                            u_prev_prev[offset_x + 2 * n]);
      res(offset_y + 2 * n) +=
        ((mu*ds) / (dt * dt)) * (u[offset_y + 2 * n] - 2 * u_prev[offset_y + 2 * n] +
                            u_prev_prev[offset_y + 2 * n]);
      res(offset_z + 2 * n) +=
        ((mu*ds) / (dt * dt)) * (u[offset_z + 2 * n] - 2 * u_prev[offset_z + 2 * n] +
                            u_prev_prev[offset_z + 2 * n]);
    }

    return res;
  }
};
}