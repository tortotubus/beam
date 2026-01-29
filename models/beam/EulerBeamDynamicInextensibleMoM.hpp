#pragma once

#include "EulerBeamStaticInextensibleMoM.hpp"

namespace ELFF {
class EulerBeamDynamicInextensibleMoM : public EulerBeamStaticInextensibleMoM
{
public:
  EulerBeamDynamicInextensibleMoM(real_t length,
                                  real_t EI,
                                  real_t mu,
                                  size_t nodes,
                                  EulerBeam::EulerBeamBCs bcs,
                                  real_t r_penalty)
    : EulerBeamStaticInextensibleMoM(length, EI, mu, nodes, bcs, r_penalty)
    , u_prev(Eigen::VectorXd::Zero(ndof))
    , u_prev_prev(Eigen::VectorXd::Zero(ndof))
    , v_prev(Eigen::VectorXd::Zero(ndof))
    , a_prev(Eigen::VectorXd::Zero(ndof))
    , mass(Eigen::MatrixXd::Zero(ndof, ndof)) {};

  virtual void solve(real_t dt, std::array<real_t, 3> load) override
  {
    const real_t alpha = 0;
    const real_t gamma = 0.5 - alpha;
    const real_t beta = 0.25 * (1 - alpha) * (1 - alpha);
    solve_newmark(dt, load, beta, gamma);
  }

  virtual void solve(real_t dt, std::vector<std::array<real_t, 3>> load)
  {
    ELFF_ASSERT(load.size() == nodes,
                "Size of load vector must equal number of nodes.");
    const real_t alpha = 0;
    const real_t gamma = 0.5 - alpha;
    const real_t beta = 0.25 * (1 - alpha) * (1 - alpha);
    solve_newmark(dt, load, beta, gamma);
  }

  void solve_newmark(real_t dt,
                     std::vector<std::array<real_t, 3>> load,
                     real_t beta,
                     real_t gamma)
  {

    // Eigen::LDLT<Eigen::MatrixXd> solver;
    // solver.setTolerance(tol_inner);
    Eigen::ConjugateGradient<Eigen::MatrixXd, Eigen::Upper | Eigen::Lower>
      solver;

    if (time_iter == 0) {
      Eigen::VectorXd R0 =
        EulerBeamStaticInextensibleMoM ::assemble_residual_template<real_t>(
          u, load);
      apply_boundary_conditions();
      for (size_t n = 0; n < nodes; ++n) {
        size_t ix = offset_x + 2 * n;
        size_t iy = offset_y + 2 * n;
        size_t iz = offset_z + 2 * n;
        a_prev(ix) = (-R0(ix)) / (mu * ds);
        a_prev(iy) = (-R0(iy)) / (mu * ds);
        a_prev(iz) = (-R0(iz)) / (mu * ds);
      }
    }

    u_prev = u;

    real_t S_norm = 0;

    for (size_t iter_outer = 0; iter_outer < max_iter_outer; iter_outer++) {
      assemble_system_newmark(dt, load, beta, gamma);
      apply_boundary_conditions();

      real_t res_norm = residual.norm();

      if (res_norm < tol_outer) {
        // std::cout << time_iter << " : " << iter_outer << " : ";
        // std::cout << "||r|| = " << res_norm << "\t";
        // std::cout << "||S|| = " << S_norm << std::endl;
        break;
      } else if (iter_outer == max_iter_outer - 1) {
        ELFF_ABORT(
          "EulerBeamDynamicInextensibleMoM::solve() did not converge.\n");
      } else {
        // std::cout << time_iter << " : " << iter_outer << " : ";
        // std::cout << "||r|| = " << res_norm << "\t";
        // std::cout << "||S|| = " << S_norm << std::endl;
      }

      solver.compute(jacobian);
      Eigen::VectorXd delta_u = solver.solve(-residual);
      u += delta_u;

      S_norm = update_lambda();
    }


    size_t nodes = mesh.get_nodes();

    for (size_t ni = 0; ni < nodes; ++ni) {
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


    update_mesh();

    time_iter++;
    t += dt;
  }

  void solve_newmark(real_t dt,
                     std::array<real_t, 3> load,
                     real_t beta,
                     real_t gamma)
  {

    // Eigen::LDLT<Eigen::MatrixXd> solver;
    // solver.setTolerance(tol_inner);
    Eigen::ConjugateGradient<Eigen::MatrixXd, Eigen::Upper | Eigen::Lower>
      solver;

    if (time_iter == 0) {
      Eigen::VectorXd R0 =
        EulerBeamStaticInextensibleMoM ::assemble_residual_template<real_t>(
          u, load);
      apply_boundary_conditions();
      for (size_t n = 0; n < nodes; ++n) {
        size_t ix = offset_x + 2 * n;
        size_t iy = offset_y + 2 * n;
        size_t iz = offset_z + 2 * n;
        a_prev(ix) = (-R0(ix)) / (mu * ds);
        a_prev(iy) = (-R0(iy)) / (mu * ds);
        a_prev(iz) = (-R0(iz)) / (mu * ds);
      }
    }

    u_prev = u;

    real_t S_norm = 0;

    for (size_t iter_outer = 0; iter_outer < max_iter_outer; iter_outer++) {
      assemble_system_newmark(dt, load, beta, gamma);
      apply_boundary_conditions();

      real_t res_norm = residual.norm();

      if (res_norm < tol_outer) {
        // std::cout << time_iter << " : " << iter_outer << " : ";
        // std::cout << "||r|| = " << res_norm << "\t";
        // std::cout << "||S|| = " << S_norm << std::endl;
        break;
      } else if (iter_outer == max_iter_outer - 1) {
        ELFF_ABORT(
          "EulerBeamDynamicInextensibleMoM::solve() did not converge.\n");
      } else {
        // std::cout << time_iter << " : " << iter_outer << " : ";
        // std::cout << "||r|| = " << res_norm << "\t";
        // std::cout << "||S|| = " << S_norm << std::endl;
      }

      solver.compute(jacobian);
      Eigen::VectorXd delta_u = solver.solve(-residual);
      u += delta_u;

      S_norm = update_lambda();
    }

    update_mesh();

    size_t nodes = mesh.get_nodes();

    for (size_t ni = 0; ni < nodes; ++ni) {
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

  void apply_initial_condition() override
  {
    EulerBeamStaticInextensibleMoM::apply_initial_condition();
    u_prev = u;
  }

  void apply_initial_condition(EulerBeamMesh& bmesh) override
  {
    EulerBeamStaticInextensibleMoM::apply_initial_condition(bmesh);
    u_prev = u;
  }

protected:
  Eigen::VectorXd v_prev;
  Eigen::VectorXd a_prev;
  Eigen::VectorXd u_prev, u_prev_prev;
  Eigen::MatrixXd mass;
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

  /**
   *
   */
  void assemble_system_newmark(real_t dt,
                               std::vector<std::array<real_t, 3>> load,
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
    real_t gamma) const
  {
    Eigen::Matrix<T, Eigen::Dynamic, 1> residual =
      EulerBeamStaticInextensibleMoM::assemble_residual_template<T>(u, load);

    if (!(dt > 0.0))
      throw std::runtime_error("Newmark: dt must be > 0");
    if (!(beta > 0.0))
      throw std::runtime_error("Newmark: beta must be > 0");

    const double inv = 1.0 / (beta * dt * dt);
    const double inv_bt = 1.0 / (beta * dt);
    const double kappa = (1.0 - 2.0 * beta) / (2.0 * beta);

    auto newmark_a = [&](Eigen::Index i) -> T {
      return inv * (u(i) - u_prev(i)) - inv_bt * v_prev(i) - kappa * a_prev(i);
    };

    for (size_t n = 0; n < nodes; ++n) {
      const Eigen::Index ix = static_cast<Eigen::Index>(offset_x + 2 * n);
      const Eigen::Index iy = static_cast<Eigen::Index>(offset_y + 2 * n);
      const Eigen::Index iz = static_cast<Eigen::Index>(offset_z + 2 * n);

      const T ax = newmark_a(ix);
      const T ay = newmark_a(iy);
      const T az = newmark_a(iz);

      residual(ix) += mu * ds * ax;
      residual(iy) += mu * ds * ay;
      residual(iz) += mu * ds * az;
    }

    return residual;
  }

  template<typename T>
  Eigen::Matrix<T, Eigen::Dynamic, 1> assemble_residual_newmark(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& u,
    real_t dt,
    std::vector<std::array<real_t, 3>> load,
    real_t beta,
    real_t gamma) const
  {
    Eigen::Matrix<T, Eigen::Dynamic, 1> residual =
      EulerBeamStaticInextensibleMoM::assemble_residual_template<T>(u, load);

    if (!(dt > 0.0))
      throw std::runtime_error("Newmark: dt must be > 0");
    if (!(beta > 0.0))
      throw std::runtime_error("Newmark: beta must be > 0");

    const double inv = 1.0 / (beta * dt * dt);
    const double inv_bt = 1.0 / (beta * dt);
    const double kappa = (1.0 - 2.0 * beta) / (2.0 * beta);

    auto newmark_a = [&](Eigen::Index i) -> T {
      return inv * (u(i) - u_prev(i)) - inv_bt * v_prev(i) - kappa * a_prev(i);
    };

    for (size_t n = 0; n < nodes; ++n) {
      const Eigen::Index ix = static_cast<Eigen::Index>(offset_x + 2 * n);
      const Eigen::Index iy = static_cast<Eigen::Index>(offset_y + 2 * n);
      const Eigen::Index iz = static_cast<Eigen::Index>(offset_z + 2 * n);

      const T ax = newmark_a(ix);
      const T ay = newmark_a(iy);
      const T az = newmark_a(iz);

      residual(ix) += mu * ds * ax;
      residual(iy) += mu * ds * ay;
      residual(iz) += mu * ds * az;
    }

    return residual;
  }

  void update_mesh()
  {

    std::vector<std::array<real_t, 3>>& centerline = mesh.get_centerline();
    std::vector<std::array<real_t, 3>>& velocity = mesh.get_centerline_velocity();
    std::vector<std::array<real_t, 3>>& slope = mesh.get_slope();
    std::vector<real_t>& s = mesh.get_curvilinear_axis();

    for (size_t i = 0; i < nodes; ++i) {
      centerline[i][0] = u(offset_x + 2 * i + 0);
      centerline[i][1] = u(offset_y + 2 * i + 0);
      centerline[i][2] = u(offset_z + 2 * i + 0);
      slope[i][0] = u(offset_x + 2 * i + 1);
      slope[i][1] = u(offset_y + 2 * i + 1);
      slope[i][2] = u(offset_z + 2 * i + 1);
      velocity[i][0] = v_prev(offset_x + 2 * i + 0);
      velocity[i][1] = v_prev(offset_y + 2 * i + 0);
      velocity[i][2] = v_prev(offset_z + 2 * i + 0);
    }
  }

};
}