#pragma once

#include "EulerBeamStaticInextensibleMoMSparse.hpp"

namespace beam {
class EulerBeamDynamicInextensibleMoMSparse : public EulerBeamStaticInextensibleMoMSparse
{
public:
  EulerBeamDynamicInextensibleMoMSparse(real_t length,
                                  real_t EI,
                                  real_t mu,
                                  size_t nodes,
                                  beam::EulerBeamBCs bcs,
                                  real_t r_penalty)
    : EulerBeamStaticInextensibleMoMSparse(length, EI, mu, nodes, bcs, r_penalty)
    , u_prev(Eigen::VectorXd::Zero(ndof))
    , u_prev_prev(Eigen::VectorXd::Zero(ndof))
    , v_prev(Eigen::VectorXd::Zero(ndof))
    , a_prev(Eigen::VectorXd::Zero(ndof)) 
    {};

  virtual void solve(real_t dt, std::array<real_t, 3> load) override
  {
    const real_t alpha = 0.0;
    // const real_t alpha = -0.1;
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
      Eigen::VectorXd R0 = EulerBeamStaticInextensibleMoMSparse ::assemble_residual_template<real_t>(u, load);

      for (size_t n = 0; n < nodes; ++n) {
        size_t ix = offset_x + 2 * n;
        size_t iy = offset_y + 2 * n;
        size_t iz = offset_z + 2 * n;
        a_prev(ix) = (-R0(ix)) / (mu*ds);
        a_prev(iy) = (-R0(iy)) / (mu*ds);
        a_prev(iz) = (-R0(iz)) / (mu*ds);
      }
    }

    u_prev = u;

    Eigen::ConjugateGradient<
      Eigen::SparseMatrix<real_t>, 
      Eigen::Lower | Eigen::Upper,
      Eigen::IncompleteCholesky<real_t>
    > solver;

    real_t S_norm = 0;

    for (size_t iter_outer = 0; iter_outer < max_iter_outer; iter_outer++) {
      assemble_system_newmark(dt, load, beta, gamma);
      apply_boundary_conditions();

      real_t res_norm = residual.norm();

      if (res_norm < tol_outer) {
        // std::cout << iter_outer << ": ||r|| = " << res_norm;
        // std::cout << "\t ||S|| = " << S_norm << std::endl;
        break;
      } else if (iter_outer == max_iter_outer - 1) {
        // BEAM_ABORT( "EulerBeamDynamicInextensibleMoMSparse::solve() did not converge.\n");
      } else {
        // std::cout << iter_outer << ": ||r|| = " << res_norm;
        // std::cout << "\t ||S|| = " << S_norm << std::endl;
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
    EulerBeamStaticInextensibleMoMSparse::apply_initial_condition();
    u_prev = u;
  }

  void apply_initial_condition(EulerBeamMesh& bmesh) override
  {
    EulerBeamStaticInextensibleMoMSparse::apply_initial_condition(bmesh);
    u_prev = u;
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
    using Tpl = Eigen::Triplet<real_t>;

    // --- 1) Build the AutoDiff input vector x_ad ---
    ADVec x_ad(ndof);
    for (int i = 0; i < ndof; ++i) {
      Eigen::VectorXd seed = Eigen::VectorXd::Zero(ndof);
      seed(i) = 1.0;
      x_ad(i) = AD(u(i), seed);
    }

    // --- 2) Compute AD residuals ---
    ADVec R_ad = assemble_residual_newmark<AD>(x_ad, dt, load, beta, gamma);

    // --- 3) Extract values + build Jacobian triplets ---
    residual.resize(ndof);

    // If you know roughly how many nonzeros per row, you can reserve:
    std::vector<Tpl> triplets;
    triplets.reserve(ndof * 5); // e.g. assume 5 nnz/row on average

    for (int i = 0; i < ndof; ++i) {
      residual(i) = R_ad(i).value();

      // Access the derivative vector for row i
      const Eigen::VectorXd& dRi = R_ad(i).derivatives();
      const int nnz = static_cast<int>(dRi.size());

      // OPTION A: Filter zeros on the fly
      for (int j = 0; j < nnz; ++j) {
        const real_t dj = dRi[j];
        if (dj != 0.0) {
          triplets.emplace_back(i, j, dj);
        }
      }

      /*
      // OPTION B: If you have a precomputed sparsity pattern
      //   (e.g. std::vector<std::vector<int>> sparsity where
      //    sparsity[i] lists the column-indices that can be nonzero in row i)
      for (int j : sparsity[i]) {
        real_t dj = dRi[j];
        if (dj != 0.0)
          triplets.emplace_back(i, j, dj);
      }
      */
    }

    // --- 4) Assemble the sparse matrix ---
    jacobian.resize(ndof, ndof);
    jacobian.setFromTriplets(triplets.begin(), triplets.end());
    jacobian.makeCompressed();
  }

  template<typename T>
  Eigen::Matrix<T, Eigen::Dynamic, 1> assemble_residual_newmark(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& u,
    real_t dt,
    const std::array<real_t, 3> load,
    real_t beta,
    real_t gamma) const
  {
    Eigen::Matrix<T, Eigen::Dynamic, 1> res = assemble_residual_template<T>(u, load);

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

    for (size_t n = 0; n < nodes; ++n) {
      const Eigen::Index ix = static_cast<Eigen::Index>(offset_x + 2 * n);
      const Eigen::Index iy = static_cast<Eigen::Index>(offset_y + 2 * n);
      const Eigen::Index iz = static_cast<Eigen::Index>(offset_z + 2 * n);

      const T ax = newmark_a(ix);
      const T ay = newmark_a(iy);
      const T az = newmark_a(iz);

      res(ix) += (mu*ds) * ax; // mu is double; AD ⊗ scalar is fine
      res(iy) += (mu*ds) * ay;
      res(iz) += (mu*ds) * az;
    }

    return res;
  }
};
}