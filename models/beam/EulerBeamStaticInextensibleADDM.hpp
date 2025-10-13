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

/**
 * @brief A class to solve the static inextensible Euler–Bernoulli beam
 * equation; this model is valid for large deflections. The strong form of the
 * system that we seek to solve is
 * \f[
 *    \frac{\partial^2}{\partial s^2} \left(EI \frac{\partial^2
 * \mathbf{r}}{\partial s^2}\right) = \mathbf{q}(s)
 * \f]
 * where \f(EI\f) is a bending stiffness, \f(\mathbf{r}(s) = (x(s),y(s))\f) is
 * the deflection of our beam. Unlike the classic Euler-Bernouli beam equation,
 * where inextensibility is implicitly enforced, in the \f(n > 1\f) dimensional
 * version, we consider the deflection of the beam in each dimension along a
 * curvilinear coordinate system, and our inextensbility must be enforced
 * explicitly. Then, in addition, we enforce on the solution a pointwise
 * constraint
 * \f[
 *    ||\mathbf{r}'(s)||^2 = 1.
 * \f]
 * To derive a weak form, we introduce a smooth test function
 * ...
 * \f[
 *    \mathcal{L}_R(x,y,p,q,\lambda,\mu) = J(x,y) + \int_0^L \left[\lambda(p-x')
 * + \mu(q-y')\right] ds + \frac{r}{2} \int_0^L \left[(p-x')^2 + (q-y')^2\right]
 * ds.
 * \f]
 */

class EulerBeamStaticInextensibleADDM : public EulerBeam
{
public:
  size_t dimension;
  size_t elements;
  size_t dof;
  real_t r_penalty, alpha;

  Eigen::MatrixXd A; 
  Eigen::VectorXd x, y;
  Eigen::VectorXd f_x, f_y;

  Eigen::LLT<Eigen::MatrixXd> llt;

  Eigen::VectorXd lambda, mu;
  Eigen::VectorXd p, q;
  Eigen::VectorXd xp, yp;

  size_t max_outer;
  real_t tol_outer;

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
    , dof((2 * nodes))

    , A(Eigen::MatrixXd::Zero(dof, dof))
    , f_x(Eigen::VectorXd::Zero(dof))
    , f_y(Eigen::VectorXd::Zero(dof))
    , x(Eigen::VectorXd::Zero(dof))
    , y(Eigen::VectorXd::Zero(dof))
    , llt()

    , xp(Eigen::VectorXd::Zero(nodes + elements))
    , yp(Eigen::VectorXd::Zero(nodes + elements))
    , lambda(Eigen::VectorXd::Zero(nodes + elements))
    , mu(Eigen::VectorXd::Zero(nodes + elements))
    , p(Eigen::VectorXd::Ones(nodes + elements))
    , q(Eigen::VectorXd::Zero(nodes + elements))

    , r_penalty(r_penalty)
    , tol_outer(1e-6)
    , max_outer(10000)
  {
    set_load({ 0., load, 0. });
  };

  ~EulerBeamStaticInextensibleADDM() {};

  const Eigen::VectorXd& get_lambda() const { return lambda; }
  const Eigen::VectorXd& get_mu() const { return mu; }
  const Eigen::MatrixXd& get_A() const { return A; }

  // void set_lambda(Eigen::VectorXd l) { this->lambda = l; }
  // void set_mu(Eigen::VectorXd m) { this->mu = m; }
  // void set_A(Eigen::MatrixXd A) { this->A = A; }

  // Update this->mesh object based on the solution
  void update_mesh() {
    size_t nodes = this->mesh.get_nodes();

    size_t ndof_x = 2 * nodes;
    size_t ndof_y = 2 * nodes;

    std::vector<std::array<real_t, 3>>& centerline = mesh.get_centerline();
    std::vector<real_t>& s = mesh.get_curvilinear_axis();

    for (size_t i = 0; i < nodes; ++i) {
      centerline[i][0] = x(2 * i);
      centerline[i][1] = y(2 * i);
      centerline[i][2] = 0.;
    }
  }

  void apply_initial_condition() {
    apply_initial_condition_xy();
    apply_initial_condition_pq();
  }

  // Set initial condition
  void apply_initial_condition_xy() {
    real_t h = mesh.get_ds();
    size_t nodes = mesh.get_nodes();

    for (size_t i = 0; i < nodes; i++) {
      x(2*i + 0) = h*i;
      x(2*i + 1) = 1.;
      y(2*i + 0) = 0.;
      y(2*i + 1) = 0.;
    }
  }

  //
  void compute_slopes_collocation() {
    xp.setZero(); yp.setZero();
    const real_t h = this->mesh.get_ds();
    const size_t nodes = this->mesh.get_nodes();

    // Node collocation: Take the slope directly from the xy-solution
    for (size_t ni = 0; ni < nodes; ni++) {
      xp[ni] = x[2*ni + 1];
      yp[ni] = y[2*ni + 1];
    }

    // Midpoint colloation: Evaulate slope via Hermite derivatives at xi = 0.5
    for (size_t ei = 0; ei < elements; ei++) {
      std::array<real_t,4> dH = CubicHermite<real_t>::derivs(0.5,h);

      size_t edofs[4] = {
        2*(ei+0)+0, 
        2*(ei+0)+1, 
        2*(ei+1)+0, 
        2*(ei+1)+1
      };

      xp[nodes + ei] = dH[0]*x[edofs[0]] + dH[1]*x[edofs[1]] + dH[2]*x[edofs[2]] + dH[3]*x[edofs[3]];
      yp[nodes + ei] = dH[0]*y[edofs[0]] + dH[1]*y[edofs[1]] + dH[2]*y[edofs[2]] + dH[3]*y[edofs[3]];
    }
  }

  void apply_initial_condition_pq() {
    const size_t nodes = this->mesh.get_nodes();
    compute_slopes_collocation();
    for (size_t ci = 0; ci < nodes + elements; ci++) {
      real_t xprime = xp[ci];
      real_t yprime = yp[ci];
      real_t den = std::max(1e-14,sqrt(xprime*xprime + yprime*yprime));
      p[ci] = xprime / den;
      q[ci] = yprime / den;
    }
    apply_boundary_condition_pq();
  }

  // Apply boundary conditions to the matrix A
  void apply_boundary_condition_A() {
    size_t nodes = mesh.get_nodes();
    real_t h = mesh.get_ds();

    for (size_t bi = 0; bi < 2; ++bi) {
      EulerBeamBCEnd bcend = boundary_conditions.end[bi];
      size_t ni = 0;

      std::vector<size_t> idx;   //idx.reserve(2);   // (<=2);
      std::vector<real_t> xvals; //xvals.reserve(2); // (<=2);
      std::vector<real_t> yvals; //yvals.reserve(2); // (<=2);

      switch (bcend)
      {
        case left:
          ni = 0;
          break;
        case right:
          ni = nodes - 1;
          break;
      }

      EulerBeamBCType bctype = boundary_conditions.type[bi];
      EulerBeamBCVals bcvals = boundary_conditions.vals[bi];

      switch(bctype) {
        case free_bc:
          idx = {};
          xvals = {};
          yvals = {};
          break;
        case simple_bc:
          idx = {2 * ni + 0};
          xvals = {bcvals.position[0]};
          yvals = {bcvals.position[1]};
          break;
        case clamped_bc:
          idx = {2 * ni + 0};//, 2 * ni + 1};
          xvals = {bcvals.position[0]}; //, bcvals.slope[0]};
          yvals = {bcvals.position[1]}; //, bcvals.slope[1]};
          break;
      }

      // Set boundary conditions in the matrix A
      for (size_t i = 0; i < xvals.size(); i++) {
        for (size_t j = 0; j < dof; ++ j) {
          A(idx[i],j) = 0.;
          A(j,idx[i]) = 0.;
        }
        A(idx[i], idx[i]) = 1.;
      }
    }
  }

  // Assemble the matrix A
  void assemble_A()
  {
    const real_t h = this->mesh.get_ds();
    const real_t xi_q[3] = { 0.1127016654, 0.5, 0.8872983346 };
    const real_t w_q[3] = { 0.2777777778, 0.4444444444, 0.2777777778 };

    // Local stiffness matrices
    Eigen::MatrixXd K4 = Eigen::MatrixXd::Zero(this->dof, this->dof);
    Eigen::MatrixXd K2 = Eigen::MatrixXd::Zero(this->dof, this->dof);

    for (size_t e = 0; e < this->elements; ++e) {
      const size_t edofs[4] = {
        2 * (e + 0) + 0, 
        2 * (e + 0) + 1, 
        2 * (e + 1) + 0, 
        2 * (e + 1) + 1
      };

      real_t K4e[4][4] = { { 0 } };
      real_t K2e[4][4] = { { 0 } };

      // Quadtrature
      for (size_t q = 0; q < 3; ++q) {
        const real_t xi = xi_q[q];
        const real_t w = w_q[q];

        auto H = CubicHermite<real_t>::values(xi, h);
        auto dH = CubicHermite<real_t>::derivs(xi, h);
        auto ddH = CubicHermite<real_t>::second_derivs(xi, h);

        // Accumulate local K4e and K2e
        for (size_t a = 0; a < 4; ++a) {
          for (size_t b = 0; b < 4; ++b) {
            K4e[a][b] += ddH[a] * ddH[b] * w * h;
            K2e[a][b] +=  dH[a] *  dH[b] * w * h;
          }
        }
      }

      // Scatter to global triplets
      for (size_t a = 0; a < 4; ++a) {
        for (size_t b = 0; b < 4; ++b) {
          K4(edofs[a], edofs[b]) += K4e[a][b];
          K2(edofs[a], edofs[b]) += K2e[a][b];
        }
      }
    }

    this->A = this->EI * K4 + this->r_penalty * K2;
  }

  // Cholesky decompose the matrix A
  void decompose_A() { 
    llt.compute(A); 
  }

  // Assemble the rhs vector for x-axis
  void assemble_f()
  {
    f_x.setZero();
    f_y.setZero();

    const real_t h = this->mesh.get_ds();
    const size_t nodes = this->mesh.get_nodes();

    const real_t xi_q[3] = { 0.1127016654, 0.5, 0.8872983346 };
    const real_t w_q[3] = { 0.2777777778, 0.4444444444, 0.2777777778 };

    for (size_t e = 0; e < elements; ++e) {
      size_t edofs[4] = {
        2 * (e + 0) + 0, 
        2 * (e + 0) + 1, 
        2 * (e + 1) + 0, 
        2 * (e + 1) + 1
      };
      real_t fxe[4] = { 0, 0, 0, 0 };
      real_t fye[4] = { 0, 0, 0, 0 };

      for (size_t qi = 0; qi < 3; ++qi) {
        real_t xi = xi_q[qi];
        real_t w = w_q[qi];

        auto L = QuadraticLagrange<real_t>::values(xi);
        auto H = CubicHermite<real_t>::values(xi, h);
        auto dH = CubicHermite<real_t>::derivs(xi, h);

        size_t l = e;
        size_t m = nodes + e;
        size_t r = e + 1;

        real_t p_val      = L[0] * p[l]      + L[1] * p[m]      + L[2] * p[r];
        real_t lambda_val = L[0] * lambda[l] + L[1] * lambda[m] + L[2] * lambda[r];
        real_t q_val      = L[0] * q[l]      + L[1] * q[m]      + L[2] * q[r];
        real_t mu_val     = L[0] * mu[l]     + L[1] * mu[m]     + L[2] * mu[r];

        for (size_t a = 0; a < 4; ++a) {
          fxe[a] += (lambda_val + r_penalty * p_val) * dH[a] * w * h;
          fye[a] += (mu_val     + r_penalty * q_val) * dH[a] * w * h;
          fxe[a] += uniform_load[0] * H[a] * w * h;
          fye[a] += uniform_load[1] * H[a] * w * h;
        }
      }
      for (size_t a = 0; a < 4; ++a) {
        f_x(edofs[a]) += fxe[a];
        f_y(edofs[a]) += fye[a];
      }
    }
  }

  // Apply the boundary condition to the f vector(s)
  void apply_boundary_condition_f() {
    size_t nodes = mesh.get_nodes();
    real_t h = mesh.get_ds();

    for (size_t bi = 0; bi < 2; ++bi) {
      EulerBeamBCEnd bcend = boundary_conditions.end[bi];
      size_t ni = 0;

      std::vector<size_t> idx;   //idx.reserve(2);   // (<=2);
      std::vector<real_t> xvals; //xvals.reserve(2); // (<=2);
      std::vector<real_t> yvals; //yvals.reserve(2); // (<=2);

      switch (bcend)
      {
        case left:
          ni = 0;
          break;
        case right:
          ni = nodes - 1;
          break;
      }

      EulerBeamBCType bctype = boundary_conditions.type[bi];
      EulerBeamBCVals bcvals = boundary_conditions.vals[bi];

      switch(bctype) {
        case free_bc:
          idx = {};
          xvals = {};
          yvals = {};
          break;
        case simple_bc:
          idx = {2 * ni + 0};
          xvals = {bcvals.position[0]};
          yvals = {bcvals.position[1]};
          break;
        case clamped_bc:
          idx = {2 * ni + 0 }; // 2 * ni + 1};
          xvals = {bcvals.position[0]}; //, bcvals.slope[0]};
          yvals = {bcvals.position[1]}; //, bcvals.slope[1]};
          break;
      }

      // Set boundary conditions in the matrix A
      for (size_t i = 0; i < xvals.size(); i++) {
        // Set boundary conditions in f_x 
        f_x(idx[i]) = xvals[i];
        // Set boundary conditions in f_y 
        f_y(idx[i]) = yvals[i];
      }
    }
  }

  // Apply the boundary condition to the p,q variables
  void apply_boundary_condition_pq() {
    const size_t nodes = mesh.get_nodes();
    for (size_t bi = 0; bi < 2; bi++) {
      if (boundary_conditions.type[bi] == clamped_bc) {
        switch (boundary_conditions.end[bi])
        {
          case left:
            p(0) = boundary_conditions.vals[bi].slope[0];
            q(0) = boundary_conditions.vals[bi].slope[1];
            break;
          case right:
            p(nodes-1) = boundary_conditions.vals[bi].slope[0];
            q(nodes-1) = boundary_conditions.vals[bi].slope[1];
            break;
        }
      }
    }
  }

  // Perform one update of the linear system
  void update_xy() {
    assemble_f();
    apply_boundary_condition_f();

    x = llt.solve(f_x);
    y = llt.solve(f_y);
  }

  // Perform a pointwise projection of (p^n, q^n) onto the unit circle to obtain (p^{n+1}, q^{n+1}) 
  void update_pq() {
    const real_t h = mesh.get_ds();
    const size_t nodes = mesh.get_nodes();
    compute_slopes_collocation();

    //#pragma omp parallel for schedule(static)
    for (size_t ci = 0; ci < nodes + elements; ci++) {
      const real_t xprime = xp[ci];
      const real_t yprime = yp[ci];
      real_t X = xprime - lambda[ci] / r_penalty; 
      real_t Y = yprime -     mu[ci] / r_penalty; 
      real_t norm = std::max(1e-6,sqrt(X*X + Y*Y));             
      p[ci] = X/norm;
      q[ci] = Y/norm;
    }

    apply_boundary_condition_pq();
  }

  // Update the lagrange multipliers \lambda and \mu
  void update_multipliers() {
    const real_t h = mesh.get_ds();
    const size_t nodes = mesh.get_nodes();
    compute_slopes_collocation();

    //#pragma omp parallel for schedule(static)
    for (size_t ci = 0; ci < nodes + elements; ci++) {
      lambda[ci] += r_penalty * (p[ci] - xp[ci]);
      mu[ci]     += r_penalty * (q[ci] - yp[ci]);
    }
  }

  // Solve 
  void solve() {
    apply_initial_condition_xy();
    apply_initial_condition_pq();

    assemble_A();
    apply_boundary_condition_A();
    decompose_A();

    for (size_t iter = 0; iter < max_outer; ++iter) {
      // std::cout << "Iter " << iter << std::endl;
      update_pq();
      update_xy();
      update_multipliers();
      if (is_converged()) { break; }
    }

    update_mesh();
  };

  // compute the inf norm ||\mathbf{p} - \mathbf{x}'||_{\infty} and ||\mathbf{q} - \mathbf{y}'||_{\infty}
bool is_converged(bool recompute_slopes = true) {
  if (recompute_slopes) {
    compute_slopes_collocation();   // fills xp, yp from current x,y
  }
  const real_t res_p = (p - xp).cwiseAbs().maxCoeff();
  const real_t res_q = (q - yp).cwiseAbs().maxCoeff();
  const real_t res   = std::max(res_p, res_q);

  // std::cout << "max|p-x'|=" << res_p << ", max|q-y'|=" << res_q << '\n';
  return res < tol_outer;
}

};

} // namespace beam