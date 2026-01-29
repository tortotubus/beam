#pragma once

#include "EulerBeamStaticInextensibleADDM.hpp"

namespace ELFF {
class EulerBeamDynamicInextensibleADDM : public EulerBeamStaticInextensibleADDM 
{
public:
  Eigen::MatrixXd M;
  
  size_t previous_times;

  Eigen::MatrixXd xt;
  Eigen::MatrixXd yt;
  Eigen::MatrixXd zt;

  real_t dt;

  // Eigen::VectorXd f_x_dyn;
  // Eigen::VectorXd f_y_dyn;
  // Eigen::VectorXd f_z_dyn;

public:
  EulerBeamDynamicInextensibleADDM(
    real_t length,
    real_t EI,
    real_t mu,
    size_t nodes,
    EulerBeam::EulerBeamBCs bcs,
    real_t r_penalty
  ) : EulerBeamStaticInextensibleADDM(
    length,
    EI,
    mu,
    nodes,
    bcs,
    r_penalty
  )
  , M(Eigen::MatrixXd::Zero(dof,dof))
  , previous_times(2)
  , xt(Eigen::MatrixXd::Zero(dof,previous_times))
  , yt(Eigen::MatrixXd::Zero(dof,previous_times))
  , zt(Eigen::MatrixXd::Zero(dof,previous_times))
  {
    assemble_M();
  };

  void apply_initial_condition() override {
    // Compute initial condition from the time-indepdent base class
    EulerBeamStaticInextensibleADDM::apply_initial_condition();

    // Copy initial values to previous time(s)
    xt.col(0) = x;
    yt.col(0) = y;
    zt.col(0) = z;

    xt.col(1) = x;
    yt.col(1) = y;
    zt.col(1) = z;    
  }

  void assemble_M() {
    const real_t h = mesh.get_ds();
    const real_t xi_q[3] = {0.1127016654, 0.5, 0.8872983346};
    const real_t w_q [3] = {0.2777777778, 0.4444444444, 0.2777777778}; 

    for (size_t ei = 0; ei < elements; ++ei)
    {
      size_t edofs[4] = {
        2*(ei+0)+0,
        2*(ei+0)+1,
        2*(ei+1)+0,
        2*(ei+1)+1,
      };

      real_t Me[4][4] = {{0}};

      // Accumulate Me over quadrature points
      for (size_t qi = 0; qi < 3; ++qi) {
        real_t xi = xi_q[qi];
        real_t w = w_q[qi];

        auto H = CubicHermite<real_t>::values(xi,h);

        for (size_t a=0; a<4; ++a) {
          for (size_t b=0; b<4; ++b) {
            Me[a][b] += H[a]*H[b]*w*h;
          }
        }
      }

      // Scatter to the global matrix
      for (size_t a = 0; a < 4; ++a) {
        for (size_t b = 0; b < 4; ++b) {
          M(edofs[a],edofs[b]) += mu * Me[a][b];
        }
      }
    }
  }

  void assemble_A_dyn_BE() {
    A = M + (dt*dt) * A;
  }

  void assemble_f_dyn_BE() {
    f_x = M * (2.0*xt.col(0) - xt.col(1)) + (dt*dt) * f_x;
    f_y = M * (2.0*yt.col(0) - yt.col(1)) + (dt*dt) * f_y;
    f_z = M * (2.0*zt.col(0) - zt.col(1)) + (dt*dt) * f_z;
  }

  void update_pq() override { 
    EulerBeamStaticInextensibleADDM::update_pq();
  }

  void update_xy(std::array<real_t, 3> load) override { 
    // EulerBeamStaticInextensibleADDM::update_xy();

    std::cout << "Dynamic::update_xy()\n"; 

    assemble_f(load);
    assemble_f_dyn_BE();    
    apply_boundary_condition_f();

    x = llt.solve(f_x);
    y = llt.solve(f_y);
    z = llt.solve(f_z);
  }

  void update_multipliers() override { 
    EulerBeamStaticInextensibleADDM::update_multipliers();
  }

  void solve(std::array<real_t,3> load) override {
    for (size_t iter = 0; iter < max_outer; ++iter) {
      std::cout << "Iter " << iter << std::endl;
      update_pq();
      update_xy(load);
      update_multipliers();
      if (is_converged()) { break; }
    }

    update_mesh();
  };

  void solve(real_t dt, std::array<real_t,3> load) override {
    // Solve the static system with modified A, f
    solve(load);

    // Copy the x^{n} terms to x^{n-1}
    xt.col(1) = xt.col(0);
    yt.col(1) = yt.col(0);
    zt.col(1) = zt.col(0);

    // Copy the x^{n+1} terms to x^{n}
    xt.col(0) = x;
    yt.col(0) = y;
    zt.col(0) = z;
  }

};
}