#pragma once

// #include <beam/LinAlg/Matrix.hpp>
// #include <beam/LinAlg/Vector.hpp>

#include "elff/models/beam/EulerBeam.hpp"
#include "elff/fem/Shapes.hpp"

#include <cstdio> // for popen, pclose, fprintf
#include <iostream>
#include <stdlib.h>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/AutoDiff>

using namespace Eigen;

namespace ELFF {
namespace Models {

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
 *    \mathcal{L}_R(x,y,p,q,\lambda_x,\mu) = J(x,y) + \int_0^L
 * \left[\lambda_x(p-x')
 * + \mu(q-y')\right] ds + \frac{r}{2} \int_0^L \left[(p-x')^2 + (q-y')^2\right]
 * ds.
 * \f]
 */

class EulerBeamStaticInextensibleADDM : public EulerBeam
{

public:
  EulerBeamStaticInextensibleADDM(real_t length,
                                  real_t EI,
                                  size_t nodes,
                                  EulerBeamBCs bcs,
                                  real_t r_penalty);

  ~EulerBeamStaticInextensibleADDM();

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

  const VectorXd& get_lambda_x() const;
  const VectorXd& get_lambda_y() const;
  const VectorXd& get_lambda_z() const;
  const MatrixXd& get_A() const;

protected:
  size_t dimension;
  size_t elements;
  size_t dof;
  real_t r_penalty, alpha;

  MatrixXd A;
  MatrixXd A_unconstrained;
  VectorXd x, y, z;
  VectorXd f_x, f_y, f_z;

  LLT<MatrixXd> llt;

  VectorXd lambda_x, lambda_y, lambda_z;
  VectorXd p, q, r;
  VectorXd xp, yp, zp;

  size_t max_outer;
  real_t tol_outer;

  EulerBeamStaticInextensibleADDM(real_t length,
                                  real_t EI,
                                  real_t mu,
                                  size_t nodes,
                                  EulerBeamBCs bcs,
                                  real_t r_penalty);

  void update_mesh();

  void apply_initial_condition_xy();

  void apply_initial_condition_xy(EulerBeamMesh& bmesh);

  void compute_slopes_collocation();

  void apply_initial_condition_pq();

  void apply_boundary_condition_pq();

  void apply_boundary_condition_lambda();

  virtual void update_pq();

  void apply_boundary_condition_A();

  void assemble_A();

  void decompose_A();

  void assemble_f(std::array<real_t, 3> load);

  void apply_boundary_condition_f();

  virtual void update_xy(std::array<real_t, 3> load);

  virtual void update_multipliers();

  bool is_converged(bool recompute_slopes = true);
};


} // namespace Models 
} // namespace ELFF 
