#pragma once

#include "elff/models/beam/EulerBeamStaticInextensibleMoMSparse.hpp"

using namespace Eigen;

namespace ELFF {
namespace Models {

/**
 * @brief Dynamic inextensible Euler beam solved with a sparse
 * method-of-multipliers formulation and Newmark time integration.
 */
class EulerBeamDynamicInextensibleMoMSparse : public EulerBeamStaticInextensibleMoMSparse
{
public:
  /**
   * @brief Constructs a dynamic inextensible sparse beam model.
   *
   * @param length Beam length
   * @param EI Flexural rigidity
   * @param mu Mass per unit length
   * @param nodes Number of discretization nodes
   * @param bcs Boundary conditions at the beam ends
   * @param r_penalty Penalty parameter used in the inextensibility constraint
   */
  EulerBeamDynamicInextensibleMoMSparse(real_t length,
                                        real_t EI,
                                        real_t mu,
                                        size_t nodes,
                                        EulerBeam::EulerBeamBCs bcs,
                                        real_t r_penalty);

  /**
   * @brief Advances the beam one time step under a uniform load.
   *
   * @param dt Time-step size
   * @param load Uniform distributed load vector
   */
  virtual void solve(real_t dt, std::array<real_t, 3> load) override;

  /**
   * @brief Advances the beam one time step under nodal loads.
   *
   * @param dt Time-step size
   * @param load Load vector specified at the mesh nodes
   */
  virtual void solve(real_t dt,
                     std::vector<std::array<real_t, 3>> load) override;

  /**
   * @brief Solves one Newmark step for uniform loading.
   *
   * @param dt Time-step size
   * @param load Uniform distributed load vector
   * @param beta Newmark beta parameter
   * @param gamma Newmark gamma parameter
   */
  void solve_newmark(real_t dt,
                     std::array<real_t, 3> load,
                     real_t beta,
                     real_t gamma);

  /**
   * @brief Solves one Newmark step for nodal loading.
   *
   * @param dt Time-step size
   * @param load Load vector specified at the mesh nodes
   * @param beta Newmark beta parameter
   * @param gamma Newmark gamma parameter
   */
  void solve_newmark(real_t dt,
                     std::vector<std::array<real_t, 3>> load,
                     real_t beta,
                     real_t gamma);

  /**
   * @brief Applies the model's default initial condition.
   */
  void apply_initial_condition() override;

  /**
   * @brief Applies initial conditions from a supplied beam mesh.
   *
   * @param bmesh Beam mesh containing the initial geometry
   */
  void apply_initial_condition(EulerBeamMesh& bmesh) override;

protected:
  /**
   * @brief Velocity degrees of freedom from the previous time step.
   */
  VectorXd v_prev;
  /**
   * @brief Acceleration degrees of freedom from the previous time step.
   */
  VectorXd a_prev;
  /**
   * @brief Displacement history used by the time integrator.
   */
  VectorXd u_prev;
  /**
   * @brief Mass-matrix storage for dynamic formulations.
   */
  MatrixXd mass;
  /**
   * @brief Cached load from the previous time step.
   */
  std::array<real_t, 3> load_prev;

  /**
   * @brief Assembles the Newmark system for uniform loading.
   *
   * @param dt Time-step size
   * @param load Uniform distributed load vector
   * @param beta Newmark beta parameter
   * @param gamma Newmark gamma parameter
   */
  void assemble_system_newmark(real_t dt,
                               std::array<real_t, 3> load,
                               real_t beta,
                               real_t gamma);

  /**
   * @brief Assembles the Newmark system for nodal loading.
   *
   * @param dt Time-step size
   * @param load Load vector specified at the mesh nodes
   * @param beta Newmark beta parameter
   * @param gamma Newmark gamma parameter
   */
  void assemble_system_newmark(real_t dt,
                               std::vector<std::array<real_t, 3>> load,
                               real_t beta,
                               real_t gamma);

  /**
   * @brief Adds the analytically assembled inertial residual and tangent
   * contributions for the Newmark update.
   *
   * @param dt Time-step size
   * @param beta Newmark beta parameter
   */
  void add_newmark_inertial_terms(real_t dt, real_t beta);

  /**
   * @brief Updates the beam mesh positions, slopes, and velocities from the
   * current state vectors.
   */
  void update_mesh();

  template<typename T>
  /**
   * @brief Assembles the Newmark residual for uniform loading.
   *
   * @tparam T Scalar type used for residual assembly
   * @param u State vector at which to evaluate the residual
   * @param dt Time-step size
   * @param load Uniform distributed load vector
   * @param beta Newmark beta parameter
   * @param gamma Newmark gamma parameter
   * @return Residual vector for the Newmark update
   */
  Matrix<T, Dynamic, 1> assemble_residual_newmark(
    const Matrix<T, Dynamic, 1>& u,
    real_t dt,
    const std::array<real_t, 3> load,
    real_t beta,
    real_t gamma) const
  {
    Matrix<T, Dynamic, 1> res =
      assemble_residual_template<T>(u, load);

    // Guards to avoid NaNs/Infs:
    if (!(dt > 0.0))
      throw std::runtime_error("Newmark: dt must be > 0");
    if (!(beta > 0.0))
      throw std::runtime_error("Newmark: beta must be > 0");

    // Keep these as doubles (scalars), not T:
    const double inv = 1.0 / (beta * dt * dt); // = 1/(β Δt²)
    const double inv_bt = 1.0 / (beta * dt);   // = 1/(β Δt)
    const double kappa = (1.0 - 2.0 * beta) / (2.0 * beta);

    auto newmark_a = [&](Index i) -> T {
      // u(i) is AD (T); u_prev/v_prev/a_prev are doubles
      return inv * (u(i) - u_prev(i)) - inv_bt * v_prev(i) - kappa * a_prev(i);
    };

    for (size_t n = 0; n < nodes; ++n) {
      const Index ix = static_cast<Index>(offset_x + 2 * n);
      const Index iy = static_cast<Index>(offset_y + 2 * n);
      const Index iz = static_cast<Index>(offset_z + 2 * n);

      const T ax = newmark_a(ix);
      const T ay = newmark_a(iy);
      const T az = newmark_a(iz);
      const real_t w = (n == 0 || n == nodes - 1) ? 0.5 * ds : ds;

      res(ix) += mu * w * ax;
      res(iy) += mu * w * ay;
      res(iz) += mu * w * az;
    }

    return res;
  }

  template<typename T>
  /**
   * @brief Assembles the Newmark residual for nodal loading.
   *
   * @tparam T Scalar type used for residual assembly
   * @param u State vector at which to evaluate the residual
   * @param dt Time-step size
   * @param load Load vector specified at the mesh nodes
   * @param beta Newmark beta parameter
   * @param gamma Newmark gamma parameter
   * @return Residual vector for the Newmark update
   */
  Matrix<T, Dynamic, 1> assemble_residual_newmark(
    const Matrix<T, Dynamic, 1>& u,
    real_t dt,
    const std::vector<std::array<real_t, 3>> load,
    real_t beta,
    real_t gamma) const
  {
    Matrix<T, Dynamic, 1> res =
      assemble_residual_template<T>(u, load);

    if (!(dt > 0.0))
      throw std::runtime_error("Newmark: dt must be > 0");
    if (!(beta > 0.0))
      throw std::runtime_error("Newmark: beta must be > 0");

    const double inv = 1.0 / (beta * dt * dt);
    const double inv_bt = 1.0 / (beta * dt);
    const double kappa = (1.0 - 2.0 * beta) / (2.0 * beta);

    auto newmark_a = [&](Index i) -> T {
      return inv * (u(i) - u_prev(i)) - inv_bt * v_prev(i) - kappa * a_prev(i);
    };

    for (size_t n = 0; n < nodes; ++n) {
      const Index ix = static_cast<Index>(offset_x + 2 * n);
      const Index iy = static_cast<Index>(offset_y + 2 * n);
      const Index iz = static_cast<Index>(offset_z + 2 * n);

      const T ax = newmark_a(ix);
      const T ay = newmark_a(iy);
      const T az = newmark_a(iz);
      const real_t w = (n == 0 || n == nodes - 1) ? 0.5 * ds : ds;

      res(ix) += mu * w * ax;
      res(iy) += mu * w * ay;
      res(iz) += mu * w * az;
    }

    return res;
  }
};

} // namespace Models 
} // namespace ELFF 
