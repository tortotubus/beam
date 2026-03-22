#pragma once

#include "models/beam/EulerBeamStaticInextensibleMoMSparse.hpp"

using namespace Eigen;

namespace ELFF {
namespace Models {

class EulerBeamDynamicInextensibleMoMSparse : public EulerBeamStaticInextensibleMoMSparse
{
public:
  EulerBeamDynamicInextensibleMoMSparse(real_t length,
                                        real_t EI,
                                        real_t mu,
                                        size_t nodes,
                                        EulerBeam::EulerBeamBCs bcs,
                                        real_t r_penalty);

  virtual void solve(real_t dt, std::array<real_t, 3> load) override;

  virtual void solve(real_t dt,
                     std::vector<std::array<real_t, 3>> load) override;

  void solve_newmark(real_t dt,
                     std::array<real_t, 3> load,
                     real_t beta,
                     real_t gamma);

  void solve_newmark(real_t dt,
                     std::vector<std::array<real_t, 3>> load,
                     real_t beta,
                     real_t gamma);

  void apply_initial_condition() override;

  void apply_initial_condition(EulerBeamMesh& bmesh) override;

protected:
  VectorXd v_prev;
  VectorXd a_prev;
  VectorXd u_prev;
  MatrixXd mass;
  std::array<real_t, 3> load_prev;

  /**
   *
   */
  void assemble_system_newmark(real_t dt,
                               std::array<real_t, 3> load,
                               real_t beta,
                               real_t gamma);

  void assemble_system_newmark(real_t dt,
                               std::vector<std::array<real_t, 3>> load,
                               real_t beta,
                               real_t gamma);

  void update_mesh();

  template<typename T>
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
