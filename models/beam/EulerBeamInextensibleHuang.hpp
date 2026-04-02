#pragma once

#include "elff/models/beam/EulerBeam.hpp"

#include <array>
#include <vector>

#include <Eigen/Dense>

using namespace Eigen;

namespace ELFF {
namespace Models {

/**
 * @brief Huang-style finite-difference prototype for the dynamic inextensible
 * beam benchmark.
 *
 * This prototype mirrors the no-fluid Huang benchmark structure and applies
 * endpoint closures directly from the supplied EulerBeamBCs on each end. The
 * currently supported endpoint types are free, simply supported, and clamped.
 * - staggered tension solve
 * - explicit Huang-style bending by default, with optional implicit bending
 *   for additional stabilization
 *
 * The implementation is intentionally compact and specialized to the 2D
 * hanging-filament benchmark rather than the full generality of the other
 * beam solvers in the codebase.
 */
class EulerBeamInextensibleHuang : public EulerBeam
{
public:
  using MatX3 = Matrix<real_t, Dynamic, 3>;
  using Vec3 = Matrix<real_t, 3, 1>;

  EulerBeamInextensibleHuang(real_t length,
                 real_t EI,
                 real_t mu,
                 size_t nodes,
                 EulerBeam::EulerBeamBCs bcs);

  void set_implicit_bending(bool enabled);
  bool get_implicit_bending() const;

  void solve(std::array<real_t, 3> load) override;
  void solve(std::vector<std::array<real_t, 3>> load) override;
  void solve(real_t dt, std::array<real_t, 3> load) override;
  void solve(real_t dt, std::vector<std::array<real_t, 3>> load) override;

  void apply_initial_condition();
  void apply_initial_condition(EulerBeamMesh& bmesh) override;

  const MatX3& get_initial_positions() const;
  const MatX3& get_positions() const;
  const VectorXd& get_last_tension() const;

protected:
  real_t ds;
  real_t tol_inner;
  real_t last_dt;
  bool initial_velocity_pending;
  bool implicit_bending;

  MatX3 X_init;
  MatX3 X_n;
  MatX3 X_nm1;
  MatX3 V_init;
  VectorXd last_T;
  MatrixXd bending_matrix;

  Index bc_index(EulerBeam::EulerBeamBCEnd end) const;
  EulerBeam::EulerBeamBCType bc_type(EulerBeam::EulerBeamBCEnd end) const;
  const EulerBeam::EulerBeamBCVals& bc_vals(EulerBeam::EulerBeamBCEnd end) const;

  bool is_free(EulerBeam::EulerBeamBCEnd end) const;
  bool is_simple(EulerBeam::EulerBeamBCEnd end) const;
  bool is_clamped(EulerBeam::EulerBeamBCEnd end) const;
  bool is_supported(EulerBeam::EulerBeamBCEnd end) const;
  bool is_supported_type(EulerBeam::EulerBeamBCType type) const;

  void validate_boundary_conditions() const;

  Vec3 endpoint_position(EulerBeam::EulerBeamBCEnd end) const;
  Vec3 endpoint_tangent(EulerBeam::EulerBeamBCEnd end) const;

  MatX3 tau_half(const MatX3& X) const;
  std::pair<real_t, real_t> compute_inextensibility_error(
    const MatX3& X) const;
  MatX3 dss_nodes(const MatX3& X) const;
  MatX3 bending_force(const MatX3& X) const;
  MatrixXd build_bending_matrix() const;
  VectorXd solve_tension(real_t dt, const Vec3& body_force) const;
  MatX3 solve_position(const VectorXd& T_half,
                       real_t dt,
                       const Vec3& body_force,
                       const MatX3& Fb_star) const;
  void update_mesh();
};

} // namespace Models
} // namespace ELFF
