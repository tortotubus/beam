#include "elff/models/beam/EulerBeamInextensibleADDMTestHarness.hpp"

#include <stdexcept>

namespace ELFF {
namespace Models {

void EulerBeamInextensibleADDMTestHarness::solve_linear(
    real_t dt, const std::vector<std::array<real_t, 3>> &averaged_load) {
  ELFF_ASSERT(averaged_load.size() == mesh.get_nodes(),
              "Size of load vector must equal number of nodes.");
  if (!(dt > 0.0)) {
    throw std::runtime_error("Linear Newmark: dt must be > 0");
  }

  have_prev_uniform_load = false;
  have_prev_nodal_load = false;
  nodal_load_prev.clear();

  x = x_prev;
  y = y_prev;
  z = z_prev;
  if (have_time_dependent_boundary_conditions) {
    apply_time_dependent_boundary_values(time_iter + 1);
    apply_prescribed_boundary_values(x, y, z);
  }

  const VectorXd x_old = x_prev;
  const VectorXd y_old = y_prev;
  const VectorXd z_old = z_prev;
  const real_t coeff = 2.0 / (dt * dt);

  A_unconstrained = 0.5 * EI * K_bending + coeff * mass_matrix;
  A_unconstrained.makeCompressed();
  A = A_unconstrained;
  apply_boundary_condition_A();
  decompose_A();

  clear_rhs();
  add_nodal_load_rhs(averaged_load);
  add_point_boundary_loads();
  apply_midpoint_bending_rhs();

  const VectorXd x_inertia = coeff * x_prev + (2.0 / dt) * vx_prev;
  const VectorXd y_inertia = coeff * y_prev + (2.0 / dt) * vy_prev;
  const VectorXd z_inertia = coeff * z_prev + (2.0 / dt) * vz_prev;

  f_x += mass_matrix * x_inertia;
  f_y += mass_matrix * y_inertia;
  f_z += mass_matrix * z_inertia;

  apply_boundary_condition_f();

  x = llt.solve(f_x);
  y = llt.solve(f_y);
  z = llt.solve(f_z);

  update_average_acceleration_state_component(x_old, x, vx_prev, ax_prev, dt);
  update_average_acceleration_state_component(y_old, y, vy_prev, ay_prev, dt);
  update_average_acceleration_state_component(z_old, z, vz_prev, az_prev, dt);
  apply_dynamic_state_boundary_conditions();

  x_prev = x;
  y_prev = y;
  z_prev = z;

  update_mesh();

  ++time_iter;
  t += dt;
}

} // namespace Models
} // namespace ELFF
