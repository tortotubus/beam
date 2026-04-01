#include "elff/models/beam/EulerBeamDynamicInextensibleAugKKT.hpp"

namespace ELFF {
namespace Models {

EulerBeamDynamicInextensibleAugKKT::
  EulerBeamDynamicInextensibleAugKKT(real_t length,
                                           real_t EI,
                                           real_t mu,
                                           size_t nodes,
                                           EulerBeam::EulerBeamBCs bcs,
                                           real_t r_penalty)
  : EulerBeamStaticInextensibleAugKKT(length, EI, mu, nodes, bcs, r_penalty)
  , v_prev(VectorXd::Zero(ndof))
  , a_prev(VectorXd::Zero(ndof))
  , u_prev(VectorXd::Zero(ndof))
  , mass(SparseMatrix<real_t>(ndof, ndof))
  , load_prev({ 0., 0., 0. })
{
  assemble_mass_matrix();
}

void
EulerBeamDynamicInextensibleAugKKT::solve(real_t dt,
                                                std::array<real_t, 3> load)
{
  const real_t alpha = 0.0;
  const real_t gamma = 0.5 - alpha;
  const real_t beta = 0.25 * (1 - alpha) * (1 - alpha);
  solve_newmark(dt, load, beta, gamma);
}

void
EulerBeamDynamicInextensibleAugKKT::solve(
  real_t dt,
  std::vector<std::array<real_t, 3>> load)
{
  ELFF_ASSERT(load.size() == nodes,
              "Size of load vector must equal number of nodes.");
  const real_t alpha = 0.0;
  const real_t gamma = 0.5 - alpha;
  const real_t beta = 0.25 * (1 - alpha) * (1 - alpha);
  solve_newmark(dt, load, beta, gamma);
}

void
EulerBeamDynamicInextensibleAugKKT::solve_newmark(
  real_t dt,
  std::array<real_t, 3> load,
  real_t beta,
  real_t gamma)
{
  SparseLU<SparseMatrix<real_t>, COLAMDOrdering<int>> solver;

  if (time_iter == 0) {
    EulerBeamStaticInextensibleAugKKT::assemble_system(load);
    apply_boundary_conditions();
    for (size_t n = 0; n < nodes; ++n) {
      const size_t ix = offset_x + 2 * n;
      const size_t iy = offset_y + 2 * n;
      const size_t iz = offset_z + 2 * n;
      const real_t w = (n == 0 || n == nodes - 1) ? 0.5 * ds : ds;
      a_prev(ix) = (-residual(ix)) / (mu * w);
      a_prev(iy) = (-residual(iy)) / (mu * w);
      a_prev(iz) = (-residual(iz)) / (mu * w);
    }
    apply_dynamic_state_boundary_conditions();
  }

  u_prev = u;
  size_t final_iter = 0;
  real_t final_res_norm = 0.0;

  for (size_t iter = 0; iter <= max_iter; ++iter) {
    assemble_system_newmark(dt, load, beta, gamma);
    apply_boundary_conditions();

    const real_t res_norm = residual.norm();
    final_iter = iter;
    final_res_norm = res_norm;

    if (res_norm < tol) {
      break;
    } else if (iter == max_iter) {
      break;
    }

    solver.analyzePattern(jacobian);
    solver.factorize(jacobian);

    if (solver.info() != Success) {
      ELFF_ABORT("EulerBeamDynamicInextensibleAugKKT::solve(): "
                 "SparseLU factorization failed.\n");
    }

    const VectorXd delta_u = solver.solve(-residual);

    if (solver.info() != Success) {
      ELFF_ABORT("EulerBeamDynamicInextensibleAugKKT::solve(): "
                 "SparseLU solve failed.\n");
    }

    u += delta_u;
  }

  ELFF_LOG(time_iter << "\t" << final_iter << "\t" << final_res_norm);

  for (size_t n = 0; n < nodes; ++n) {
    const size_t ix = offset_x + 2 * n;
    const size_t iy = offset_y + 2 * n;
    const size_t iz = offset_z + 2 * n;

    auto update_state = [&](size_t i) {
      const real_t a_new =
        (u(i) - u_prev(i) - dt * v_prev(i)) / (beta * dt * dt) -
        ((1.0 - 2.0 * beta) / (2.0 * beta)) * a_prev(i);
      const real_t v_new =
        v_prev(i) + dt * ((1.0 - gamma) * a_prev(i) + gamma * a_new);
      a_prev(i) = a_new;
      v_prev(i) = v_new;
    };

    update_state(ix);
    update_state(iy);
    update_state(iz);
  }

  apply_dynamic_state_boundary_conditions();
  u_prev = u;
  update_mesh();

  time_iter++;
  t += dt;
}

void
EulerBeamDynamicInextensibleAugKKT::solve_newmark(
  real_t dt,
  std::vector<std::array<real_t, 3>> load,
  real_t beta,
  real_t gamma)
{
  SparseLU<SparseMatrix<real_t>, COLAMDOrdering<int>> solver;

  if (time_iter == 0) {
    assemble_system(load);
    apply_boundary_conditions();
    for (size_t n = 0; n < nodes; ++n) {
      const size_t ix = offset_x + 2 * n;
      const size_t iy = offset_y + 2 * n;
      const size_t iz = offset_z + 2 * n;
      const real_t w = (n == 0 || n == nodes - 1) ? 0.5 * ds : ds;
      a_prev(ix) = (-residual(ix)) / (mu * w);
      a_prev(iy) = (-residual(iy)) / (mu * w);
      a_prev(iz) = (-residual(iz)) / (mu * w);
    }
    apply_dynamic_state_boundary_conditions();
  }

  u_prev = u;
  size_t final_iter = 0;
  real_t final_res_norm = 0.0;

  for (size_t iter = 0; iter <= max_iter; ++iter) {
    assemble_system_newmark(dt, load, beta, gamma);
    apply_boundary_conditions();

    const real_t res_norm = residual.norm();
    final_iter = iter;
    final_res_norm = res_norm;

    if (res_norm < tol) {
      break;
    } else if (iter == max_iter) {
      break;
    }

    solver.analyzePattern(jacobian);
    solver.factorize(jacobian);

    if (solver.info() != Success) {
      ELFF_ABORT("EulerBeamDynamicInextensibleAugKKT::solve(): "
                 "SparseLU factorization failed.\n");
    }

    const VectorXd delta_u = solver.solve(-residual);

    if (solver.info() != Success) {
      ELFF_ABORT("EulerBeamDynamicInextensibleAugKKT::solve(): "
                 "SparseLU solve failed.\n");
    }

    u += delta_u;
  }

  ELFF_LOG(time_iter << "\t" << final_iter << "\t" << final_res_norm);

  for (size_t n = 0; n < nodes; ++n) {
    const size_t ix = offset_x + 2 * n;
    const size_t iy = offset_y + 2 * n;
    const size_t iz = offset_z + 2 * n;

    auto update_state = [&](size_t i) {
      const real_t a_new =
        (u(i) - u_prev(i) - dt * v_prev(i)) / (beta * dt * dt) -
        ((1.0 - 2.0 * beta) / (2.0 * beta)) * a_prev(i);
      const real_t v_new =
        v_prev(i) + dt * ((1.0 - gamma) * a_prev(i) + gamma * a_new);
      a_prev(i) = a_new;
      v_prev(i) = v_new;
    };

    update_state(ix);
    update_state(iy);
    update_state(iz);
  }

  apply_dynamic_state_boundary_conditions();
  u_prev = u;
  update_mesh();

  time_iter++;
  t += dt;
}

void
EulerBeamDynamicInextensibleAugKKT::assemble_mass_matrix()
{
  using Tpl = Triplet<real_t>;

  std::vector<Tpl> triplets;
  triplets.reserve(elements * 3 * 16);

  const real_t xi_q[] = { 0.1127016654, 0.5, 0.8872983346 };
  const real_t w_q[] = { 0.2777777778, 0.4444444444, 0.2777777778 };

  for (size_t e = 0; e < elements; ++e) {
    const size_t n0 = e;
    const size_t n1 = e + 1;

    const std::array<size_t, 4> idx_x = { offset_x + 2 * n0 + 0,
                                          offset_x + 2 * n0 + 1,
                                          offset_x + 2 * n1 + 0,
                                          offset_x + 2 * n1 + 1 };
    const std::array<size_t, 4> idx_y = { offset_y + 2 * n0 + 0,
                                          offset_y + 2 * n0 + 1,
                                          offset_y + 2 * n1 + 0,
                                          offset_y + 2 * n1 + 1 };
    const std::array<size_t, 4> idx_z = { offset_z + 2 * n0 + 0,
                                          offset_z + 2 * n0 + 1,
                                          offset_z + 2 * n1 + 0,
                                          offset_z + 2 * n1 + 1 };

    real_t me[4][4] = { { 0.0 } };

    for (size_t qi = 0; qi < 3; ++qi) {
      const real_t xi = xi_q[qi];
      const real_t w = w_q[qi];
      const auto H = ELFF::FEM::CubicHermite<real_t>::values(xi, ds);

      for (size_t a = 0; a < 4; ++a) {
        for (size_t b = 0; b < 4; ++b) {
          me[a][b] += mu * H[a] * H[b] * w * ds;
        }
      }
    }

    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = 0; b < 4; ++b) {
        if (me[a][b] == 0.0) {
          continue;
        }
        triplets.emplace_back(idx_x[a], idx_x[b], me[a][b]);
        triplets.emplace_back(idx_y[a], idx_y[b], me[a][b]);
        triplets.emplace_back(idx_z[a], idx_z[b], me[a][b]);
      }
    }
  }

  mass.resize(ndof, ndof);
  mass.setFromTriplets(triplets.begin(), triplets.end());
  mass.makeCompressed();
}

void
EulerBeamDynamicInextensibleAugKKT::apply_initial_condition()
{
  EulerBeamStaticInextensibleAugKKT::apply_initial_condition();
  u_prev = u;
}

void
EulerBeamDynamicInextensibleAugKKT::apply_initial_condition(
  EulerBeamMesh& bmesh)
{
  EulerBeamStaticInextensibleAugKKT::apply_initial_condition(bmesh);
  u_prev = u;
}

void
EulerBeamDynamicInextensibleAugKKT::assemble_system_newmark(
  real_t dt,
  std::array<real_t, 3> load,
  real_t beta,
  real_t gamma)
{
  static_cast<void>(gamma);
  EulerBeamStaticInextensibleAugKKT::assemble_system(load);
  add_newmark_inertial_terms(dt, beta);
}

void
EulerBeamDynamicInextensibleAugKKT::assemble_system_newmark(
  real_t dt,
  std::vector<std::array<real_t, 3>> load,
  real_t beta,
  real_t gamma)
{
  static_cast<void>(gamma);
  assemble_system(load);
  add_newmark_inertial_terms(dt, beta);
}

void
EulerBeamDynamicInextensibleAugKKT::assemble_system(
  std::vector<std::array<real_t, 3>> load)
{
  using ADDeriv = Matrix<real_t, 14, 1>;
  using AD = AutoDiffScalar<ADDeriv>;
  using ADVec = Matrix<AD, 14, 1>;
  using Tpl = Triplet<real_t>;

  residual = VectorXd::Zero(ndof);
  std::vector<Tpl> triplets;
  triplets.reserve(elements * 14 * 14);

  for (size_t e = 0; e < elements; ++e) {
    const auto idx = get_element_dof_indices(e);
    const Matrix<real_t, 14, 1> u_elem = get_element_state(idx);
    const std::array<std::array<real_t, 3>, 2> load_elem = { load[e],
                                                             load[e + 1] };

    ADVec u_ad;
    for (int a = 0; a < 14; ++a) {
      ADDeriv seed = ADDeriv::Zero();
      seed(a) = 1.0;
      u_ad(a) = AD(u_elem(a), seed);
    }

    const ADVec r_elem_ad =
      assemble_element_residual_template<AD>(u_ad, load_elem);

    for (int a = 0; a < 14; ++a) {
      residual(idx[a]) += r_elem_ad(a).value();

      const ADDeriv& dRa = r_elem_ad(a).derivatives();
      for (int b = 0; b < 14; ++b) {
        const real_t value = dRa(b);
        if (value != 0.0) {
          triplets.emplace_back(idx[a], idx[b], value);
        }
      }
    }
  }

  jacobian.resize(ndof, ndof);
  jacobian.setFromTriplets(triplets.begin(), triplets.end());
  jacobian.makeCompressed();
}

void
EulerBeamDynamicInextensibleAugKKT::add_newmark_inertial_terms(real_t dt,
                                                                     real_t beta)
{
  if (!(dt > 0.0)) {
    throw std::runtime_error("Newmark: dt must be > 0");
  }
  if (!(beta > 0.0)) {
    throw std::runtime_error("Newmark: beta must be > 0");
  }

  const real_t inv = 1.0 / (beta * dt * dt);
  const real_t inv_bt = 1.0 / (beta * dt);
  const real_t kappa = (1.0 - 2.0 * beta) / (2.0 * beta);
  VectorXd a_new = VectorXd::Zero(ndof);

  for (size_t i = 0; i < offset_l; ++i) {
    a_new(i) = inv * (u(i) - u_prev(i)) - inv_bt * v_prev(i) - kappa * a_prev(i);
  }

  residual += mass * a_new;

  for (int col = 0; col < mass.outerSize(); ++col) {
    for (SparseMatrix<real_t>::InnerIterator it(mass, col); it; ++it) {
      jacobian.coeffRef(it.row(), it.col()) += inv * it.value();
    }
  }

  jacobian.makeCompressed();
}

void
EulerBeamDynamicInextensibleAugKKT::apply_dynamic_state_boundary_conditions()
{
  for (size_t bi = 0; bi < 2; ++bi) {
    const EulerBeamBCType bctype = boundary_conditions.type[bi];
    if (bctype != simple_bc && bctype != clamped_bc) {
      continue;
    }

    size_t ni = 0;
    switch (boundary_conditions.end[bi]) {
      case left:
        ni = 0;
        break;
      case right:
        ni = nodes - 1;
        break;
    }

    const size_t ix = offset_x + 2 * ni;
    const size_t iy = offset_y + 2 * ni;
    const size_t iz = offset_z + 2 * ni;
    const size_t sx = offset_x + 2 * ni + 1;
    const size_t sy = offset_y + 2 * ni + 1;
    const size_t sz = offset_z + 2 * ni + 1;

    v_prev(ix) = 0.0;
    v_prev(iy) = 0.0;
    v_prev(iz) = 0.0;
    a_prev(ix) = 0.0;
    a_prev(iy) = 0.0;
    a_prev(iz) = 0.0;

    if (bctype == clamped_bc) {
      v_prev(sx) = 0.0;
      v_prev(sy) = 0.0;
      v_prev(sz) = 0.0;
      a_prev(sx) = 0.0;
      a_prev(sy) = 0.0;
      a_prev(sz) = 0.0;
    }
  }
}

void
EulerBeamDynamicInextensibleAugKKT::update_mesh()
{
  std::vector<std::array<real_t, 3>>& centerline = mesh.get_centerline();
  std::vector<std::array<real_t, 3>>& velocity = mesh.get_centerline_velocity();
  std::vector<std::array<real_t, 3>>& slope = mesh.get_slope();

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

} // namespace Models
} // namespace ELFF
