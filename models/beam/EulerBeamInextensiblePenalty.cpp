#include "elff/general/error.hpp"
#include "elff/models/beam/EulerBeamInextensiblePenalty.hpp"

namespace ELFF {
namespace Models {

EulerBeamInextensiblePenalty::EulerBeamInextensiblePenalty(
  real_t                  length,
  real_t                  EI,
  size_t                  nodes,
  EulerBeam::EulerBeamBCs bcs,
  real_t                  r_penalty)
  : EulerBeam(length, EI, nodes, bcs)
  , dimension(3)
  , elements(nodes - 1)
  , nodes(nodes)
  , ds(mesh.get_ds())
  , ndof_x(2 * nodes)
  , ndof_y(2 * nodes)
  , ndof_z(2 * nodes)
  , offset_x(0)
  , offset_y(ndof_x)
  , offset_z(ndof_x + ndof_y)
  , ndof(ndof_x + ndof_y + ndof_z)
  , r_penalty(r_penalty)
  , max_iter_nonlinear(1000)
  , tol_primal(1e-5)
  , residual(VectorXd::Zero(ndof))
  , jacobian(SparseMatrix<real_t>(ndof, ndof))
  , u(VectorXd::Zero(ndof))
  , u_prev(VectorXd::Zero(ndof))
  , v_prev(VectorXd::Zero(ndof))
  , a_prev(VectorXd::Zero(ndof))
{
  apply_initial_condition(mesh);
}

EulerBeamInextensiblePenalty::EulerBeamInextensiblePenalty(
  real_t                  length,
  real_t                  EI,
  real_t                  mu,
  size_t                  nodes,
  EulerBeam::EulerBeamBCs bcs,
  real_t                  r_penalty)
  : EulerBeam(length, EI, mu, nodes, bcs)
  , dimension(3)
  , elements(nodes - 1)
  , nodes(nodes)
  , ds(mesh.get_ds())
  , ndof_x(2 * nodes)
  , ndof_y(2 * nodes)
  , ndof_z(2 * nodes)
  , offset_x(0)
  , offset_y(ndof_x)
  , offset_z(ndof_x + ndof_y)
  , ndof(ndof_x + ndof_y + ndof_z)
  , r_penalty(r_penalty)
  , max_iter_nonlinear(1000)
  , tol_primal(1e-5)
  , residual(VectorXd::Zero(ndof))
  , jacobian(SparseMatrix<real_t>(ndof, ndof))
  , u(VectorXd::Zero(ndof))
  , u_prev(VectorXd::Zero(ndof))
  , v_prev(VectorXd::Zero(ndof))
  , a_prev(VectorXd::Zero(ndof))
{
  apply_initial_condition(mesh);
}

void
EulerBeamInextensiblePenalty::solve()
{
  solve({ 0.0, 0.0, 0.0 });
}

void
EulerBeamInextensiblePenalty::solve(std::array<real_t, 3> load)
{
  SparseLU<SparseMatrix<real_t>, COLAMDOrdering<int>> solver;

  for (size_t iter = 0; iter < max_iter_nonlinear; ++iter) {
    assemble_system(load);
    apply_boundary_conditions();

    const real_t res_norm = residual.norm();
    if (res_norm < tol_primal) {
      ELFF_LOG(iter << ": ||r|| = " << res_norm);
      break;
    }

    if (iter == max_iter_nonlinear - 1) {
      ELFF_ABORT("EulerBeamInextensiblePenalty::solve() did not converge.\n");
    }

    ELFF_LOG(iter << ": ||r|| = " << res_norm);

    solver.compute(jacobian);
    if (solver.info() != Success) {
      ELFF_ABORT("EulerBeamInextensiblePenalty::solve(): "
                 "SparseLU factorization failed.\n");
    }

    const VectorXd delta_u = solver.solve(-residual);
    if (solver.info() != Success) {
      ELFF_ABORT("EulerBeamInextensiblePenalty::solve(): "
                 "SparseLU solve failed.\n");
    }

    u += delta_u;
  }

  update_mesh();
}

void
EulerBeamInextensiblePenalty::solve(std::vector<std::array<real_t, 3>> load)
{
  ELFF_ASSERT(load.size() == nodes,
              "Size of load vector must equal number of nodes.");

  SparseLU<SparseMatrix<real_t>, COLAMDOrdering<int>> solver;

  for (size_t iter = 0; iter < max_iter_nonlinear; ++iter) {
    assemble_system(load);
    apply_boundary_conditions();

    const real_t res_norm = residual.norm();
    if (res_norm < tol_primal) {
      ELFF_LOG(iter << ": ||r|| = " << res_norm);
      break;
    }

    if (iter == max_iter_nonlinear - 1) {
      ELFF_ABORT("EulerBeamInextensiblePenalty::solve() did not converge.\n");
    }

    ELFF_LOG(iter << ": ||r|| = " << res_norm);

    solver.compute(jacobian);
    if (solver.info() != Success) {
      ELFF_ABORT("EulerBeamInextensiblePenalty::solve(): "
                 "SparseLU factorization failed.\n");
    }

    const VectorXd delta_u = solver.solve(-residual);
    if (solver.info() != Success) {
      ELFF_ABORT("EulerBeamInextensiblePenalty::solve(): "
                 "SparseLU solve failed.\n");
    }

    u += delta_u;
  }

  update_mesh();
}

void
EulerBeamInextensiblePenalty::solve(real_t dt, std::array<real_t, 3> load)
{
  ELFF_ASSERT(mu > 0.0,
              "Dynamic EulerBeamInextensiblePenalty requires mu > 0.\n");

  const real_t alpha = 0.0;
  const real_t gamma = 0.5 - alpha;
  const real_t beta = 0.25 * (1.0 - alpha) * (1.0 - alpha);
  solve_newmark(dt, load, beta, gamma);
}

void
EulerBeamInextensiblePenalty::solve(real_t                              dt,
                                    std::vector<std::array<real_t, 3>> load)
{
  ELFF_ASSERT(mu > 0.0,
              "Dynamic EulerBeamInextensiblePenalty requires mu > 0.\n");
  ELFF_ASSERT(load.size() == nodes,
              "Size of load vector must equal number of nodes.");

  const real_t alpha = 0.0;
  const real_t gamma = 0.5 - alpha;
  const real_t beta = 0.25 * (1.0 - alpha) * (1.0 - alpha);
  solve_newmark(dt, load, beta, gamma);
}

void
EulerBeamInextensiblePenalty::apply_initial_condition()
{
  apply_initial_condition(mesh);
}

void
EulerBeamInextensiblePenalty::apply_initial_condition(EulerBeamMesh& bmesh)
{
  ELFF_ASSERT(nodes == bmesh.get_nodes(),
              "Provided mesh must have same node count as the beam mesh.\n");

  const auto centerline = bmesh.get_centerline();
  const auto slopes = bmesh.get_slope();
  const auto velocity = bmesh.get_centerline_velocity();

  for (size_t i = 0; i < nodes; ++i) {
    u(offset_x + 2 * i + 0) = centerline[i][0];
    u(offset_x + 2 * i + 1) = slopes[i][0];
    u(offset_y + 2 * i + 0) = centerline[i][1];
    u(offset_y + 2 * i + 1) = slopes[i][1];
    u(offset_z + 2 * i + 0) = centerline[i][2];
    u(offset_z + 2 * i + 1) = slopes[i][2];

    v_prev(offset_x + 2 * i + 0) = velocity[i][0];
    v_prev(offset_y + 2 * i + 0) = velocity[i][1];
    v_prev(offset_z + 2 * i + 0) = velocity[i][2];
  }

  u_prev = u;
  a_prev.setZero();
  update_mesh();
}

std::array<size_t, 12>
EulerBeamInextensiblePenalty::get_element_dof_indices(size_t e) const
{
  const size_t n0 = e;
  const size_t n1 = e + 1;

  return {
    offset_x + 2 * n0 + 0, offset_x + 2 * n0 + 1, offset_x + 2 * n1 + 0,
    offset_x + 2 * n1 + 1, offset_y + 2 * n0 + 0, offset_y + 2 * n0 + 1,
    offset_y + 2 * n1 + 0, offset_y + 2 * n1 + 1, offset_z + 2 * n0 + 0,
    offset_z + 2 * n0 + 1, offset_z + 2 * n1 + 0, offset_z + 2 * n1 + 1
  };
}

Matrix<real_t, 12, 1>
EulerBeamInextensiblePenalty::get_element_state(
  const std::array<size_t, 12>& idx) const
{
  Matrix<real_t, 12, 1> u_elem;
  for (int i = 0; i < 12; ++i) {
    u_elem(i) = u(idx[i]);
  }
  return u_elem;
}

void
EulerBeamInextensiblePenalty::assemble_system(std::array<real_t, 3> load)
{
  using ADDeriv = Matrix<real_t, 12, 1>;
  using AD = AutoDiffScalar<ADDeriv>;
  using ADVec = Matrix<AD, 12, 1>;
  using Tpl = Triplet<real_t>;

  residual.setZero();
  std::vector<Tpl> triplets;
  triplets.reserve(elements * 12 * 12);

  for (size_t e = 0; e < elements; ++e) {
    const auto idx = get_element_dof_indices(e);
    const Matrix<real_t, 12, 1> u_elem = get_element_state(idx);

    ADVec u_ad;
    for (int a = 0; a < 12; ++a) {
      ADDeriv seed = ADDeriv::Zero();
      seed(a) = 1.0;
      u_ad(a) = AD(u_elem(a), seed);
    }

    const ADVec R_loc_ad = assemble_element_residual_template<AD>(u_ad, load);

    for (int a = 0; a < 12; ++a) {
      residual(idx[a]) += R_loc_ad(a).value();

      const ADDeriv& dRa = R_loc_ad(a).derivatives();
      for (int b = 0; b < 12; ++b) {
        const real_t dj = dRa(b);
        if (dj != 0.0) {
          triplets.emplace_back(idx[a], idx[b], dj);
        }
      }
    }
  }

  jacobian.resize(ndof, ndof);
  jacobian.setFromTriplets(triplets.begin(), triplets.end());
  jacobian.makeCompressed();
}

void
EulerBeamInextensiblePenalty::assemble_system(
  std::vector<std::array<real_t, 3>> load)
{
  using ADDeriv = Matrix<real_t, 12, 1>;
  using AD = AutoDiffScalar<ADDeriv>;
  using ADVec = Matrix<AD, 12, 1>;
  using Tpl = Triplet<real_t>;

  residual.setZero();
  std::vector<Tpl> triplets;
  triplets.reserve(elements * 12 * 12);

  for (size_t e = 0; e < elements; ++e) {
    const auto idx = get_element_dof_indices(e);
    const Matrix<real_t, 12, 1> u_elem = get_element_state(idx);
    const std::array<std::array<real_t, 3>, 2> load_elem = { load[e],
                                                             load[e + 1] };

    ADVec u_ad;
    for (int a = 0; a < 12; ++a) {
      ADDeriv seed = ADDeriv::Zero();
      seed(a) = 1.0;
      u_ad(a) = AD(u_elem(a), seed);
    }

    const ADVec R_loc_ad =
      assemble_element_residual_template<AD>(u_ad, load_elem);

    for (int a = 0; a < 12; ++a) {
      residual(idx[a]) += R_loc_ad(a).value();

      const ADDeriv& dRa = R_loc_ad(a).derivatives();
      for (int b = 0; b < 12; ++b) {
        const real_t dj = dRa(b);
        if (dj != 0.0) {
          triplets.emplace_back(idx[a], idx[b], dj);
        }
      }
    }
  }

  jacobian.resize(ndof, ndof);
  jacobian.setFromTriplets(triplets.begin(), triplets.end());
  jacobian.makeCompressed();
}

std::pair<real_t, real_t>
EulerBeamInextensiblePenalty::compute_inextensibility_error(
  const VectorXd& state) const
{
  const real_t xi_q[] = { 0.1127016654, 0.5, 0.8872983346 };

  real_t max_abs = 0.0;
  real_t sum_sq = 0.0;
  size_t sample_count = 0;

  for (size_t e = 0; e < elements; ++e) {
    const auto idx = get_element_dof_indices(e);

    Matrix<real_t, 12, 1> u_elem;
    for (int i = 0; i < 12; ++i) {
      u_elem(i) = state(idx[i]);
    }

    const Matrix<real_t, 4, 1> ux = u_elem.segment<4>(0);
    const Matrix<real_t, 4, 1> uy = u_elem.segment<4>(4);
    const Matrix<real_t, 4, 1> uz = u_elem.segment<4>(8);

    for (real_t xi : xi_q) {
      const auto dH = ELFF::FEM::CubicHermite<real_t>::derivs(xi, ds);

      real_t xp = 0.0;
      real_t yp = 0.0;
      real_t zp = 0.0;
      for (size_t i = 0; i < 4; ++i) {
        xp += dH[i] * ux(i);
        yp += dH[i] * uy(i);
        zp += dH[i] * uz(i);
      }

      const real_t g = xp * xp + yp * yp + zp * zp - 1.0;
      const real_t abs_g = std::abs(g);
      max_abs = std::max(max_abs, abs_g);
      sum_sq += g * g;
      ++sample_count;
    }
  }

  if (sample_count == 0) {
    return { 0.0, 0.0 };
  }

  return { std::sqrt(sum_sq / static_cast<real_t>(sample_count)), max_abs };
}

void
EulerBeamInextensiblePenalty::apply_boundary_conditions()
{
  for (size_t bi = 0; bi < 2; ++bi) {
    const EulerBeamBCEnd bcend = boundary_conditions.end[bi];
    size_t ni = 0;
    std::vector<size_t> idx(6);

    switch (bcend) {
      case left:
        ni = 0;
        break;
      case right:
        ni = nodes - 1;
        break;
    }

    const EulerBeamBCType bctype = boundary_conditions.type[bi];
    const EulerBeamBCVals bcvals = boundary_conditions.vals[bi];
    std::vector<real_t> vals(6);

    switch (bctype) {
      case free_bc:
        idx = {};
        vals = {};
        break;
      case simple_bc:
        idx = { offset_x + 2 * ni + 0,
                offset_y + 2 * ni + 0,
                offset_z + 2 * ni + 0 };
        vals = { bcvals.position[0], bcvals.position[1], bcvals.position[2] };
        break;
      case clamped_bc:
        idx = { offset_x + 2 * ni + 0, offset_x + 2 * ni + 1,
                offset_y + 2 * ni + 0, offset_y + 2 * ni + 1,
                offset_z + 2 * ni + 0, offset_z + 2 * ni + 1 };
        vals = { bcvals.position[0], bcvals.slope[0],    bcvals.position[1],
                 bcvals.slope[1],    bcvals.position[2], bcvals.slope[2] };
        break;
      case point_force_bc:
        idx = { offset_x + 2 * ni + 0,
                offset_y + 2 * ni + 0,
                offset_z + 2 * ni + 0 };
        vals = { bcvals.force[0], bcvals.force[1], bcvals.force[2] };
        break;
      case point_torque_bc:
        idx = { offset_x + 2 * ni + 1,
                offset_y + 2 * ni + 1,
                offset_z + 2 * ni + 1 };
        vals = { bcvals.torque[0], bcvals.torque[1], bcvals.torque[2] };
        break;
      default:
        break;
    }

    switch (bctype) {
      case point_force_bc:
      case point_torque_bc:
        for (size_t i = 0; i < vals.size(); ++i) {
          residual(idx[i]) -= vals[i];
        }
        break;
      default:
        for (auto i : idx) {
          for (SparseMatrix<real_t>::InnerIterator it(jacobian, i); it; ++it) {
            it.valueRef() = 0.0;
          }

          for (int col = 0; col < jacobian.outerSize(); ++col) {
            for (SparseMatrix<real_t>::InnerIterator it(jacobian, col); it;
                 ++it) {
              if (it.row() == static_cast<int>(i)) {
                it.valueRef() = 0.0;
              }
            }
          }

          jacobian.coeffRef(i, i) = 1.0;
        }

        for (size_t i = 0; i < vals.size(); ++i) {
          residual(idx[i]) = u(idx[i]) - vals[i];
        }
        break;
    }
  }
}

void
EulerBeamInextensiblePenalty::update_mesh()
{
  auto& centerline = mesh.get_centerline();
  auto& velocity = mesh.get_centerline_velocity();
  auto& slope = mesh.get_slope();

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

void
EulerBeamInextensiblePenalty::solve_newmark(real_t dt,
                                            std::array<real_t, 3> load,
                                            real_t beta,
                                            real_t gamma)
{
  SparseLU<SparseMatrix<real_t>, COLAMDOrdering<int>> solver;

  if (time_iter == 0) {
    initialize_acceleration(load);
  }

  u_prev = u;
  const VectorXd v_old = v_prev;
  const VectorXd a_old = a_prev;
  size_t final_iter = 0;
  real_t final_res_norm = 0.0;

  for (size_t iter = 0; iter < max_iter_nonlinear; ++iter) {
    assemble_system(load);
    add_newmark_inertial_terms(dt, beta);
    apply_boundary_conditions();

    const real_t res_norm = residual.norm();
    final_iter = iter;
    final_res_norm = res_norm;
    if (res_norm < tol_primal) {
      break;
    }

    if (iter == max_iter_nonlinear - 1) {
      ELFF_ABORT(
        "EulerBeamInextensiblePenalty::solve_newmark() did not converge.\n");
    }

    solver.compute(jacobian);
    if (solver.info() != Success) {
      ELFF_ABORT("EulerBeamInextensiblePenalty::solve_newmark(): "
                 "SparseLU factorization failed.\n");
    }

    const VectorXd delta_u = solver.solve(-residual);
    if (solver.info() != Success) {
      ELFF_ABORT("EulerBeamInextensiblePenalty::solve_newmark(): "
                 "SparseLU solve failed.\n");
    }

    u += delta_u;
  }

  const auto [inext_l2, inext_lmax] = compute_inextensibility_error(u);
  ELFF_LOG(time_iter << "\t" << final_iter << "\t" << final_res_norm
                     << "\t|g|_l2=" << inext_l2
                     << "\t|g|_max=" << inext_lmax);

  update_newmark_state_from_displacement(dt, beta, gamma, v_old, a_old);
  apply_dynamic_state_boundary_conditions();
  u_prev = u;
  update_mesh();
  ++time_iter;
  t += dt;
}

void
EulerBeamInextensiblePenalty::solve_newmark(
  real_t                              dt,
  std::vector<std::array<real_t, 3>> load,
  real_t                              beta,
  real_t                              gamma)
{
  ELFF_ASSERT(load.size() == nodes,
              "Size of load vector must equal number of nodes.");

  SparseLU<SparseMatrix<real_t>, COLAMDOrdering<int>> solver;

  if (time_iter == 0) {
    initialize_acceleration(load);
  }

  u_prev = u;
  const VectorXd v_old = v_prev;
  const VectorXd a_old = a_prev;
  size_t final_iter = 0;
  real_t final_res_norm = 0.0;

  for (size_t iter = 0; iter < max_iter_nonlinear; ++iter) {
    assemble_system(load);
    add_newmark_inertial_terms(dt, beta);
    apply_boundary_conditions();

    const real_t res_norm = residual.norm();
    final_iter = iter;
    final_res_norm = res_norm;
    if (res_norm < tol_primal) {
      break;
    }

    if (iter == max_iter_nonlinear - 1) {
      ELFF_ABORT(
        "EulerBeamInextensiblePenalty::solve_newmark() did not converge.\n");
    }

    solver.compute(jacobian);
    if (solver.info() != Success) {
      ELFF_ABORT("EulerBeamInextensiblePenalty::solve_newmark(): "
                 "SparseLU factorization failed.\n");
    }

    const VectorXd delta_u = solver.solve(-residual);
    if (solver.info() != Success) {
      ELFF_ABORT("EulerBeamInextensiblePenalty::solve_newmark(): "
                 "SparseLU solve failed.\n");
    }

    u += delta_u;
  }

  const auto [inext_l2, inext_lmax] = compute_inextensibility_error(u);
  ELFF_LOG(time_iter << "\t" << final_iter << "\t" << final_res_norm
                     << "\t|g|_l2=" << inext_l2
                     << "\t|g|_max=" << inext_lmax);

  update_newmark_state_from_displacement(dt, beta, gamma, v_old, a_old);
  apply_dynamic_state_boundary_conditions();
  u_prev = u;
  update_mesh();
  ++time_iter;
  t += dt;
}

void
EulerBeamInextensiblePenalty::initialize_acceleration(std::array<real_t, 3> load)
{
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
}

void
EulerBeamInextensiblePenalty::initialize_acceleration(
  std::vector<std::array<real_t, 3>> load)
{
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
}

void
EulerBeamInextensiblePenalty::add_newmark_inertial_terms(real_t dt, real_t beta)
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

  for (size_t n = 0; n < nodes; ++n) {
    const size_t ix = offset_x + 2 * n;
    const size_t iy = offset_y + 2 * n;
    const size_t iz = offset_z + 2 * n;
    const real_t w = (n == 0 || n == nodes - 1) ? 0.5 * ds : ds;
    const real_t mass_coeff = mu * w;
    const real_t tangent_coeff = mass_coeff * inv;

    auto add_newmark = [&](size_t i) {
      const real_t a_new =
        inv * (u(i) - u_prev(i)) - inv_bt * v_prev(i) - kappa * a_prev(i);
      residual(i) += mass_coeff * a_new;
      jacobian.coeffRef(i, i) += tangent_coeff;
    };

    add_newmark(ix);
    add_newmark(iy);
    add_newmark(iz);
  }

  jacobian.makeCompressed();
}

void
EulerBeamInextensiblePenalty::update_newmark_state_from_displacement(
  real_t          dt,
  real_t          beta,
  real_t          gamma,
  const VectorXd& v_old,
  const VectorXd& a_old)
{
  if (!(dt > 0.0)) {
    throw std::runtime_error("Newmark: dt must be > 0");
  }
  if (!(beta > 0.0)) {
    throw std::runtime_error("Newmark: beta must be > 0");
  }
  if (!(gamma > 0.0)) {
    throw std::runtime_error("Newmark: gamma must be > 0");
  }

  const real_t inv = 1.0 / (beta * dt * dt);
  const real_t kappa = (1.0 - 2.0 * beta) / (2.0 * beta);

  for (Eigen::Index i = 0; i < u.size(); ++i) {
    const real_t a_new =
      inv * (u(i) - u_prev(i) - dt * v_old(i)) - kappa * a_old(i);
    const real_t v_new =
      v_old(i) + dt * ((1.0 - gamma) * a_old(i) + gamma * a_new);
    a_prev(i) = a_new;
    v_prev(i) = v_new;
  }
}

void
EulerBeamInextensiblePenalty::apply_dynamic_state_boundary_conditions()
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

} // namespace Models
} // namespace ELFF
