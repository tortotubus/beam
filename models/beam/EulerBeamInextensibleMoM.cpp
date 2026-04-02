#include "elff/models/beam/EulerBeamInextensibleMoM.hpp"

#include <cmath>
#include <stdexcept>

namespace ELFF {
namespace Models {

EulerBeamInextensibleMoM::EulerBeamInextensibleMoM(
  real_t length,
  real_t EI,
  size_t nodes,
  EulerBeam::EulerBeamBCs bcs,
  real_t r_penalty)
  : EulerBeam(length, EI, nodes, bcs)
  , dimension()
  , elements(nodes - 1)
  , nodes(nodes)
  , ds(mesh.get_ds())
  , ndof_x(2 * nodes)
  , ndof_y(2 * nodes)
  , ndof_z(2 * nodes)
  , ndof_l(1 * nodes)
  , offset_x(0)
  , offset_y(ndof_x)
  , offset_z(ndof_x + ndof_y)
  , offset_l()
  , ndof(ndof_x + ndof_y + ndof_z)
  , r_penalty(r_penalty)
  , max_iter_inner(100)
  , min_iter_inner(1)
  , max_iter_outer(1000)
  , min_iter_outer(1)
  , tol_linear(5e-5)
  , tol_primal(1e-5)
  , tol_constraint(5e-5)
  , residual(VectorXd::Zero(ndof))
  , lambda(VectorXd::Zero(ndof_l))
  , jacobian(SparseMatrix<real_t>(ndof, ndof))
  , mass(SparseMatrix<real_t>(ndof, ndof))
  , bending_matrix(SparseMatrix<real_t>(ndof, ndof))
  , u(VectorXd::Zero(ndof))
  , u_prev(VectorXd::Zero(ndof))
  , v_prev(VectorXd::Zero(ndof))
  , a_prev(VectorXd::Zero(ndof))
  , load_prev({ 0.0, 0.0, 0.0 })
  , nodal_load_prev(nodes, { 0.0, 0.0, 0.0 })
  , have_prev_uniform_load(false)
  , have_prev_nodal_load(false)
{
  assemble_mass_matrix();
  assemble_bending_matrix();
  apply_initial_condition(mesh);
}

EulerBeamInextensibleMoM::EulerBeamInextensibleMoM(
  real_t length,
  real_t EI,
  real_t mu,
  size_t nodes,
  EulerBeam::EulerBeamBCs bcs,
  real_t r_penalty)
  : EulerBeam(length, EI, mu, nodes, bcs)
  , dimension()
  , elements(nodes - 1)
  , nodes(nodes)
  , ds(mesh.get_ds())
  , ndof_x(2 * nodes)
  , ndof_y(2 * nodes)
  , ndof_z(2 * nodes)
  , ndof_l(1 * nodes)
  , offset_x(0)
  , offset_y(ndof_x)
  , offset_z(ndof_x + ndof_y)
  , offset_l(ndof_x + ndof_y + ndof_z)
  , ndof(ndof_x + ndof_y + ndof_z)
  , r_penalty(r_penalty)
  , max_iter_inner(2)
  , min_iter_inner(1)
  , max_iter_outer(1000)
  , min_iter_outer(1)
  , tol_linear(5e-12)
  , tol_primal(1e-5)
  , tol_constraint(100)
  , residual(VectorXd::Zero(ndof))
  , lambda(VectorXd::Zero(ndof_l))
  , jacobian(SparseMatrix<real_t>(ndof, ndof))
  , mass(SparseMatrix<real_t>(ndof, ndof))
  , bending_matrix(SparseMatrix<real_t>(ndof, ndof))
  , u(VectorXd::Zero(ndof))
  , u_prev(VectorXd::Zero(ndof))
  , v_prev(VectorXd::Zero(ndof))
  , a_prev(VectorXd::Zero(ndof))
  , load_prev({ 0.0, 0.0, 0.0 })
  , nodal_load_prev(nodes, { 0.0, 0.0, 0.0 })
  , have_prev_uniform_load(false)
  , have_prev_nodal_load(false)
{
  assemble_mass_matrix();
  assemble_bending_matrix();
  apply_initial_condition(mesh);
}

EulerBeamInextensibleMoM::~EulerBeamInextensibleMoM() =
  default;

void
EulerBeamInextensibleMoM::solve()
{
  solve({ 0., 0., 0. });
}

void
EulerBeamInextensibleMoM::solve(std::array<real_t, 3> load)
{
  real_t dlambda_norm = 0, res_norm = 0, constr_norm = 0.;
  size_t iter_outer = 0, iter_inner = 0;

  SparseLU<SparseMatrix<real_t>, COLAMDOrdering<int>> solver;

  for (iter_outer = 0; iter_outer < max_iter_outer; iter_outer++) {
    for (iter_inner = 0; iter_inner < max_iter_inner; iter_inner++) {
      assemble_system(load);
      apply_boundary_conditions();
      res_norm = residual.norm();

      if (iter_inner + 1 >= min_iter_inner)
        if (res_norm < tol_primal)
          break;
      solver.compute(jacobian);
      if (solver.info() != Success) {
        ELFF_ABORT("EulerBeamInextensibleMoM::solve(): "
                   "SparseLU factorization failed.\n");
      }
      const VectorXd delta_u = solver.solve(-residual);
      if (solver.info() != Success) {
        ELFF_ABORT("EulerBeamInextensibleMoM::solve(): "
                   "SparseLU solve failed.\n");
      }
      u += delta_u;
    }
    // Update lambda
    dlambda_norm = update_lambda();
    constr_norm = compute_inextensibility_error_l2();

    ELFF_LOG(iter_outer << "\t" << iter_inner << "\t" << res_norm << "\t"
                        << dlambda_norm << "\t" << constr_norm);

    if (iter_outer + 1 >= min_iter_outer)
      if (dlambda_norm < tol_constraint)
        break;
  }

  if (iter_outer == max_iter_outer) {
    ELFF_LOG(iter_outer << "\t" << iter_inner << "\t" << res_norm << "\t"
                        << dlambda_norm << "\t" << constr_norm);
    ELFF_ABORT("EulerBeamInextensibleMoM::solve(): "
               "Maximum outer iterations reached.\n");
  }

  update_mesh();
}

void
EulerBeamInextensibleMoM::solve(
  std::vector<std::array<real_t, 3>> load)
{
  ELFF_ASSERT(load.size() == nodes,
              "Size of load vector must equal number of nodes.");

  real_t dlambda_norm = 0, res_norm = 0, constr_norm = 0.;
  size_t iter_outer = 0, iter_inner = 0;

  ConjugateGradient<SparseMatrix<real_t>,
                    Lower | Upper,
                    IncompleteCholesky<real_t>>
    solver;

  for (iter_outer = 0; iter_outer < max_iter_outer; iter_outer++) {
    for (iter_inner = 0; iter_inner < max_iter_inner; iter_inner++) {
      assemble_system(load);
      apply_boundary_conditions();
      res_norm = residual.norm();

      if (iter_inner + 1 >= min_iter_inner)
        if (res_norm < tol_primal)
          break;
      solver.setTolerance(tol_linear);
      solver.compute(jacobian);

      if (solver.info() != Success) {
        ELFF_ABORT("EulerBeamInextensibleMoM::solve(): "
                   "Preconditioner failed.\n");
      }

      const VectorXd delta_u = solver.solve(-residual);
      u += delta_u;
    }

    dlambda_norm = update_lambda();
    constr_norm = compute_inextensibility_error_l2();

    ELFF_LOG(iter_outer << "\t" << iter_inner << "\t" << res_norm << "\t"
                        << dlambda_norm << "\t" << constr_norm);

    if (iter_outer + 1 >= min_iter_outer)
      if (dlambda_norm < tol_constraint)
        break;
  }

  if (iter_outer == max_iter_outer) {
    ELFF_LOG(iter_outer << "\t" << iter_inner << "\t" << res_norm << "\t"
                        << dlambda_norm << "\t" << constr_norm);
    ELFF_ABORT("EulerBeamInextensibleMoM::solve(): "
               "Maximum outer iterations reached.\n");
  }

  update_mesh();
}

void
EulerBeamInextensibleMoM::solve(real_t dt,
                                      std::array<real_t, 3> load)
{
  ELFF_ASSERT(mu > 0.0,
              "Dynamic EulerBeamInextensibleMoM requires mu > 0.\n");

  const real_t alpha = 0.0;
  const real_t gamma = 0.5 - alpha;
  const real_t beta = 0.25 * (1.0 - alpha) * (1.0 - alpha);
  solve_newmark(dt, load, beta, gamma);
}

void
EulerBeamInextensibleMoM::solve(
  real_t dt,
  std::vector<std::array<real_t, 3>> load)
{
  ELFF_ASSERT(mu > 0.0,
              "Dynamic EulerBeamInextensibleMoM requires mu > 0.\n");
  ELFF_ASSERT(load.size() == nodes,
              "Size of load vector must equal number of nodes.");

  const real_t alpha = 0.0;
  const real_t gamma = 0.5 - alpha;
  const real_t beta = 0.25 * (1.0 - alpha) * (1.0 - alpha);
  solve_newmark(dt, load, beta, gamma);
}

void
EulerBeamInextensibleMoM::apply_initial_condition(
  EulerBeamMesh& bmesh)
{
  ELFF_ASSERT(
    nodes == bmesh.get_nodes(),
    "Provided mesh must have same node count as the previous mesh.\n");

  const auto centerline = bmesh.get_centerline();
  const auto slopes = bmesh.get_slope();
  const auto velocity = bmesh.get_centerline_velocity();

  for (size_t ni = 0; ni < nodes; ni++) {
    u(offset_x + 2 * ni + 0) = centerline[ni][0];
    u(offset_x + 2 * ni + 1) = slopes[ni][0];
    u(offset_y + 2 * ni + 0) = centerline[ni][1];
    u(offset_y + 2 * ni + 1) = slopes[ni][1];
    u(offset_z + 2 * ni + 0) = centerline[ni][2];
    u(offset_z + 2 * ni + 1) = slopes[ni][2];

    v_prev(offset_x + 2 * ni + 0) = velocity[ni][0];
    v_prev(offset_y + 2 * ni + 0) = velocity[ni][1];
    v_prev(offset_z + 2 * ni + 0) = velocity[ni][2];
  }

  u_prev = u;
  a_prev.setZero();
  have_prev_uniform_load = false;
  have_prev_nodal_load = false;
  update_mesh();
}

void
EulerBeamInextensibleMoM::apply_initial_condition()
{
  for (size_t i = 0; i < nodes; i++) {
    u(offset_x + 2 * i + 0) = ds * i;
    u(offset_x + 2 * i + 1) = 1.;
    u(offset_y + 2 * i + 0) = 0.;
    u(offset_y + 2 * i + 1) = 0.;
    u(offset_z + 2 * i + 0) = 0.;
    u(offset_z + 2 * i + 1) = 0.;
  }

  u_prev = u;
  v_prev.setZero();
  a_prev.setZero();
  have_prev_uniform_load = false;
  have_prev_nodal_load = false;
  update_mesh();
}

void
EulerBeamInextensibleMoM::assemble_residual(
  std::array<real_t, 3> load)
{
  residual = assemble_residual_template<real_t>(u, load);
}

void
EulerBeamInextensibleMoM::assemble_residual(
  std::vector<std::array<real_t, 3>> load)
{
  residual = assemble_residual_template<real_t>(u, load);
}

std::array<size_t, 12>
EulerBeamInextensibleMoM::get_element_dof_indices(size_t e) const
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
EulerBeamInextensibleMoM::get_element_state(
  const std::array<size_t, 12>& idx) const
{
  Matrix<real_t, 12, 1> u_elem;
  for (int i = 0; i < 12; ++i) {
    u_elem(i) = u(idx[i]);
  }
  return u_elem;
}

std::array<real_t, 2>
EulerBeamInextensibleMoM::get_element_lambda(size_t e) const
{
  return { lambda(e), lambda(e + 1) };
}

void
EulerBeamInextensibleMoM::assemble_mass_matrix()
{
  std::vector<Triplet<real_t>> triplets;
  triplets.reserve(3 * nodes);

  for (size_t n = 0; n < nodes; ++n) {
    const real_t w = (n == 0 || n == nodes - 1) ? 0.5 * ds : ds;
    const real_t entry = mu * w;

    triplets.emplace_back(offset_x + 2 * n + 0, offset_x + 2 * n + 0, entry);
    triplets.emplace_back(offset_y + 2 * n + 0, offset_y + 2 * n + 0, entry);
    triplets.emplace_back(offset_z + 2 * n + 0, offset_z + 2 * n + 0, entry);
  }

  mass.resize(ndof, ndof);
  mass.setFromTriplets(triplets.begin(), triplets.end());
  mass.makeCompressed();
}

void
EulerBeamInextensibleMoM::assemble_bending_matrix()
{
  std::vector<Triplet<real_t>> triplets;
  triplets.reserve(elements * 3 * 16);

  const real_t xi_q[] = { 0.1127016654, 0.5, 0.8872983346 };
  const real_t w_q[] = { 0.2777777778, 0.4444444444, 0.2777777778 };

  for (size_t e = 0; e < elements; ++e) {
    const auto idx = get_element_dof_indices(e);
    Matrix<real_t, 4, 4> K_local = Matrix<real_t, 4, 4>::Zero();

    for (size_t qi = 0; qi < 3; ++qi) {
      const auto ddH =
        ELFF::FEM::CubicHermite<real_t>::second_derivs(xi_q[qi], ds);

      for (size_t a = 0; a < 4; ++a) {
        for (size_t b = 0; b < 4; ++b) {
          K_local(a, b) += EI * ddH[a] * ddH[b] * w_q[qi] * ds;
        }
      }
    }

    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = 0; b < 4; ++b) {
        const real_t value = K_local(a, b);
        if (value == 0.0) {
          continue;
        }

        triplets.emplace_back(idx[a], idx[b], value);
        triplets.emplace_back(idx[4 + a], idx[4 + b], value);
        triplets.emplace_back(idx[8 + a], idx[8 + b], value);
      }
    }
  }

  bending_matrix.resize(ndof, ndof);
  bending_matrix.setFromTriplets(triplets.begin(), triplets.end());
  bending_matrix.makeCompressed();
}

real_t
EulerBeamInextensibleMoM::update_lambda(real_t omega)
{
  static constexpr std::array<real_t, 3> xi_q = { 0.1127016654,
                                                  0.5,
                                                  0.8872983346 };
  static constexpr std::array<real_t, 3> w_q = { 0.2777777778,
                                                 0.4444444444,
                                                 0.2777777778 };

  Matrix<real_t, Dynamic, 1> lambda_rhs =
    Matrix<real_t, Dynamic, 1>::Zero(ndof_l);
  Matrix<real_t, Dynamic, 1> lambda_mass_lumped =
    Matrix<real_t, Dynamic, 1>::Zero(ndof_l);

  for (size_t e = 0; e < elements; ++e) {
    const std::vector<size_t> elem_nodes = { e, e + 1 };
    const std::vector<size_t> idx_x = { offset_x + 2 * elem_nodes[0],
                                        offset_x + 2 * elem_nodes[0] + 1,
                                        offset_x + 2 * elem_nodes[1],
                                        offset_x + 2 * elem_nodes[1] + 1 };
    const std::vector<size_t> idx_y = { offset_y + 2 * elem_nodes[0],
                                        offset_y + 2 * elem_nodes[0] + 1,
                                        offset_y + 2 * elem_nodes[1],
                                        offset_y + 2 * elem_nodes[1] + 1 };
    const std::vector<size_t> idx_z = { offset_z + 2 * elem_nodes[0],
                                        offset_z + 2 * elem_nodes[0] + 1,
                                        offset_z + 2 * elem_nodes[1],
                                        offset_z + 2 * elem_nodes[1] + 1 };
    const std::vector<size_t> idx_l = { elem_nodes[0], elem_nodes[1] };

    const std::array<real_t, 4> ux = {
      u[idx_x[0]], u[idx_x[1]], u[idx_x[2]], u[idx_x[3]]
    };
    const std::array<real_t, 4> uy = {
      u[idx_y[0]], u[idx_y[1]], u[idx_y[2]], u[idx_y[3]]
    };
    const std::array<real_t, 4> uz = {
      u[idx_z[0]], u[idx_z[1]], u[idx_z[2]], u[idx_z[3]]
    };

    std::array<real_t, 2> R_loc_l = { 0.0, 0.0 };
    std::array<real_t, 2> M_loc_lumped = { 0.0, 0.0 };

    for (size_t qi = 0; qi < 3; ++qi) {
      const real_t xi = xi_q[qi];
      const real_t w = w_q[qi];

      const auto dH = ELFF::FEM::CubicHermite<real_t>::derivs(xi, ds);
      const auto M = ELFF::FEM::LinearShape<real_t>::values(xi);

      real_t xp = 0.0;
      real_t yp = 0.0;
      real_t zp = 0.0;

      for (size_t i = 0; i < 4; i++) {
        xp += dH[i] * ux[i];
        yp += dH[i] * uy[i];
        zp += dH[i] * uz[i];
      }

      const real_t S = xp * xp + yp * yp + zp * zp - 1.0;

      for (size_t a = 0; a < 2; ++a) {
        R_loc_l[a] += S * M[a] * w * ds;
        M_loc_lumped[a] += M[a] * w * ds;
      }
    }

    for (size_t i = 0; i < 2; ++i) {
      lambda_rhs[idx_l[i]] += R_loc_l[i];
      lambda_mass_lumped[idx_l[i]] += M_loc_lumped[i];
    }
  }

  const Matrix<real_t, Dynamic, 1> delta_lambda =
    lambda_rhs.cwiseQuotient(lambda_mass_lumped);

  lambda += omega * r_penalty * delta_lambda;
  // lambda  += omega * delta_lambda;
  apply_lambda_boundary_conditions();
  return delta_lambda.norm();
}

real_t
EulerBeamInextensibleMoM::update_lambda_newmark(real_t dt, real_t omega)
{

  return update_lambda();

  ELFF_ASSERT(dt > 0.0, "Time step must be positive.\n");

  static constexpr std::array<real_t, 3> xi_q = { 0.1127016654,
                                                  0.5,
                                                  0.8872983346 };
  static constexpr std::array<real_t, 3> w_q = { 0.2777777778,
                                                 0.4444444444,
                                                 0.2777777778 };

  VectorXd lambda_rhs = VectorXd::Zero(ndof_l);
  VectorXd lambda_mass_lumped = VectorXd::Zero(ndof_l);
  VectorXd schur_rhs_correction = VectorXd::Zero(ndof_l);
  VectorXd schur_diag = VectorXd::Zero(ndof_l);
  const VectorXd primal_diag =
    jacobian.diagonal().cwiseAbs().cwiseMax(real_t(1e-12));

  for (size_t e = 0; e < elements; ++e) {
    const std::vector<size_t> elem_nodes = { e, e + 1 };
    const std::vector<size_t> idx_x = { offset_x + 2 * elem_nodes[0],
                                        offset_x + 2 * elem_nodes[0] + 1,
                                        offset_x + 2 * elem_nodes[1],
                                        offset_x + 2 * elem_nodes[1] + 1 };
    const std::vector<size_t> idx_y = { offset_y + 2 * elem_nodes[0],
                                        offset_y + 2 * elem_nodes[0] + 1,
                                        offset_y + 2 * elem_nodes[1],
                                        offset_y + 2 * elem_nodes[1] + 1 };
    const std::vector<size_t> idx_z = { offset_z + 2 * elem_nodes[0],
                                        offset_z + 2 * elem_nodes[0] + 1,
                                        offset_z + 2 * elem_nodes[1],
                                        offset_z + 2 * elem_nodes[1] + 1 };
    const std::vector<size_t> idx_l = { elem_nodes[0], elem_nodes[1] };

    const std::array<real_t, 4> ux = {
      u[idx_x[0]], u[idx_x[1]], u[idx_x[2]], u[idx_x[3]]
    };
    const std::array<real_t, 4> uy = {
      u[idx_y[0]], u[idx_y[1]], u[idx_y[2]], u[idx_y[3]]
    };
    const std::array<real_t, 4> uz = {
      u[idx_z[0]], u[idx_z[1]], u[idx_z[2]], u[idx_z[3]]
    };

    std::array<real_t, 2> R_loc_l = { 0.0, 0.0 };
    std::array<real_t, 2> M_loc_lumped = { 0.0, 0.0 };

    for (size_t qi = 0; qi < 3; ++qi) {
      const real_t xi = xi_q[qi];
      const real_t w = w_q[qi];

      const auto dH = ELFF::FEM::CubicHermite<real_t>::derivs(xi, ds);
      const auto M = ELFF::FEM::LinearShape<real_t>::values(xi);

      real_t xp = 0.0;
      real_t yp = 0.0;
      real_t zp = 0.0;

      for (size_t i = 0; i < 4; ++i) {
        xp += dH[i] * ux[i];
        yp += dH[i] * uy[i];
        zp += dH[i] * uz[i];
      }

      const real_t S = xp * xp + yp * yp + zp * zp - 1.0;

      for (size_t a = 0; a < 2; ++a) {
        const real_t weight = M[a] * w * ds;
        R_loc_l[a] += S * weight;
        M_loc_lumped[a] += weight;

        for (size_t i = 0; i < 4; ++i) {
          const real_t gx = 2.0 * xp * dH[i] * weight;
          const real_t gy = 2.0 * yp * dH[i] * weight;
          const real_t gz = 2.0 * zp * dH[i] * weight;

          schur_rhs_correction(idx_l[a]) +=
            gx * residual(idx_x[i]) / primal_diag(idx_x[i]);
          schur_rhs_correction(idx_l[a]) +=
            gy * residual(idx_y[i]) / primal_diag(idx_y[i]);
          schur_rhs_correction(idx_l[a]) +=
            gz * residual(idx_z[i]) / primal_diag(idx_z[i]);

          schur_diag(idx_l[a]) += gx * gx / primal_diag(idx_x[i]);
          schur_diag(idx_l[a]) += gy * gy / primal_diag(idx_y[i]);
          schur_diag(idx_l[a]) += gz * gz / primal_diag(idx_z[i]);
        }
      }
    }

    for (size_t i = 0; i < 2; ++i) {
      lambda_rhs[idx_l[i]] += R_loc_l[i];
      lambda_mass_lumped[idx_l[i]] += M_loc_lumped[i];
    }
  }

  const real_t dual_shift = dt;
  const VectorXd delta_lambda =
    (lambda_rhs - schur_rhs_correction)
      .cwiseQuotient(
        (schur_diag.array() + lambda_mass_lumped.array() + dual_shift)
          .matrix());

  lambda += omega * r_penalty * delta_lambda;
  apply_lambda_boundary_conditions();
  return delta_lambda.norm();
}

void
EulerBeamInextensibleMoM::solve_newmark(real_t dt,
                                        std::array<real_t, 3> load,
                                        real_t                beta,
                                        real_t                gamma)
{
  SparseLU<SparseMatrix<real_t>, COLAMDOrdering<int>> solver;

  if (time_iter == 0) {
    initialize_acceleration(load);
  }
  const std::array<real_t, 3> load_current = load;
  if (!have_prev_uniform_load) {
    load_prev = load;
    have_prev_uniform_load = true;
  }

  u_prev = u;
  const VectorXd v_old = v_prev;
  const VectorXd a_old = a_prev;

  real_t dlambda_norm = 0.0, constr_norm = 0.0, res_norm = 0.0;
  size_t iter_outer = 0, iter_inner = 0;
 
  for (iter_outer = 0; iter_outer < max_iter_outer; ++iter_outer) {
    for (iter_inner = 0; iter_inner < max_iter_inner; ++iter_inner) {
      assemble_system_newmark(dt, load, beta, gamma);
      apply_boundary_conditions();

      res_norm = residual.norm();
      if (iter_inner + 1 >= min_iter_inner && res_norm < tol_primal) {
        break;
      }

      solver.compute(jacobian);

      if (solver.info() != Success) {
        ELFF_ABORT("EulerBeamInextensibleMoM::solve_newmark(): "
                   "SparseLU factorization failed.\n");
      }

      const VectorXd delta_u = solver.solve(-residual);

      if (solver.info() != Success) {
        ELFF_ABORT("EulerBeamInextensibleMoM::solve_newmark(): "
                   "SparseLU solve failed.\n");
      }
      u += delta_u;
    }

    dlambda_norm = update_lambda_newmark(dt);
    constr_norm = compute_inextensibility_error_l2();

    ELFF_LOG(time_iter << "\t" << iter_outer << "\t" << iter_inner << "\t"
                       << res_norm << "\t" << dlambda_norm << "\t"
                       << constr_norm);

    if (iter_outer >= min_iter_outer - 1 && iter_inner >= min_iter_inner - 1) {
      if (res_norm < tol_primal && constr_norm < tol_constraint) {
        break;
      }
    }

    if (iter_outer == max_iter_outer - 1) {
      ELFF_LOG(time_iter << "\t" << iter_outer << "\t" << iter_inner << "\t"
                         << res_norm << "\t" << dlambda_norm << "\t"
                         << constr_norm);
      ELFF_ABORT("EulerBeamInextensibleMoM::solve_newmark() did not "
                 "converge.\n");
    }
  }

  update_newmark_state_from_displacement(dt, beta, gamma, v_old, a_old);
  apply_dynamic_state_boundary_conditions();
  u_prev = u;
  load_prev = load_current;
  have_prev_uniform_load = true;
  update_mesh();
  ++time_iter;
  t += dt;
}

void
EulerBeamInextensibleMoM::solve_newmark(
  real_t dt,
  std::vector<std::array<real_t, 3>> load,
  real_t beta,
  real_t gamma)
{
  ELFF_ASSERT(load.size() == nodes,
              "Size of load vector must equal number of nodes.");

  SparseLU<SparseMatrix<real_t>, COLAMDOrdering<int>> solver;

  if (time_iter == 0) {
    initialize_acceleration(load);
  }
  const std::vector<std::array<real_t, 3>> load_current = load;
  if (!have_prev_nodal_load) {
    nodal_load_prev = load;
    have_prev_nodal_load = true;
  }

  u_prev = u;
  const VectorXd v_old = v_prev;
  const VectorXd a_old = a_prev;

  real_t dlambda_norm = 0.0, constr_norm = 0.0, res_norm = 0.0;
  size_t iter_outer = 0, iter_inner = 0;

  for (iter_outer = 0; iter_outer < max_iter_outer; ++iter_outer) {
    for (iter_inner = 0; iter_inner < max_iter_inner; ++iter_inner) {
      assemble_system_newmark(dt, load, beta, gamma);
      apply_boundary_conditions();

      res_norm = residual.norm();
      if (iter_inner + 1 >= min_iter_inner && res_norm < tol_primal) {
        break;
      }
      solver.compute(jacobian);

      if (solver.info() != Success) {
        ELFF_ABORT("EulerBeamInextensibleMoM::solve_newmark(): "
                   "SparseLU factorization failed.\n");
      }

      const VectorXd delta_u = solver.solve(-residual);

      if (solver.info() != Success) {
        ELFF_ABORT("EulerBeamInextensibleMoM::solve_newmark(): "
                   "SparseLU solve failed.\n");
      }

      u += delta_u;
    }

    dlambda_norm = update_lambda_newmark(dt);
    constr_norm = compute_inextensibility_error_l2();

    ELFF_LOG(time_iter << "\t" << iter_outer << "\t" << iter_inner << "\t"
                       << res_norm << "\t" << dlambda_norm << "\t"
                       << constr_norm);

    if (iter_outer >= min_iter_outer - 1 && iter_inner >= min_iter_inner - 1) {
      if (res_norm < tol_primal && constr_norm < tol_constraint) {
        break;
      }
    }

    if (iter_outer == max_iter_outer - 1) {
      ELFF_ABORT("EulerBeamInextensibleMoM::solve_newmark() did not "
                 "converge.\n");
    }
  }

  update_newmark_state_from_displacement(dt, beta, gamma, v_old, a_old);
  apply_dynamic_state_boundary_conditions();
  u_prev = u;
  nodal_load_prev = load_current;
  have_prev_nodal_load = true;
  update_mesh();
  ++time_iter;
  t += dt;
}

void
EulerBeamInextensibleMoM::initialize_acceleration(
  std::array<real_t, 3> load)
{
  residual = assemble_residual_template<real_t>(u, load);
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
EulerBeamInextensibleMoM::initialize_acceleration(
  std::vector<std::array<real_t, 3>> load)
{
  residual = assemble_residual_template<real_t>(u, load);
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

real_t
EulerBeamInextensibleMoM::compute_inextensibility_error_l2() const
{
  const real_t xi_q[] = { 0.1127016654, 0.5, 0.8872983346 };
  const real_t w_q[] = { 0.2777777778, 0.4444444444, 0.2777777778 };

  real_t error_sq = 0.0;

  for (size_t e = 0; e < elements; ++e) {
    const std::vector<size_t> elem_nodes = { e, e + 1 };
    const std::vector<size_t> idx_x = { offset_x + 2 * elem_nodes[0],
                                        offset_x + 2 * elem_nodes[0] + 1,
                                        offset_x + 2 * elem_nodes[1],
                                        offset_x + 2 * elem_nodes[1] + 1 };
    const std::vector<size_t> idx_y = { offset_y + 2 * elem_nodes[0],
                                        offset_y + 2 * elem_nodes[0] + 1,
                                        offset_y + 2 * elem_nodes[1],
                                        offset_y + 2 * elem_nodes[1] + 1 };
    const std::vector<size_t> idx_z = { offset_z + 2 * elem_nodes[0],
                                        offset_z + 2 * elem_nodes[0] + 1,
                                        offset_z + 2 * elem_nodes[1],
                                        offset_z + 2 * elem_nodes[1] + 1 };

    const std::array<real_t, 4> ux = {
      u[idx_x[0]], u[idx_x[1]], u[idx_x[2]], u[idx_x[3]]
    };
    const std::array<real_t, 4> uy = {
      u[idx_y[0]], u[idx_y[1]], u[idx_y[2]], u[idx_y[3]]
    };
    const std::array<real_t, 4> uz = {
      u[idx_z[0]], u[idx_z[1]], u[idx_z[2]], u[idx_z[3]]
    };

    for (size_t qi = 0; qi < 3; ++qi) {
      const real_t xi = xi_q[qi];
      const real_t w = w_q[qi];
      const auto dH = ELFF::FEM::CubicHermite<real_t>::derivs(xi, ds);

      real_t xp = 0.0;
      real_t yp = 0.0;
      real_t zp = 0.0;

      for (size_t i = 0; i < 4; ++i) {
        xp += dH[i] * ux[i];
        yp += dH[i] * uy[i];
        zp += dH[i] * uz[i];
      }

      const real_t S = xp * xp + yp * yp + zp * zp - 1.0;
      error_sq += S * S * w * ds;
    }
  }

  return std::sqrt(error_sq);
}

real_t
EulerBeamInextensibleMoM::compute_inextensibility_error_linf() const
{
  const real_t xi_q[] = { 0.1127016654, 0.5, 0.8872983346 };

  real_t max_error = 0.0;

  for (size_t e = 0; e < elements; ++e) {
    const std::vector<size_t> elem_nodes = { e, e + 1 };
    const std::vector<size_t> idx_x = { offset_x + 2 * elem_nodes[0],
                                        offset_x + 2 * elem_nodes[0] + 1,
                                        offset_x + 2 * elem_nodes[1],
                                        offset_x + 2 * elem_nodes[1] + 1 };
    const std::vector<size_t> idx_y = { offset_y + 2 * elem_nodes[0],
                                        offset_y + 2 * elem_nodes[0] + 1,
                                        offset_y + 2 * elem_nodes[1],
                                        offset_y + 2 * elem_nodes[1] + 1 };
    const std::vector<size_t> idx_z = { offset_z + 2 * elem_nodes[0],
                                        offset_z + 2 * elem_nodes[0] + 1,
                                        offset_z + 2 * elem_nodes[1],
                                        offset_z + 2 * elem_nodes[1] + 1 };

    const std::array<real_t, 4> ux = {
      u[idx_x[0]], u[idx_x[1]], u[idx_x[2]], u[idx_x[3]]
    };
    const std::array<real_t, 4> uy = {
      u[idx_y[0]], u[idx_y[1]], u[idx_y[2]], u[idx_y[3]]
    };
    const std::array<real_t, 4> uz = {
      u[idx_z[0]], u[idx_z[1]], u[idx_z[2]], u[idx_z[3]]
    };

    for (size_t qi = 0; qi < 3; ++qi) {
      const real_t xi = xi_q[qi];
      const auto dH = ELFF::FEM::CubicHermite<real_t>::derivs(xi, ds);

      real_t xp = 0.0;
      real_t yp = 0.0;
      real_t zp = 0.0;

      for (size_t i = 0; i < 4; ++i) {
        xp += dH[i] * ux[i];
        yp += dH[i] * uy[i];
        zp += dH[i] * uz[i];
      }

      const real_t S = xp * xp + yp * yp + zp * zp - 1.0;
      max_error = std::max(max_error, std::abs(S));
    }
  }

  return max_error;
}

void
EulerBeamInextensibleMoM::update_mesh()
{
  const size_t nodes = this->mesh.get_nodes();

  std::vector<std::array<real_t, 3>>& centerline = mesh.get_centerline();
  std::vector<std::array<real_t, 3>>& velocity = mesh.get_centerline_velocity();
  std::vector<std::array<real_t, 3>>& slope = mesh.get_slope();
  std::vector<real_t>& s = mesh.get_curvilinear_axis();
  (void)s;

  for (size_t i = 0; i < nodes; ++i) {
    centerline[i][0] = u(offset_x + 2 * i);
    centerline[i][1] = u(offset_y + 2 * i);
    centerline[i][2] = u(offset_z + 2 * i);
    slope[i][0] = u(offset_x + 2 * i + 1);
    slope[i][1] = u(offset_y + 2 * i + 1);
    slope[i][2] = u(offset_z + 2 * i + 1);
    velocity[i][0] = v_prev(offset_x + 2 * i);
    velocity[i][1] = v_prev(offset_y + 2 * i);
    velocity[i][2] = v_prev(offset_z + 2 * i);
  }
}

void
EulerBeamInextensibleMoM::assemble_system(
  std::array<real_t, 3> load,
  real_t                bending_scale)
{
  using ADDeriv = Matrix<real_t, 12, 1>;
  using AD = AutoDiffScalar<ADDeriv>;
  using ADVec = Matrix<AD, 12, 1>;
  using Tpl = Triplet<real_t>;

  residual = VectorXd::Zero(ndof);
  std::vector<Tpl> triplets;
  triplets.reserve(elements * 12 * 12);

  for (size_t e = 0; e < elements; ++e) {
    const auto idx = get_element_dof_indices(e);
    const auto lambda_elem = get_element_lambda(e);
    const Matrix<real_t, 12, 1> u_elem = get_element_state(idx);

    ADVec u_ad;
    for (int a = 0; a < 12; ++a) {
      ADDeriv seed = ADDeriv::Zero();
      seed(a) = 1.0;
      u_ad(a) = AD(u_elem(a), seed);
    }

    const ADVec R_loc_ad =
      assemble_element_residual_template<AD>(u_ad, lambda_elem, load, bending_scale);

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
EulerBeamInextensibleMoM::assemble_system(
  std::vector<std::array<real_t, 3>> load,
  real_t                              bending_scale)
{
  using ADDeriv = Matrix<real_t, 12, 1>;
  using AD = AutoDiffScalar<ADDeriv>;
  using ADVec = Matrix<AD, 12, 1>;
  using Tpl = Triplet<real_t>;

  residual = VectorXd::Zero(ndof);
  std::vector<Tpl> triplets;
  triplets.reserve(elements * 12 * 12);

  for (size_t e = 0; e < elements; ++e) {
    const auto idx = get_element_dof_indices(e);
    const auto lambda_elem = get_element_lambda(e);
    const Matrix<real_t, 12, 1> u_elem = get_element_state(idx);
    const std::array<std::array<real_t, 3>, 2> load_elem = { load[e],
                                                             load[e + 1] };

    ADVec u_ad;
    for (int a = 0; a < 12; ++a) {
      ADDeriv seed = ADDeriv::Zero();
      seed(a) = 1.0;
      u_ad(a) = AD(u_elem(a), seed);
    }

    const ADVec R_loc_ad = assemble_element_residual_template<AD>(
      u_ad, lambda_elem, load_elem, bending_scale);

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
EulerBeamInextensibleMoM::apply_boundary_conditions()
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
        idx = {
          offset_x + 2 * ni + 0,
          offset_y + 2 * ni + 0,
          offset_z + 2 * ni + 0,
        };
        vals = { bcvals.force[0], bcvals.force[1], bcvals.force[2] };
        break;
      case point_torque_bc:
        idx = {
          offset_x + 2 * ni + 1,
          offset_y + 2 * ni + 1,
          offset_z + 2 * ni + 1,
        };
        vals = { bcvals.torque[0], bcvals.torque[1], bcvals.torque[2] };
        break;
      default:
        break;
    }

    switch (bctype) {
      case point_force_bc:
      case point_torque_bc:
        for (size_t i = 0; i < vals.size(); i++) {
          residual[idx[i]] -= vals[i];
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
              if (it.row() == i) {
                it.valueRef() = 0.0;
              }
            }
          }

          jacobian.coeffRef(i, i) = 1.0;
        }

        for (size_t i = 0; i < vals.size(); i++) {
          residual[idx[i]] = u[idx[i]] - vals[i];
        }
        break;
    }
  }
}

void
EulerBeamInextensibleMoM::apply_lambda_boundary_conditions()
{
  for (size_t bi = 0; bi < 2; ++bi) {
    if (boundary_conditions.type[bi] != free_bc ) {
      continue;
    }

    switch (boundary_conditions.end[bi]) {
      case left:
        lambda(0) = 0.0;
        break;
      case right:
        lambda(nodes - 1) = 0.0;
        break;
    }
  }
}

void
EulerBeamInextensibleMoM::assemble_system_newmark(
  real_t dt,
  std::array<real_t, 3> load,
  real_t beta,
  real_t gamma)
{
  static_cast<void>(gamma);
  add_averaged_uniform_load(load);
  assemble_system(load, 0.0);
  add_midpoint_bending_terms();
  add_newmark_inertial_terms(dt, beta);
}

void
EulerBeamInextensibleMoM::assemble_system_newmark(
  real_t dt,
  std::vector<std::array<real_t, 3>> load,
  real_t beta,
  real_t gamma)
{
  static_cast<void>(gamma);
  add_averaged_nodal_load(load);
  assemble_system(load, 0.0);
  add_midpoint_bending_terms();
  add_newmark_inertial_terms(dt, beta);
}

void
EulerBeamInextensibleMoM::add_midpoint_bending_terms()
{
  residual.noalias() += 0.5 * (bending_matrix * u);
  residual.noalias() += 0.5 * (bending_matrix * u_prev);
  jacobian += 0.5 * bending_matrix;
  jacobian.makeCompressed();
}

void
EulerBeamInextensibleMoM::add_averaged_uniform_load(
  std::array<real_t, 3>& load) const
{
  if (!have_prev_uniform_load) {
    return;
  }

  for (size_t i = 0; i < 3; ++i) {
    load[i] = 0.5 * (load_prev[i] + load[i]);
  }
}

void
EulerBeamInextensibleMoM::add_averaged_nodal_load(
  std::vector<std::array<real_t, 3>>& load) const
{
  if (!have_prev_nodal_load) {
    return;
  }

  ELFF_ASSERT(load.size() == nodal_load_prev.size(),
              "Stored nodal load history must match node count.");

  for (size_t i = 0; i < load.size(); ++i) {
    for (size_t d = 0; d < 3; ++d) {
      load[i][d] = 0.5 * (nodal_load_prev[i][d] + load[i][d]);
    }
  }
}

void
EulerBeamInextensibleMoM::add_newmark_inertial_terms(real_t dt,
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
EulerBeamInextensibleMoM::update_newmark_state_from_displacement(
  real_t dt,
  real_t beta,
  real_t gamma,
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
EulerBeamInextensibleMoM::apply_dynamic_state_boundary_conditions()
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
EulerBeamInextensibleMoM::project_velocity_onto_constraint_manifold()
{
  static constexpr std::array<real_t, 3> xi_q = { 0.1127016654,
                                                  0.5,
                                                  0.8872983346 };
  static constexpr std::array<real_t, 3> w_q = { 0.2777777778,
                                                 0.4444444444,
                                                 0.2777777778 };

  VectorXd eta_rhs = VectorXd::Zero(ndof_l);
  VectorXd eta_mass_lumped = VectorXd::Zero(ndof_l);
  VectorXd velocity_mass_diag = VectorXd::Zero(ndof);
  std::vector<bool> skip_projection_dof(ndof, false);

  for (size_t bi = 0; bi < 2; ++bi) {
    const EulerBeamBCType bctype = boundary_conditions.type[bi];
    size_t ni = 0;

    switch (boundary_conditions.end[bi]) {
      case left:
        ni = 0;
        break;
      case right:
        ni = nodes - 1;
        break;
    }

    if (bctype == free_bc || bctype == simple_bc || bctype == clamped_bc) {
      skip_projection_dof[offset_x + 2 * ni + 0] = true;
      skip_projection_dof[offset_y + 2 * ni + 0] = true;
      skip_projection_dof[offset_z + 2 * ni + 0] = true;
    }
  }

  for (size_t e = 0; e < elements; ++e) {
    const size_t n0 = e;
    const size_t n1 = e + 1;
    const std::array<size_t, 4> idx_x = { offset_x + 2 * n0,
                                          offset_x + 2 * n0 + 1,
                                          offset_x + 2 * n1,
                                          offset_x + 2 * n1 + 1 };
    const std::array<size_t, 4> idx_y = { offset_y + 2 * n0,
                                          offset_y + 2 * n0 + 1,
                                          offset_y + 2 * n1,
                                          offset_y + 2 * n1 + 1 };
    const std::array<size_t, 4> idx_z = { offset_z + 2 * n0,
                                          offset_z + 2 * n0 + 1,
                                          offset_z + 2 * n1,
                                          offset_z + 2 * n1 + 1 };
    const std::array<size_t, 2> idx_l = { n0, n1 };

    const std::array<real_t, 4> ux = {
      u(idx_x[0]), u(idx_x[1]), u(idx_x[2]), u(idx_x[3])
    };
    const std::array<real_t, 4> uy = {
      u(idx_y[0]), u(idx_y[1]), u(idx_y[2]), u(idx_y[3])
    };
    const std::array<real_t, 4> uz = {
      u(idx_z[0]), u(idx_z[1]), u(idx_z[2]), u(idx_z[3])
    };

    const std::array<real_t, 4> vx = {
      v_prev(idx_x[0]), v_prev(idx_x[1]), v_prev(idx_x[2]), v_prev(idx_x[3])
    };
    const std::array<real_t, 4> vy = {
      v_prev(idx_y[0]), v_prev(idx_y[1]), v_prev(idx_y[2]), v_prev(idx_y[3])
    };
    const std::array<real_t, 4> vz = {
      v_prev(idx_z[0]), v_prev(idx_z[1]), v_prev(idx_z[2]), v_prev(idx_z[3])
    };

    std::array<real_t, 2> R_loc = { 0.0, 0.0 };
    std::array<real_t, 2> M_loc_lumped = { 0.0, 0.0 };

    for (size_t qi = 0; qi < xi_q.size(); ++qi) {
      const real_t xi = xi_q[qi];
      const real_t w = w_q[qi];

      const auto H = ELFF::FEM::CubicHermite<real_t>::values(xi, ds);
      const auto dH = ELFF::FEM::CubicHermite<real_t>::derivs(xi, ds);
      const auto M = ELFF::FEM::LinearShape<real_t>::values(xi);

      real_t xp = 0.0, yp = 0.0, zp = 0.0;
      real_t vxp = 0.0, vyp = 0.0, vzp = 0.0;

      for (size_t i = 0; i < 4; ++i) {
        xp += dH[i] * ux[i];
        yp += dH[i] * uy[i];
        zp += dH[i] * uz[i];
        vxp += dH[i] * vx[i];
        vyp += dH[i] * vy[i];
        vzp += dH[i] * vz[i];
      }
      for (size_t a = 0; a < 2; ++a) {
        const size_t i = 2 * a;
        const real_t mass_q = H[i] * H[i] * w * ds;
        velocity_mass_diag(idx_x[i]) += mass_q;
        velocity_mass_diag(idx_y[i]) += mass_q;
        velocity_mass_diag(idx_z[i]) += mass_q;
      }

      const real_t Sv = 2.0 * (xp * vxp + yp * vyp + zp * vzp);

      for (size_t a = 0; a < 2; ++a) {
        R_loc[a] += Sv * M[a] * w * ds;
        M_loc_lumped[a] += M[a] * w * ds;
      }
    }

    for (size_t a = 0; a < 2; ++a) {
      eta_rhs(idx_l[a]) += R_loc[a];
      eta_mass_lumped(idx_l[a]) += M_loc_lumped[a];
    }
  }

  const VectorXd eta =
    eta_rhs.cwiseQuotient(eta_mass_lumped.array().max(1e-14).matrix());

  VectorXd velocity_correction = VectorXd::Zero(ndof);

  for (size_t e = 0; e < elements; ++e) {
    const size_t n0 = e;
    const size_t n1 = e + 1;
    const std::array<size_t, 4> idx_x = { offset_x + 2 * n0,
                                          offset_x + 2 * n0 + 1,
                                          offset_x + 2 * n1,
                                          offset_x + 2 * n1 + 1 };
    const std::array<size_t, 4> idx_y = { offset_y + 2 * n0,
                                          offset_y + 2 * n0 + 1,
                                          offset_y + 2 * n1,
                                          offset_y + 2 * n1 + 1 };
    const std::array<size_t, 4> idx_z = { offset_z + 2 * n0,
                                          offset_z + 2 * n0 + 1,
                                          offset_z + 2 * n1,
                                          offset_z + 2 * n1 + 1 };
    const std::array<size_t, 2> idx_l = { n0, n1 };
    const std::array<real_t, 4> ux = {
      u(idx_x[0]), u(idx_x[1]), u(idx_x[2]), u(idx_x[3])
    };
    const std::array<real_t, 4> uy = {
      u(idx_y[0]), u(idx_y[1]), u(idx_y[2]), u(idx_y[3])
    };
    const std::array<real_t, 4> uz = {
      u(idx_z[0]), u(idx_z[1]), u(idx_z[2]), u(idx_z[3])
    };

    for (size_t qi = 0; qi < xi_q.size(); ++qi) {
      const real_t xi = xi_q[qi];
      const real_t w = w_q[qi];

      const auto dH = ELFF::FEM::CubicHermite<real_t>::derivs(xi, ds);
      const auto M = ELFF::FEM::LinearShape<real_t>::values(xi);

      real_t xp = 0.0, yp = 0.0, zp = 0.0;
      for (size_t i = 0; i < 4; ++i) {
        xp += dH[i] * ux[i];
        yp += dH[i] * uy[i];
        zp += dH[i] * uz[i];
      }

      real_t eta_q = 0.0;
      for (size_t a = 0; a < 2; ++a) {
        eta_q += M[a] * eta(idx_l[a]);
      }

      for (size_t a = 0; a < 2; ++a) {
        const size_t i = 2 * a;
        const real_t coeff = 2.0 * eta_q * dH[i] * w * ds;
        velocity_correction(idx_x[i]) += coeff * xp;
        velocity_correction(idx_y[i]) += coeff * yp;
        velocity_correction(idx_z[i]) += coeff * zp;
      }
    }
  }

  for (Eigen::Index i = 0; i < v_prev.size(); ++i) {
    if (skip_projection_dof[static_cast<size_t>(i)]) {
      continue;
    }

    const real_t lumped_mass =
      mu * std::max<real_t>(velocity_mass_diag(i), 1e-14);
    v_prev(i) -= velocity_correction(i) / lumped_mass;
  }
}

} // namespace Models
} // namespace ELFF
