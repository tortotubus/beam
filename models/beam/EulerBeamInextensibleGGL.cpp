
#include "elff/general/error.hpp"
#include "elff/models/beam/EulerBeamInextensibleGGL.hpp"

namespace ELFF {
namespace Models {

namespace {
constexpr std::array<real_t, 3> xi_q = { 0.1127016654,
                                         0.5,
                                         0.8872983346 };
constexpr std::array<real_t, 3> w_q  = { 0.2777777778,
                                         0.4444444444,
                                         0.2777777778 };

struct GGLLinearizedElementData
{
  Matrix<real_t, 16, 1> residual;
  Matrix<real_t, 12, 1> velocity_residual;
  Matrix<real_t, 12, 12> displacement_jacobian;
  Matrix<real_t, 12, 12> velocity_displacement_jacobian;
  Matrix<real_t, 2, 12> lambda_displacement_jacobian;
  Matrix<real_t, 2, 12> mu_displacement_jacobian;
  Matrix<real_t, 2, 12> mu_velocity_jacobian;

  void reset(const Matrix<real_t, 12, 12>& local_bending_jacobian)
  {
    residual.setZero();
    velocity_residual.setZero();
    displacement_jacobian = local_bending_jacobian;
    velocity_displacement_jacobian.setZero();
    lambda_displacement_jacobian.setZero();
    mu_displacement_jacobian.setZero();
    mu_velocity_jacobian.setZero();
  }
};

void assemble_ggl_element_linearized_data(
  real_t                                      EI,
  real_t                                      ds,
  const Matrix<real_t, 12, 1>&                u_elem,
  const std::array<real_t, 2>&                lambda_elem,
  const std::array<real_t, 2>&                mu_elem,
  const std::array<real_t, 4>&                vx_elem,
  const std::array<real_t, 4>&                vy_elem,
  const std::array<real_t, 4>&                vz_elem,
  const std::array<std::array<real_t, 3>, 2>& load_elem,
  const std::array<std::array<real_t, 4>, 3>& quad_H,
  const std::array<std::array<real_t, 4>, 3>& quad_dH,
  const std::array<std::array<real_t, 4>, 3>& quad_ddH,
  const std::array<std::array<real_t, 2>, 3>& quad_M,
  const Matrix<real_t, 12, 12>&               local_bending_jacobian,
  GGLLinearizedElementData&                    data)
{
  data.reset(local_bending_jacobian);

  for (size_t qi = 0; qi < xi_q.size(); ++qi) {
    const real_t wds = w_q[qi] * ds;
    const auto& H   = quad_H[qi];
    const auto& dH  = quad_dH[qi];
    const auto& ddH = quad_ddH[qi];
    const auto& M   = quad_M[qi];

    real_t xp = 0.0, yp = 0.0, zp = 0.0;
    real_t xpp = 0.0, ypp = 0.0, zpp = 0.0;
    for (size_t i = 0; i < 4; ++i) {
      xp  += dH[i]  * u_elem(i);
      yp  += dH[i]  * u_elem(i + 4);
      zp  += dH[i]  * u_elem(i + 8);
      xpp += ddH[i] * u_elem(i);
      ypp += ddH[i] * u_elem(i + 4);
      zpp += ddH[i] * u_elem(i + 8);
    }

    real_t vxp = 0, vyp = 0, vzp = 0;
    for (size_t i = 0; i < 4; ++i) {
      vxp += dH[i] * vx_elem[i];
      vyp += dH[i] * vy_elem[i];
      vzp += dH[i] * vz_elem[i];
    }

    const real_t lambda_q = M[0] * lambda_elem[0] + M[1] * lambda_elem[1];
    const real_t mu_q     = M[0] * mu_elem[0]     + M[1] * mu_elem[1];
    const real_t fx       = M[0] * load_elem[0][0] + M[1] * load_elem[1][0];
    const real_t fy       = M[0] * load_elem[0][1] + M[1] * load_elem[1][1];
    const real_t fz       = M[0] * load_elem[0][2] + M[1] * load_elem[1][2];

    const std::array<real_t, 3> slope = { xp, yp, zp };
    const std::array<real_t, 3> curvature = { xpp, ypp, zpp };
    const std::array<real_t, 3> velocity_slope = { vxp, vyp, vzp };
    const std::array<real_t, 3> load_q = { fx, fy, fz };

    for (size_t axis = 0; axis < 3; ++axis) {
      const size_t row_offset = 4 * axis;
      for (size_t a = 0; a < 4; ++a) {
        const real_t weighted_dH = dH[a] * wds;
        data.residual(row_offset + a) += EI * curvature[axis] * ddH[a] * wds;
        data.residual(row_offset + a) -=
          2.0 * lambda_q * slope[axis] * weighted_dH;
        data.residual(row_offset + a) -= load_q[axis] * H[a] * wds;

        data.velocity_residual(row_offset + a) -=
          mu_q * slope[axis] * weighted_dH;

        for (size_t b = 0; b < 4; ++b) {
          const real_t dHdH = dH[a] * dH[b] * wds;
          data.displacement_jacobian(row_offset + a, row_offset + b) -=
            2.0 * lambda_q * dHdH;
          data.velocity_displacement_jacobian(row_offset + a,
                                              row_offset + b) -= mu_q * dHdH;
        }
      }
    }

    const real_t g  = xp * xp + yp * yp + zp * zp - 1.0;
    const real_t Gv = xp * vxp + yp * vyp + zp * vzp;
    for (size_t a = 0; a < 2; ++a) {
      const real_t weighted_M = M[a] * wds;
      data.residual(12 + a) += g * weighted_M;
      data.residual(14 + a) += Gv * weighted_M;

      for (size_t axis = 0; axis < 3; ++axis) {
        const size_t row_offset = 4 * axis;
        for (size_t i = 0; i < 4; ++i) {
          const real_t coeff = weighted_M * dH[i];
          data.lambda_displacement_jacobian(a, row_offset + i) +=
            2.0 * slope[axis] * coeff;
          data.mu_displacement_jacobian(a, row_offset + i) +=
            velocity_slope[axis] * coeff;
          data.mu_velocity_jacobian(a, row_offset + i) +=
            slope[axis] * coeff;
        }
      }
    }
  }
}

} // namespace

// -----------------------------------------------------------------------
// Construction and initial conditions
// -----------------------------------------------------------------------
EulerBeamInextensibleGGL::EulerBeamInextensibleGGL(
  real_t                  length,
  real_t                  EI,
  real_t                  mu_density,
  size_t                  n_nodes,
  EulerBeam::EulerBeamBCs bcs,
  real_t                  r_penalty)
  : EulerBeam(length, EI, mu_density, n_nodes, bcs)
  , elements(n_nodes - 1)
  , nodes(n_nodes)
  , ds(mesh.get_ds())
  , ndof_l(nodes)
  , offset_x(0)
  , offset_y(2 * nodes)
  , offset_z(4 * nodes)
  , ndof(6 * nodes)
  , r_penalty(r_penalty)
  , u(VectorXd::Zero(ndof))
  , u_prev(VectorXd::Zero(ndof))
  , v_prev(VectorXd::Zero(ndof))
  , a_prev(VectorXd::Zero(ndof))
  , lambda(VectorXd::Zero(ndof_l))
  , lambda_prev(VectorXd::Zero(ndof_l))
  , mu_prev(VectorXd::Zero(ndof_l))
  , tension(VectorXd::Zero(ndof_l))
  , mass(SparseMatrix<real_t>(ndof, ndof))
  , offset_v(ndof)
  , offset_lambda(2 * ndof)
  , offset_mu(2 * ndof + ndof_l)
  , ggl_dof(2 * ndof + 2 * ndof_l)
  , ggl_residual(VectorXd::Zero(ggl_dof))
  , ggl_jacobian(SparseMatrix<real_t>(ggl_dof, ggl_dof))
  , newmark_coeff_a(0.0)
  , newmark_coeff_v(0.0)
  , newmark_gamma(0.5)
  , newmark_beta(0.25)
  , u_tilde(VectorXd::Zero(ndof))
  , v_tilde(VectorXd::Zero(ndof))
  , quad_H()
  , quad_dH()
  , quad_ddH()
  , quad_M()
  , element_disp_dof_indices_cache(elements)
  , local_bending_jacobian(Matrix<real_t, 12, 12>::Zero())
  , max_newton(25)
  , tol_newton(5e-5)
{
  dimension = 3;
  apply_initial_condition(mesh);
  u_prev = u;
  initialize_quadrature_cache();
  initialize_element_dof_cache();
  initialize_constant_element_matrices();
  assemble_mass_matrix();
}

void
EulerBeamInextensibleGGL::apply_initial_condition()
{
  for (size_t i = 0; i < nodes; ++i) {
    u(offset_x + 2 * i + 0) = ds * i;
    u(offset_x + 2 * i + 1) = 1.0;
    u(offset_y + 2 * i + 0) = 0.0;
    u(offset_y + 2 * i + 1) = 0.0;
    u(offset_z + 2 * i + 0) = 0.0;
    u(offset_z + 2 * i + 1) = 0.0;
    lambda(i) = 0.0;
  }

  update_mesh();
  u_prev = u;
  v_prev.setZero();
  a_prev.setZero();
  lambda_prev.setZero();
  mu_prev.setZero();
  tension.setZero();
  time_iter = 0;
  t = 0.0;
}

void
EulerBeamInextensibleGGL::apply_initial_condition(EulerBeamMesh& bmesh)
{
  ELFF_ASSERT(nodes == bmesh.get_nodes(),
              "Provided mesh must have same node count as the beam mesh.\n");

  const auto centerline = bmesh.get_centerline();
  const auto slopes = bmesh.get_slope();

  for (size_t ni = 0; ni < nodes; ++ni) {
    u(offset_x + 2 * ni + 0) = centerline[ni][0];
    u(offset_x + 2 * ni + 1) = slopes[ni][0];
    u(offset_y + 2 * ni + 0) = centerline[ni][1];
    u(offset_y + 2 * ni + 1) = slopes[ni][1];
    u(offset_z + 2 * ni + 0) = centerline[ni][2];
    u(offset_z + 2 * ni + 1) = slopes[ni][2];
    lambda(ni) = 0.0;
  }

  update_mesh();
  u_prev = u;
  v_prev.setZero();

  const auto& vel = bmesh.get_centerline_velocity();
  for (size_t i = 0; i < nodes; ++i) {
    v_prev(offset_x + 2 * i) = vel[i][0];
    v_prev(offset_y + 2 * i) = vel[i][1];
    v_prev(offset_z + 2 * i) = vel[i][2];
  }

  a_prev.setZero();
  lambda_prev.setZero();
  mu_prev.setZero();
  tension.setZero();
  time_iter = 0;
  t = 0.0;
}

// -----------------------------------------------------------------------
// Public time stepping
// -----------------------------------------------------------------------
void
EulerBeamInextensibleGGL::solve(real_t dt, std::array<real_t, 3> load)
{
  if (!(dt > 0.0)) {
    throw std::runtime_error("GGL Newmark: dt must be > 0");
  }

  const real_t alpha = 0.0;
  newmark_gamma = 0.5 - alpha;
  newmark_beta = 0.25 * (1.0 - alpha) * (1.0 - alpha);

  update_newmark_predictors(dt, newmark_beta, newmark_gamma);

  VectorXd u_cur = u_prev;
  VectorXd v_cur = v_prev;
  VectorXd lambda_cur = lambda_prev;
  VectorXd mu_cur = mu_prev;
  size_t final_iter = 0;
  real_t final_res_norm = 0.0;

  SparseLU<SparseMatrix<real_t>, COLAMDOrdering<int>> solver;
  bool pattern_analyzed = false;

  // Start Newton from the previous converged state.
  for (size_t iter = 0; iter < max_newton; ++iter) {
    assemble_ggl_system(u_cur, v_cur, lambda_cur, mu_cur, load);

    const real_t res_norm = ggl_residual.norm();
    final_iter = iter;
    final_res_norm = res_norm;

    if (res_norm < tol_newton) {
      break;
    }

    if (iter == max_newton - 1) {
      ELFF_WARNING("GGL Newton did not converge at step "
                   << time_iter << ": residual = " << res_norm);
    }

    if (!pattern_analyzed) {
      // The graph is fixed across Newton iterations, so only the numeric
      // factorization needs to be refreshed after the first pass.
      solver.analyzePattern(ggl_jacobian);
      pattern_analyzed = true;
    }

    solver.factorize(ggl_jacobian);
    if (solver.info() != Success) {
      ELFF_ABORT("GGL: LU factorization failed at step "
                 << time_iter << "\n");
    }

    const VectorXd delta = solver.solve(-ggl_residual);
    if (solver.info() != Success) {
      ELFF_ABORT("GGL: back-solve failed at step " << time_iter << "\n");
    }

    u_cur += delta.head(ndof);
    v_cur += delta.segment(offset_v, ndof);
    lambda_cur += delta.segment(offset_lambda, ndof_l);
    mu_cur += delta.segment(offset_mu, ndof_l);
  }

  tension = lambda_cur;
  lambda = lambda_cur;

  const auto [inext_l2, inext_lmax] = compute_inextensibility_error(u_cur);
  ELFF_LOG(time_iter << "\t" << final_iter
                     << "\t" << final_res_norm
                     << "\t|g|_l2=" << inext_l2
                     << "\t|g|_max=" << inext_lmax);

  // Commit the converged state and refresh mesh-level output fields.
  update_newmark_state(u_cur, v_cur);
  apply_dynamic_state_boundary_conditions();

  u_prev = u;
  lambda_prev = lambda_cur;
  mu_prev = mu_cur;

  update_mesh();
  ++time_iter;
  t += dt;
}

void
EulerBeamInextensibleGGL::solve(
  real_t                              dt,
  std::vector<std::array<real_t, 3>> load)
{
  ELFF_ASSERT(load.size() == nodes,
              "Load vector size must equal node count.\n");

  if (!(dt > 0.0)) {
    throw std::runtime_error("GGL Newmark: dt must be > 0");
  }

  const real_t alpha = 0.0;
  newmark_gamma = 0.5 - alpha;
  newmark_beta = 0.25 * (1.0 - alpha) * (1.0 - alpha);

  update_newmark_predictors(dt, newmark_beta, newmark_gamma);

  VectorXd u_cur = u_prev;
  VectorXd v_cur = v_prev;
  VectorXd lambda_cur = lambda_prev;
  VectorXd mu_cur = mu_prev;
  size_t final_iter = 0;
  real_t final_res_norm = 0.0;

  SparseLU<SparseMatrix<real_t>, COLAMDOrdering<int>> solver;
  bool pattern_analyzed = false;

  for (size_t iter = 0; iter < max_newton; ++iter) {
    assemble_ggl_system(u_cur, v_cur, lambda_cur, mu_cur, load);

    const real_t res_norm = ggl_residual.norm();
    final_iter = iter;
    final_res_norm = res_norm;

    if (res_norm < tol_newton) {
      break;
    }

    if (iter == max_newton - 1) {
      ELFF_WARNING("GGL Newton did not converge at step "
                   << time_iter << ": residual = " << res_norm);
    }

    if (!pattern_analyzed) {
      solver.analyzePattern(ggl_jacobian);
      pattern_analyzed = true;
    }

    solver.factorize(ggl_jacobian);
    if (solver.info() != Success) {
      ELFF_ABORT("GGL: LU factorization failed at step "
                 << time_iter << "\n");
    }

    const VectorXd delta = solver.solve(-ggl_residual);
    if (solver.info() != Success) {
      ELFF_ABORT("GGL: back-solve failed at step " << time_iter << "\n");
    }

    u_cur += delta.head(ndof);
    v_cur += delta.segment(offset_v, ndof);
    lambda_cur += delta.segment(offset_lambda, ndof_l);
    mu_cur += delta.segment(offset_mu, ndof_l);
  }

  tension = lambda_cur;
  lambda = lambda_cur;

  const auto [inext_l2, inext_lmax] = compute_inextensibility_error(u_cur);
  ELFF_LOG(time_iter << "\t" << final_iter
                     << "\t" << final_res_norm
                     << "\t|g|_l2=" << inext_l2
                     << "\t|g|_max=" << inext_lmax);

  update_newmark_state(u_cur, v_cur);
  apply_dynamic_state_boundary_conditions();

  u_prev = u;
  lambda_prev = lambda_cur;
  mu_prev = mu_cur;

  update_mesh();
  ++time_iter;
  t += dt;
}

// -----------------------------------------------------------------------
// Newmark predictor / recovery helpers
// -----------------------------------------------------------------------
void
EulerBeamInextensibleGGL::update_newmark_predictors(
  real_t dt, real_t beta, real_t gamma)
{
  newmark_coeff_a = 1.0 / (beta * dt * dt);
  newmark_coeff_v = gamma / (beta * dt);

  u_tilde = u_prev
          + dt * v_prev
          + dt * dt * (0.5 - beta) * a_prev;
  v_tilde = v_prev + dt * (1.0 - gamma) * a_prev;
}

VectorXd
EulerBeamInextensibleGGL::compute_velocity(const VectorXd& u_cur) const
{
  // Recover the velocity implied by the current displacement iterate.
  return v_tilde + newmark_coeff_v * (u_cur - u_tilde);
}

VectorXd
EulerBeamInextensibleGGL::compute_acceleration(
  const VectorXd& u_cur) const
{
  return newmark_coeff_a * (u_cur - u_tilde);
}

std::pair<real_t, real_t>
EulerBeamInextensibleGGL::compute_inextensibility_error(
  const VectorXd& u_cur) const
{
  real_t max_abs = 0.0;
  real_t sum_sq = 0.0;
  size_t sample_count = 0;

  for (size_t e = 0; e < elements; ++e) {
    const auto& idx = get_element_disp_dof_indices(e);

    Matrix<real_t, 12, 1> u_elem;
    for (int i = 0; i < 12; ++i) {
      u_elem(i) = u_cur(idx[i]);
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

// -----------------------------------------------------------------------
// Cached quadrature and constant matrices
// -----------------------------------------------------------------------
void
EulerBeamInextensibleGGL::initialize_quadrature_cache()
{
  for (size_t qi = 0; qi < xi_q.size(); ++qi) {
    quad_H[qi] = ELFF::FEM::CubicHermite<real_t>::values(xi_q[qi], ds);
    quad_dH[qi] = ELFF::FEM::CubicHermite<real_t>::derivs(xi_q[qi], ds);
    quad_ddH[qi] =
      ELFF::FEM::CubicHermite<real_t>::second_derivs(xi_q[qi], ds);
    quad_M[qi] = ELFF::FEM::LinearShape<real_t>::values(xi_q[qi]);
  }
}

void
EulerBeamInextensibleGGL::initialize_element_dof_cache()
{
  for (size_t e = 0; e < elements; ++e) {
    const size_t n0 = e;
    const size_t n1 = e + 1;
    element_disp_dof_indices_cache[e] = {
      offset_x + 2 * n0,     offset_x + 2 * n0 + 1,
      offset_x + 2 * n1,     offset_x + 2 * n1 + 1,
      offset_y + 2 * n0,     offset_y + 2 * n0 + 1,
      offset_y + 2 * n1,     offset_y + 2 * n1 + 1,
      offset_z + 2 * n0,     offset_z + 2 * n0 + 1,
      offset_z + 2 * n1,     offset_z + 2 * n1 + 1
    };
  }
}

void
EulerBeamInextensibleGGL::initialize_constant_element_matrices()
{
  local_bending_jacobian.setZero();

  for (size_t qi = 0; qi < xi_q.size(); ++qi) {
    const auto& ddH = quad_ddH[qi];
    const real_t wds = w_q[qi] * ds;

    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = 0; b < 4; ++b) {
        const real_t value = EI * ddH[a] * ddH[b] * wds;
        local_bending_jacobian(a, b) += value;
        local_bending_jacobian(a + 4, b + 4) += value;
        local_bending_jacobian(a + 8, b + 8) += value;
      }
    }
  }
}

void
EulerBeamInextensibleGGL::assemble_mass_matrix()
{
  using Tpl = Triplet<real_t>;

  std::vector<Tpl> triplets;
  triplets.reserve(elements * 3 * 16);

  // Each element contributes the same 4x4 Hermite mass block in x, y, and z.
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

    for (size_t qi = 0; qi < xi_q.size(); ++qi) {
      const auto& H = quad_H[qi];

      for (size_t a = 0; a < 4; ++a) {
        for (size_t b = 0; b < 4; ++b) {
          me[a][b] += mu * H[a] * H[b] * w_q[qi] * ds;
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

template<typename Scalar>
Matrix<Scalar, 16, 1>
EulerBeamInextensibleGGL::assemble_ggl_element(
  const Matrix<Scalar, 12, 1>& u_elem,
  const std::array<real_t, 2>& lambda_elem,
  const std::array<real_t, 4>& vx_elem,
  const std::array<real_t, 4>& vy_elem,
  const std::array<real_t, 4>& vz_elem,
  std::array<real_t, 3>        load) const
{
  // Readable standalone element residual for a uniform load. The optimized
  // assembly path uses the fused helper above.
  Matrix<Scalar, 16, 1> R = Matrix<Scalar, 16, 1>::Zero();

  for (size_t qi = 0; qi < xi_q.size(); ++qi) {
    const real_t xi = xi_q[qi];
    const real_t w = w_q[qi];

    const auto H = ELFF::FEM::CubicHermite<real_t>::values(xi, ds);
    const auto dH = ELFF::FEM::CubicHermite<real_t>::derivs(xi, ds);
    const auto ddH = ELFF::FEM::CubicHermite<real_t>::second_derivs(xi, ds);
    const auto M = ELFF::FEM::LinearShape<real_t>::values(xi);

    Scalar xp = 0, yp = 0, zp = 0;
    Scalar xpp = 0, ypp = 0, zpp = 0;
    for (size_t i = 0; i < 4; ++i) {
      xp += dH[i] * u_elem(i);
      yp += dH[i] * u_elem(i + 4);
      zp += dH[i] * u_elem(i + 8);
      xpp += ddH[i] * u_elem(i);
      ypp += ddH[i] * u_elem(i + 4);
      zpp += ddH[i] * u_elem(i + 8);
    }

    real_t vxp = 0, vyp = 0, vzp = 0;
    for (size_t i = 0; i < 4; ++i) {
      vxp += dH[i] * vx_elem[i];
      vyp += dH[i] * vy_elem[i];
      vzp += dH[i] * vz_elem[i];
    }

    const Scalar lambda_q = M[0] * lambda_elem[0] + M[1] * lambda_elem[1];

    for (size_t a = 0; a < 4; ++a) {
      R(a) += EI * xpp * ddH[a] * w * ds;
      R(a + 4) += EI * ypp * ddH[a] * w * ds;
      R(a + 8) += EI * zpp * ddH[a] * w * ds;

      const Scalar coeff = 2.0 * lambda_q * dH[a] * w * ds;
      R(a) -= coeff * xp;
      R(a + 4) -= coeff * yp;
      R(a + 8) -= coeff * zp;

      R(a) -= load[0] * H[a] * w * ds;
      R(a + 4) -= load[1] * H[a] * w * ds;
      R(a + 8) -= load[2] * H[a] * w * ds;
    }

    const Scalar g = xp * xp + yp * yp + zp * zp - Scalar(1.0);
    const Scalar Gv = xp * vxp + yp * vyp + zp * vzp;
    for (size_t a = 0; a < 2; ++a) {
      R(12 + a) += g * M[a] * w * ds;
      R(14 + a) += Gv * M[a] * w * ds;
    }
  }

  return R;
}

template<typename Scalar>
Matrix<Scalar, 16, 1>
EulerBeamInextensibleGGL::assemble_ggl_element(
  const Matrix<Scalar, 12, 1>&                u_elem,
  const std::array<real_t, 2>&                lambda_elem,
  const std::array<real_t, 4>&                vx_elem,
  const std::array<real_t, 4>&                vy_elem,
  const std::array<real_t, 4>&                vz_elem,
  const std::array<std::array<real_t, 3>, 2>& load_elem) const
{
  // Same residual as above, but with elementwise interpolated nodal loading.
  Matrix<Scalar, 16, 1> R = Matrix<Scalar, 16, 1>::Zero();

  for (size_t qi = 0; qi < xi_q.size(); ++qi) {
    const real_t xi = xi_q[qi];
    const real_t w = w_q[qi];

    const auto H = ELFF::FEM::CubicHermite<real_t>::values(xi, ds);
    const auto dH = ELFF::FEM::CubicHermite<real_t>::derivs(xi, ds);
    const auto ddH = ELFF::FEM::CubicHermite<real_t>::second_derivs(xi, ds);
    const auto M = ELFF::FEM::LinearShape<real_t>::values(xi);

    Scalar xp = 0, yp = 0, zp = 0;
    Scalar xpp = 0, ypp = 0, zpp = 0;
    for (size_t i = 0; i < 4; ++i) {
      xp += dH[i] * u_elem(i);
      yp += dH[i] * u_elem(i + 4);
      zp += dH[i] * u_elem(i + 8);
      xpp += ddH[i] * u_elem(i);
      ypp += ddH[i] * u_elem(i + 4);
      zpp += ddH[i] * u_elem(i + 8);
    }

    real_t vxp = 0, vyp = 0, vzp = 0;
    for (size_t i = 0; i < 4; ++i) {
      vxp += dH[i] * vx_elem[i];
      vyp += dH[i] * vy_elem[i];
      vzp += dH[i] * vz_elem[i];
    }

    const Scalar lambda_q = M[0] * lambda_elem[0] + M[1] * lambda_elem[1];
    const real_t fx = M[0] * load_elem[0][0] + M[1] * load_elem[1][0];
    const real_t fy = M[0] * load_elem[0][1] + M[1] * load_elem[1][1];
    const real_t fz = M[0] * load_elem[0][2] + M[1] * load_elem[1][2];

    for (size_t a = 0; a < 4; ++a) {
      R(a) += EI * xpp * ddH[a] * w * ds;
      R(a + 4) += EI * ypp * ddH[a] * w * ds;
      R(a + 8) += EI * zpp * ddH[a] * w * ds;

      const Scalar coeff = 2.0 * lambda_q * dH[a] * w * ds;
      R(a) -= coeff * xp;
      R(a + 4) -= coeff * yp;
      R(a + 8) -= coeff * zp;

      R(a) -= fx * H[a] * w * ds;
      R(a + 4) -= fy * H[a] * w * ds;
      R(a + 8) -= fz * H[a] * w * ds;
    }

    const Scalar g = xp * xp + yp * yp + zp * zp - Scalar(1.0);
    const Scalar Gv = xp * vxp + yp * vyp + zp * vzp;
    for (size_t a = 0; a < 2; ++a) {
      R(12 + a) += g * M[a] * w * ds;
      R(14 + a) += Gv * M[a] * w * ds;
    }
  }

  return R;
}

template<typename Scalar>
Matrix<Scalar, 12, 1>
EulerBeamInextensibleGGL::assemble_velocity_multiplier_residual(
  const Matrix<Scalar, 12, 1>& u_elem,
  const std::array<real_t, 2>& mu_elem) const
{
  // Isolate the -C(u)^T mu contribution used in the velocity block.
  Matrix<Scalar, 12, 1> Rv = Matrix<Scalar, 12, 1>::Zero();

  for (size_t qi = 0; qi < xi_q.size(); ++qi) {
    const real_t xi = xi_q[qi];
    const real_t w = w_q[qi];

    const auto dH = ELFF::FEM::CubicHermite<real_t>::derivs(xi, ds);
    const auto M = ELFF::FEM::LinearShape<real_t>::values(xi);

    Scalar xp = 0, yp = 0, zp = 0;
    for (size_t i = 0; i < 4; ++i) {
      xp += dH[i] * u_elem(i);
      yp += dH[i] * u_elem(i + 4);
      zp += dH[i] * u_elem(i + 8);
    }

    const real_t mu_q = M[0] * mu_elem[0] + M[1] * mu_elem[1];
    for (size_t a = 0; a < 4; ++a) {
      const Scalar coeff = mu_q * dH[a] * w * ds;
      Rv(a) -= coeff * xp;
      Rv(a + 4) -= coeff * yp;
      Rv(a + 8) -= coeff * zp;
    }
  }

  return Rv;
}

const std::array<size_t, 12>&
EulerBeamInextensibleGGL::get_element_disp_dof_indices(size_t e) const
{
  return element_disp_dof_indices_cache[e];
}

// -----------------------------------------------------------------------
// Local DOF extraction helpers
// -----------------------------------------------------------------------
std::array<real_t, 4>
EulerBeamInextensibleGGL::get_element_velocity_x(
  size_t e, const VectorXd& v) const
{
  const size_t n0 = e;
  const size_t n1 = e + 1;
  return { v(offset_x + 2 * n0),     v(offset_x + 2 * n0 + 1),
           v(offset_x + 2 * n1),     v(offset_x + 2 * n1 + 1) };
}

std::array<real_t, 4>
EulerBeamInextensibleGGL::get_element_velocity_y(
  size_t e, const VectorXd& v) const
{
  const size_t n0 = e;
  const size_t n1 = e + 1;
  return { v(offset_y + 2 * n0),     v(offset_y + 2 * n0 + 1),
           v(offset_y + 2 * n1),     v(offset_y + 2 * n1 + 1) };
}

std::array<real_t, 4>
EulerBeamInextensibleGGL::get_element_velocity_z(
  size_t e, const VectorXd& v) const
{
  const size_t n0 = e;
  const size_t n1 = e + 1;
  return { v(offset_z + 2 * n0),     v(offset_z + 2 * n0 + 1),
           v(offset_z + 2 * n1),     v(offset_z + 2 * n1 + 1) };
}

// -----------------------------------------------------------------------
// Full GGL assembly
// -----------------------------------------------------------------------
void
EulerBeamInextensibleGGL::assemble_ggl_system(
  const VectorXd&       u_cur,
  const VectorXd&       v_cur,
  const VectorXd&       lambda_cur,
  const VectorXd&       mu_cur,
  std::array<real_t, 3> load)
{
  using Tpl     = Triplet<real_t>;

  const VectorXd a_cur = compute_acceleration(u_cur);
  const VectorXd v_newmark = compute_velocity(u_cur);

  ggl_residual.setZero();
  ggl_residual.head(ndof) += mass * a_cur;
  ggl_residual.segment(offset_v, ndof) += mass * (v_cur - v_newmark);

  std::vector<char> constrained(ggl_dof, 0);
  VectorXd          constrained_value = VectorXd::Zero(ggl_dof);
  collect_ggl_constraints(constrained, constrained_value);

  std::vector<Tpl> triplets;
  triplets.reserve(3 * mass.nonZeros() + elements * 192 + ggl_dof / 8);

  auto add_jacobian_entry = [&](size_t row, size_t col, real_t value) {
    if (value == 0.0) {
      return;
    }
    if (constrained[row]) {
      return;
    }
    if (constrained[col]) {
      ggl_residual(row) -= value * constrained_value(col);
      return;
    }

    triplets.emplace_back(static_cast<int>(row),
                          static_cast<int>(col),
                          value);
  };

  for (int col = 0; col < mass.outerSize(); ++col) {
    for (SparseMatrix<real_t>::InnerIterator it(mass, col); it; ++it) {
      add_jacobian_entry(it.row(), it.col(), newmark_coeff_a * it.value());
      add_jacobian_entry(offset_v + it.row(),
                         it.col(),
                         -newmark_coeff_v * it.value());
      add_jacobian_entry(offset_v + it.row(),
                         offset_v + it.col(),
                         it.value());
    }
  }

  for (size_t e = 0; e < elements; ++e) {
    const auto& idx = get_element_disp_dof_indices(e);
    const std::array<real_t, 2> lambda_elem = { lambda_cur(e),
                                                lambda_cur(e + 1) };
    const std::array<real_t, 2> mu_elem = { mu_cur(e), mu_cur(e + 1) };
    const auto vx_elem = get_element_velocity_x(e, v_cur);
    const auto vy_elem = get_element_velocity_y(e, v_cur);
    const auto vz_elem = get_element_velocity_z(e, v_cur);
    Matrix<real_t, 12, 1> u_elem;
    for (int a = 0; a < 12; ++a) {
      u_elem(a) = u_cur(idx[a]);
    }

    const std::array<std::array<real_t, 3>, 2> load_elem = { load, load };
    GGLLinearizedElementData elem_data;

    assemble_ggl_element_linearized_data(
      EI, ds, u_elem, lambda_elem, mu_elem,
      vx_elem, vy_elem, vz_elem, load_elem,
      quad_H, quad_dH, quad_ddH, quad_M,
      local_bending_jacobian, elem_data);

    for (int a = 0; a < 12; ++a) {
      ggl_residual(idx[a]) += elem_data.residual(a);
      ggl_residual(offset_v + idx[a]) += elem_data.velocity_residual(a);

      for (int b = 0; b < 12; ++b) {
        add_jacobian_entry(idx[a],
                           idx[b],
                           elem_data.displacement_jacobian(a, b));
        add_jacobian_entry(offset_v + idx[a],
                           idx[b],
                           elem_data.velocity_displacement_jacobian(a, b));
      }
    }

    for (int a = 0; a < 2; ++a) {
      const size_t global_l = offset_lambda + e + a;
      const size_t global_mu = offset_mu + e + a;

      ggl_residual(global_l) += elem_data.residual(12 + a);
      ggl_residual(global_mu) += elem_data.residual(14 + a);

      for (int b = 0; b < 12; ++b) {
        const real_t dRlambda =
          elem_data.lambda_displacement_jacobian(a, b);
        add_jacobian_entry(global_l, idx[b], dRlambda);
        add_jacobian_entry(idx[b], global_l, -dRlambda);
        add_jacobian_entry(global_mu,
                           idx[b],
                           elem_data.mu_displacement_jacobian(a, b));

        const real_t dRmu_dv = elem_data.mu_velocity_jacobian(a, b);
        add_jacobian_entry(global_mu, offset_v + idx[b], dRmu_dv);
        add_jacobian_entry(offset_v + idx[b], global_mu, -dRmu_dv);
      }
    }
  }

  for (size_t bi = 0; bi < 2; ++bi) {
    const size_t          ni     = (boundary_conditions.end[bi] == left)
                                   ? 0 : nodes - 1;
    const EulerBeamBCType bctype = boundary_conditions.type[bi];
    const EulerBeamBCVals bcvals = boundary_conditions.vals[bi];

    if (bctype == point_force_bc) {
      ggl_residual(offset_x + 2 * ni) -= bcvals.force[0];
      ggl_residual(offset_y + 2 * ni) -= bcvals.force[1];
      ggl_residual(offset_z + 2 * ni) -= bcvals.force[2];
    } else if (bctype == point_torque_bc) {
      ggl_residual(offset_x + 2 * ni + 1) -= bcvals.torque[0];
      ggl_residual(offset_y + 2 * ni + 1) -= bcvals.torque[1];
      ggl_residual(offset_z + 2 * ni + 1) -= bcvals.torque[2];
    }
  }

  for (size_t dof = 0; dof < ggl_dof; ++dof) {
    if (!constrained[dof]) {
      continue;
    }

    triplets.emplace_back(static_cast<int>(dof),
                          static_cast<int>(dof),
                          1.0);
    if (dof < offset_v) {
      ggl_residual(dof) = u_cur(dof) - constrained_value(dof);
    } else if (dof < offset_lambda) {
      ggl_residual(dof) = v_cur(dof - offset_v) - constrained_value(dof);
    } else {
      ggl_residual(dof) = -constrained_value(dof);
    }
  }

  ggl_jacobian.resize(ggl_dof, ggl_dof);
  ggl_jacobian.setFromTriplets(triplets.begin(), triplets.end());
  ggl_jacobian.makeCompressed();
}

void
EulerBeamInextensibleGGL::assemble_ggl_system(
  const VectorXd&                            u_cur,
  const VectorXd&                            v_cur,
  const VectorXd&                            lambda_cur,
  const VectorXd&                            mu_cur,
  const std::vector<std::array<real_t, 3>>& load)
{
  using Tpl     = Triplet<real_t>;

  const VectorXd a_cur = compute_acceleration(u_cur);
  const VectorXd v_newmark = compute_velocity(u_cur);

  ggl_residual.setZero();
  ggl_residual.head(ndof) += mass * a_cur;
  ggl_residual.segment(offset_v, ndof) += mass * (v_cur - v_newmark);

  std::vector<char> constrained(ggl_dof, 0);
  VectorXd          constrained_value = VectorXd::Zero(ggl_dof);
  collect_ggl_constraints(constrained, constrained_value);

  std::vector<Tpl> triplets;
  triplets.reserve(3 * mass.nonZeros() + elements * 192 + ggl_dof / 8);

  auto add_jacobian_entry = [&](size_t row, size_t col, real_t value) {
    if (value == 0.0) {
      return;
    }
    if (constrained[row]) {
      return;
    }
    if (constrained[col]) {
      ggl_residual(row) -= value * constrained_value(col);
      return;
    }

    triplets.emplace_back(static_cast<int>(row),
                          static_cast<int>(col),
                          value);
  };

  for (int col = 0; col < mass.outerSize(); ++col) {
    for (SparseMatrix<real_t>::InnerIterator it(mass, col); it; ++it) {
      add_jacobian_entry(it.row(), it.col(), newmark_coeff_a * it.value());
      add_jacobian_entry(offset_v + it.row(),
                         it.col(),
                         -newmark_coeff_v * it.value());
      add_jacobian_entry(offset_v + it.row(),
                         offset_v + it.col(),
                         it.value());
    }
  }

  for (size_t e = 0; e < elements; ++e) {
    const auto& idx = get_element_disp_dof_indices(e);
    const std::array<real_t, 2> lambda_elem = { lambda_cur(e),
                                                lambda_cur(e + 1) };
    const std::array<real_t, 2> mu_elem = { mu_cur(e), mu_cur(e + 1) };
    const auto vx_elem = get_element_velocity_x(e, v_cur);
    const auto vy_elem = get_element_velocity_y(e, v_cur);
    const auto vz_elem = get_element_velocity_z(e, v_cur);
    const std::array<std::array<real_t, 3>, 2> load_elem = { load[e],
                                                             load[e + 1] };
    Matrix<real_t, 12, 1> u_elem;
    for (int a = 0; a < 12; ++a) {
      u_elem(a) = u_cur(idx[a]);
    }

    GGLLinearizedElementData elem_data;

    assemble_ggl_element_linearized_data(
      EI, ds, u_elem, lambda_elem, mu_elem,
      vx_elem, vy_elem, vz_elem, load_elem,
      quad_H, quad_dH, quad_ddH, quad_M,
      local_bending_jacobian, elem_data);

    for (int a = 0; a < 12; ++a) {
      ggl_residual(idx[a]) += elem_data.residual(a);
      ggl_residual(offset_v + idx[a]) += elem_data.velocity_residual(a);

      for (int b = 0; b < 12; ++b) {
        add_jacobian_entry(idx[a],
                           idx[b],
                           elem_data.displacement_jacobian(a, b));
        add_jacobian_entry(offset_v + idx[a],
                           idx[b],
                           elem_data.velocity_displacement_jacobian(a, b));
      }
    }

    for (int a = 0; a < 2; ++a) {
      const size_t global_l = offset_lambda + e + a;
      const size_t global_mu = offset_mu + e + a;

      ggl_residual(global_l) += elem_data.residual(12 + a);
      ggl_residual(global_mu) += elem_data.residual(14 + a);

      for (int b = 0; b < 12; ++b) {
        const real_t dRlambda =
          elem_data.lambda_displacement_jacobian(a, b);
        add_jacobian_entry(global_l, idx[b], dRlambda);
        add_jacobian_entry(idx[b], global_l, -dRlambda);
        add_jacobian_entry(global_mu,
                           idx[b],
                           elem_data.mu_displacement_jacobian(a, b));

        const real_t dRmu_dv = elem_data.mu_velocity_jacobian(a, b);
        add_jacobian_entry(global_mu, offset_v + idx[b], dRmu_dv);
        add_jacobian_entry(offset_v + idx[b], global_mu, -dRmu_dv);
      }
    }
  }

  for (size_t bi = 0; bi < 2; ++bi) {
    const size_t          ni     = (boundary_conditions.end[bi] == left)
                                   ? 0 : nodes - 1;
    const EulerBeamBCType bctype = boundary_conditions.type[bi];
    const EulerBeamBCVals bcvals = boundary_conditions.vals[bi];

    if (bctype == point_force_bc) {
      ggl_residual(offset_x + 2 * ni) -= bcvals.force[0];
      ggl_residual(offset_y + 2 * ni) -= bcvals.force[1];
      ggl_residual(offset_z + 2 * ni) -= bcvals.force[2];
    } else if (bctype == point_torque_bc) {
      ggl_residual(offset_x + 2 * ni + 1) -= bcvals.torque[0];
      ggl_residual(offset_y + 2 * ni + 1) -= bcvals.torque[1];
      ggl_residual(offset_z + 2 * ni + 1) -= bcvals.torque[2];
    }
  }

  for (size_t dof = 0; dof < ggl_dof; ++dof) {
    if (!constrained[dof]) {
      continue;
    }

    triplets.emplace_back(static_cast<int>(dof),
                          static_cast<int>(dof),
                          1.0);
    if (dof < offset_v) {
      ggl_residual(dof) = u_cur(dof) - constrained_value(dof);
    } else if (dof < offset_lambda) {
      ggl_residual(dof) = v_cur(dof - offset_v) - constrained_value(dof);
    } else {
      ggl_residual(dof) = -constrained_value(dof);
    }
  }

  ggl_jacobian.resize(ggl_dof, ggl_dof);
  ggl_jacobian.setFromTriplets(triplets.begin(), triplets.end());
  ggl_jacobian.makeCompressed();
}

void
EulerBeamInextensibleGGL::collect_ggl_constraints(
  std::vector<char>& constrained,
  VectorXd&          constrained_value) const
{
  constrained.assign(ggl_dof, 0);
  constrained_value.setZero(ggl_dof);

  auto add_constraint = [&](size_t dof, real_t value) {
    constrained[dof] = 1;
    constrained_value(dof) = value;
  };

  for (size_t bi = 0; bi < 2; ++bi) {
    const size_t          ni     = (boundary_conditions.end[bi] == left)
                                   ? 0 : nodes - 1;
    const EulerBeamBCType bctype = boundary_conditions.type[bi];
    const EulerBeamBCVals bcvals = boundary_conditions.vals[bi];

    std::vector<std::pair<size_t, real_t>> disp_bcs;

    switch (bctype) {
      case simple_bc:
        disp_bcs.push_back({ offset_x + 2 * ni, bcvals.position[0] });
        disp_bcs.push_back({ offset_y + 2 * ni, bcvals.position[1] });
        disp_bcs.push_back({ offset_z + 2 * ni, bcvals.position[2] });
        break;
      case clamped_bc:
        disp_bcs.push_back({ offset_x + 2 * ni,     bcvals.position[0] });
        disp_bcs.push_back({ offset_x + 2 * ni + 1, bcvals.slope[0] });
        disp_bcs.push_back({ offset_y + 2 * ni,     bcvals.position[1] });
        disp_bcs.push_back({ offset_y + 2 * ni + 1, bcvals.slope[1] });
        disp_bcs.push_back({ offset_z + 2 * ni,     bcvals.position[2] });
        disp_bcs.push_back({ offset_z + 2 * ni + 1, bcvals.slope[2] });
        break;
      default:
        break;
    }

    for (const auto& [dof, value] : disp_bcs) {
      add_constraint(dof, value);
      add_constraint(offset_v + dof, 0.0);
    }

    if (bctype == free_bc) {
      add_constraint(offset_lambda + ni, 0.0);
      add_constraint(offset_mu + ni, 0.0);
    }
  }
}

// -----------------------------------------------------------------------
// State commit and output
// -----------------------------------------------------------------------
void
EulerBeamInextensibleGGL::update_newmark_state(
  const VectorXd& u_new,
  const VectorXd& v_new)
{
  // Acceleration remains derived from displacement, while velocity is solved
  // as an independent state variable.
  a_prev = compute_acceleration(u_new);
  v_prev = v_new;
  u = u_new;
}

void
EulerBeamInextensibleGGL::apply_dynamic_state_boundary_conditions()
{
  for (size_t bi = 0; bi < 2; ++bi) {
    const EulerBeamBCType bctype = boundary_conditions.type[bi];
    if (bctype != simple_bc && bctype != clamped_bc) {
      continue;
    }

    const size_t ni = (boundary_conditions.end[bi] == left) ? 0 : nodes - 1;
    const size_t ix = offset_x + 2 * ni;
    const size_t iy = offset_y + 2 * ni;
    const size_t iz = offset_z + 2 * ni;

    v_prev(ix) = 0.0;
    v_prev(iy) = 0.0;
    v_prev(iz) = 0.0;
    a_prev(ix) = 0.0;
    a_prev(iy) = 0.0;
    a_prev(iz) = 0.0;

    if (bctype == clamped_bc) {
      const size_t sx = offset_x + 2 * ni + 1;
      const size_t sy = offset_y + 2 * ni + 1;
      const size_t sz = offset_z + 2 * ni + 1;
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
EulerBeamInextensibleGGL::update_mesh()
{
  std::vector<std::array<real_t, 3>>& centerline = mesh.get_centerline();
  std::vector<std::array<real_t, 3>>& slope = mesh.get_slope();
  std::vector<std::array<real_t, 3>>& vel = mesh.get_centerline_velocity();
  for (size_t i = 0; i < nodes; ++i) {
    centerline[i][0] = u(offset_x + 2 * i);
    centerline[i][1] = u(offset_y + 2 * i);
    centerline[i][2] = u(offset_z + 2 * i);
    slope[i][0] = u(offset_x + 2 * i + 1);
    slope[i][1] = u(offset_y + 2 * i + 1);
    slope[i][2] = u(offset_z + 2 * i + 1);
    vel[i][0] = v_prev(offset_x + 2 * i);
    vel[i][1] = v_prev(offset_y + 2 * i);
    vel[i][2] = v_prev(offset_z + 2 * i);
  }
}

} // namespace Models
} // namespace ELFF
