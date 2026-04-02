#pragma once

#include "elff/models/beam/EulerBeam.hpp"
#include "elff/fem/Shapes.hpp"

#include <cstdio>
#include <stdlib.h> 

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky> // for SimplicialLLT
#include <Eigen/SparseLU>
#include <unsupported/Eigen/AutoDiff>

using namespace Eigen;

namespace ELFF {
namespace Models {

/**
 * @brief Inextensible Euler beam solved with a sparse
 * method-of-multipliers formulation.
 *
 * The same finite-element formulation is used for both static equilibrium
 * solves and dynamic Newmark time stepping. The dynamic overloads are enabled
 * when the beam is constructed with a positive mass density `mu`.
 */
class EulerBeamInextensibleMoM : public EulerBeam
{
public:
  /**
   * @brief Constructs a static inextensible sparse beam model.
   *
   * @param length Beam length
   * @param EI Flexural rigidity
   * @param nodes Number of discretization nodes
   * @param bcs Boundary conditions at the beam ends
   * @param r_penalty Penalty parameter used in the inextensibility constraint
   */
  EulerBeamInextensibleMoM(real_t length,
                           real_t EI,
                           size_t nodes,
                           EulerBeam::EulerBeamBCs bcs,
                           real_t r_penalty);

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
  EulerBeamInextensibleMoM(real_t length,
                           real_t EI,
                           real_t mu,
                           size_t nodes,
                           EulerBeam::EulerBeamBCs bcs,
                           real_t r_penalty);

  /**
   * @brief Destroys the sparse beam model and releases any owned resources.
   */
  ~EulerBeamInextensibleMoM();

  /**
   * @brief Solves the static beam problem with no external load.
   */
  void solve() override;

  /**
   * @brief Solves the static beam problem under uniform loading.
   *
   * @param load Uniform distributed load vector
   */
  void solve(std::array<real_t, 3> load) override;

  /**
   * @brief Solves the static beam problem under nodal loading.
   *
   * @param load Load vector specified at the mesh nodes
   */
  void solve(std::vector<std::array<real_t, 3>> load) override;

  /**
   * @brief Advances the beam one time step under uniform loading.
   *
   * @param dt Time-step size
   * @param load Uniform distributed load vector
   */
  void solve(real_t dt, std::array<real_t, 3> load) override;

  /**
   * @brief Advances the beam one time step under nodal loading.
   *
   * @param dt Time-step size
   * @param load Load vector specified at the mesh nodes
   */
  void solve(real_t dt, std::vector<std::array<real_t, 3>> load) override;

  /**
   * @brief Initializes the solver state from a supplied beam mesh.
   *
   * @param bmesh Beam mesh containing the initial geometry
   */
  virtual void apply_initial_condition(EulerBeamMesh& bmesh) override;

  /**
   * @brief Applies the model's default initial condition.
   */
  virtual void apply_initial_condition();

  /**
   * @brief Returns the current Lagrange multiplier state.
   */
  const VectorXd& get_lambda() const { return lambda; }

  /**
   * @brief Replaces the current Lagrange multiplier state.
   *
   * @param lambda_in Multiplier vector to use as the current iterate
   */
  void set_lambda(const VectorXd& lambda_in)
  {
    ELFF_ASSERT(lambda_in.size() == lambda.size(),
                "Size of lambda vector must equal ndof_l.\n");
    lambda = lambda_in;
    apply_lambda_boundary_conditions();
  }

protected:
  /**
   * @brief Spatial dimension of the beam formulation.
   */
  size_t dimension;
  /**
   * @brief Element and node counts used by the discretization.
   */
  size_t elements, nodes;
  /**
   * @brief Uniform spacing between nodes along the beam centerline.
   */
  real_t ds;
  /**
   * @brief Degrees of freedom in each Cartesian block and in the constraint
   * block.
   */
  size_t ndof_x, ndof_y, ndof_z, ndof_l;
  /**
   * @brief Offsets into the global state vector for each degree-of-freedom
   * block.
   */
  size_t offset_x, offset_y, offset_z, offset_l;
  /**
   * @brief Total number of algebraic degrees of freedom.
   */
  size_t ndof;
  /**
   * @brief Penalty parameter used to enforce inextensibility.
   */
  real_t r_penalty;
  /**
   * @brief Inner and outer nonlinear-iteration limits.
   */
  size_t max_iter_inner, max_iter_outer, min_iter_inner, min_iter_outer;
  /**
   * @brief Inner and outer nonlinear solver tolerances.
   */
  real_t tol_linear, tol_primal, tol_constraint;

  /**
   * @brief Current residual vector and Lagrange multiplier state.
   */
  VectorXd residual, lambda;
  /**
   * @brief Sparse global tangent matrix.
   */
  SparseMatrix<real_t> jacobian;
  SparseMatrix<real_t> mass;
  SparseMatrix<real_t> bending_matrix;
  /**
   * @brief Global state vector for beam displacements and slopes.
   */
  VectorXd u;
  /**
   * @brief State from the previous dynamic time step.
   */
  VectorXd u_prev, v_prev, a_prev;
  std::array<real_t, 3> load_prev;
  std::vector<std::array<real_t, 3>> nodal_load_prev;
  bool have_prev_uniform_load;
  bool have_prev_nodal_load;

  /**
   * @brief Assembles the residual for uniform loading.
   *
   * @param load Uniform distributed load vector
   */
  void assemble_residual(std::array<real_t, 3> load);

  /**
   * @brief Assembles the residual for nodal loading.
   *
   * @param load Load vector specified at the mesh nodes
   */
  void assemble_residual(std::vector<std::array<real_t, 3>> load);

  /**
   * @brief Returns the global degree-of-freedom indices for a two-node beam
   * element.
   *
   * @param e Element index
   * @return Global indices for the element displacement and slope unknowns
   */
  std::array<size_t, 12> get_element_dof_indices(size_t e) const;

  /**
   * @brief Extracts the local element state from the global state vector.
   *
   * @param idx Global degree-of-freedom indices for the element
   * @return Local 12-entry state vector for the element
   */
  Matrix<real_t, 12, 1> get_element_state(const std::array<size_t, 12>& idx) const;

  /**
   * @brief Extracts the local Lagrange-multiplier values for an element.
   *
   * @param e Element index
   * @return Multiplier values associated with the element endpoints
   */
  std::array<real_t, 2> get_element_lambda(size_t e) const;

  template<typename T>
  /**
   * @brief Assembles a local static element residual for uniform loading.
   *
   * @tparam T Scalar type used for residual assembly
   * @param u_elem Local element state vector
   * @param lambda_elem Element-end multiplier values
   * @param load Uniform distributed load vector
   * @return Local residual contribution for the element
   */
  Matrix<T, 12, 1> assemble_element_residual_template(
    const Matrix<T, 12, 1>& u_elem,
    const std::array<real_t, 2>& lambda_elem,
    const std::array<real_t, 3>& load,
    real_t                       bending_scale = 1.0) const
  {
    Matrix<T, 12, 1> residual = Matrix<T, 12, 1>::Zero();

    const Matrix<T, 4, 1> ux = u_elem.template segment<4>(0);
    const Matrix<T, 4, 1> uy = u_elem.template segment<4>(4);
    const Matrix<T, 4, 1> uz = u_elem.template segment<4>(8);

    const real_t xi_q[] = { 0.1127016654, 0.5, 0.8872983346 };
    const real_t w_q[] = { 0.2777777778, 0.4444444444, 0.2777777778 };

    for (size_t qi = 0; qi < 3; ++qi) {
      const real_t xi = xi_q[qi];
      const real_t w = w_q[qi];

      const auto H = ELFF::FEM::CubicHermite<real_t>::values(xi, ds);
      const auto dH = ELFF::FEM::CubicHermite<real_t>::derivs(xi, ds);
      const auto ddH = ELFF::FEM::CubicHermite<real_t>::second_derivs(xi, ds);
      const auto M = ELFF::FEM::LinearShape<real_t>::values(xi);

      T xp = 0, xpp = 0;
      T yp = 0, ypp = 0;
      T zp = 0, zpp = 0;

      for (size_t i = 0; i < 4; ++i) {
        xp += dH[i] * ux(i);
        yp += dH[i] * uy(i);
        zp += dH[i] * uz(i);
        xpp += ddH[i] * ux(i);
        ypp += ddH[i] * uy(i);
        zpp += ddH[i] * uz(i);
      }

      const T l = M[0] * lambda_elem[0] + M[1] * lambda_elem[1];
      const T S = xp * xp + yp * yp + zp * zp - 1.0;

      for (size_t a = 0; a < 4; ++a) {
        residual(a) += bending_scale * EI * xpp * ddH[a] * w * ds;
        residual(4 + a) += bending_scale * EI * ypp * ddH[a] * w * ds;
        residual(8 + a) += bending_scale * EI * zpp * ddH[a] * w * ds;

        residual(a) -= load[0] * H[a] * w * ds;
        residual(4 + a) -= load[1] * H[a] * w * ds;
        residual(8 + a) -= load[2] * H[a] * w * ds;

        const T coeff = 2 * (l + r_penalty * S) * dH[a] * w * ds;
        residual(a) += xp * coeff;
        residual(4 + a) += yp * coeff;
        residual(8 + a) += zp * coeff;
      }
    }

    return residual;
  }

  template<typename T>
  /**
   * @brief Assembles a local static element residual for nodal loading.
   *
   * @tparam T Scalar type used for residual assembly
   * @param u_elem Local element state vector
   * @param lambda_elem Element-end multiplier values
   * @param load_elem Nodal load values at the two element endpoints
   * @return Local residual contribution for the element
   */
  Matrix<T, 12, 1> assemble_element_residual_template(
    const Matrix<T, 12, 1>& u_elem,
    const std::array<real_t, 2>& lambda_elem,
    const std::array<std::array<real_t, 3>, 2>& load_elem,
    real_t                                      bending_scale = 1.0) const
  {
    Matrix<T, 12, 1> residual = Matrix<T, 12, 1>::Zero();

    const Matrix<T, 4, 1> ux = u_elem.template segment<4>(0);
    const Matrix<T, 4, 1> uy = u_elem.template segment<4>(4);
    const Matrix<T, 4, 1> uz = u_elem.template segment<4>(8);

    const real_t xi_q[] = { 0.1127016654, 0.5, 0.8872983346 };
    const real_t w_q[] = { 0.2777777778, 0.4444444444, 0.2777777778 };

    for (size_t qi = 0; qi < 3; ++qi) {
      const real_t xi = xi_q[qi];
      const real_t w = w_q[qi];

      const auto H = ELFF::FEM::CubicHermite<real_t>::values(xi, ds);
      const auto dH = ELFF::FEM::CubicHermite<real_t>::derivs(xi, ds);
      const auto ddH = ELFF::FEM::CubicHermite<real_t>::second_derivs(xi, ds);
      const auto M = ELFF::FEM::LinearShape<real_t>::values(xi);

      T xp = 0, xpp = 0;
      T yp = 0, ypp = 0;
      T zp = 0, zpp = 0;

      for (size_t i = 0; i < 4; ++i) {
        xp += dH[i] * ux(i);
        yp += dH[i] * uy(i);
        zp += dH[i] * uz(i);
        xpp += ddH[i] * ux(i);
        ypp += ddH[i] * uy(i);
        zpp += ddH[i] * uz(i);
      }

      const T l = M[0] * lambda_elem[0] + M[1] * lambda_elem[1];
      const real_t fx = M[0] * load_elem[0][0] + M[1] * load_elem[1][0];
      const real_t fy = M[0] * load_elem[0][1] + M[1] * load_elem[1][1];
      const real_t fz = M[0] * load_elem[0][2] + M[1] * load_elem[1][2];
      const T S = xp * xp + yp * yp + zp * zp - 1.0;

      for (size_t a = 0; a < 4; ++a) {
        residual(a) += bending_scale * EI * xpp * ddH[a] * w * ds;
        residual(4 + a) += bending_scale * EI * ypp * ddH[a] * w * ds;
        residual(8 + a) += bending_scale * EI * zpp * ddH[a] * w * ds;

        residual(a) -= fx * H[a] * w * ds;
        residual(4 + a) -= fy * H[a] * w * ds;
        residual(8 + a) -= fz * H[a] * w * ds;

        const T coeff = 2 * (l + r_penalty * S) * dH[a] * w * ds;
        residual(a) += xp * coeff;
        residual(4 + a) += yp * coeff;
        residual(8 + a) += zp * coeff;
      }
    }

    return residual;
  }

  /**
   * @brief Updates the Lagrange multiplier iterate.
   *
   * @param omega Relaxation parameter for the multiplier update
   * @return Norm of the multiplier correction
   */
  real_t update_lambda(real_t omega = 1.);

  /**
   * @brief Updates the multiplier iterate for dynamic Newmark solves.
   *
   * Uses the same defect projection as the static update but with an added
   * pseudo-time regularization tied to the step size so the dual correction is
   * less aggressive during transient solves.
   *
   * @param dt Time-step size
   * @param omega Relaxation parameter for the multiplier update
   * @return Norm of the multiplier correction
   */
  real_t update_lambda_newmark(real_t dt, real_t omega = 0.1);


  /**
   * @brief Computes the L2 norm of the inextensibility defect.
   *
   * This evaluates \f$S = |\partial_s \mathbf{r}|^2 - 1\f$ at the element
   * quadrature points and returns \f$\sqrt{\int S^2\,ds}\f$.
   */
  real_t compute_inextensibility_error_l2() const;

  /**
   * @brief Computes the maximum quadrature-point inextensibility defect.
   *
   * This evaluates \f$|S|\f$ at the element quadrature points and returns the
   * maximum observed value.
   */
  real_t compute_inextensibility_error_linf() const;

  /**
   * @brief Updates the beam mesh from the current solution state.
   */
  void update_mesh();

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
   * @brief Initializes the acceleration state at the start of a simulation.
   *
   * @param load Uniform distributed load vector
   */
  void initialize_acceleration(std::array<real_t, 3> load);

  /**
   * @brief Initializes the acceleration state at the start of a simulation.
   *
   * @param load Load vector specified at the mesh nodes
   */
  void initialize_acceleration(std::vector<std::array<real_t, 3>> load);
  void assemble_mass_matrix();
  void assemble_bending_matrix();

  /**
   * @brief Assembles the nonlinear system for uniform loading.
   *
   * @param load Uniform distributed load vector
   */
  void assemble_system(std::array<real_t, 3> load,
                       real_t                bending_scale = 1.0);

  /**
   * @brief Assembles the nonlinear system for nodal loading.
   *
   * @param load Load vector specified at the mesh nodes
   */
  void assemble_system(std::vector<std::array<real_t, 3>> load,
                       real_t                              bending_scale = 1.0);

  /**
   * @brief Apply boundary conditions to the residual and jacobian
   *
   * Modifies the residual and jacobian according to boundary conditions:
   * - free_bc: No constraints
   * - simple_bc: Position constraints only
   * - clamped_bc: Position and slope constraints
   */

  void apply_boundary_conditions();

  /**
   * @brief Applies boundary conditions directly to the multiplier field.
   *
   * For now, free ends enforce zero endpoint multiplier values while leaving
   * any enriched interior multiplier DOFs untouched.
   */
  void apply_lambda_boundary_conditions();

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
  void add_midpoint_bending_terms();
  void add_averaged_uniform_load(std::array<real_t, 3>& load) const;
  void add_averaged_nodal_load(std::vector<std::array<real_t, 3>>& load) const;

  /**
   * @brief Adds the Newmark inertial residual and tangent contributions.
   *
   * @param dt Time-step size
   * @param beta Newmark beta parameter
   */
  void add_newmark_inertial_terms(real_t dt, real_t beta);

  /**
   * @brief Updates the stored Newmark velocity and acceleration from the
   * converged displacement.
   *
   * @param dt Time-step size
   * @param beta Newmark beta parameter
   * @param gamma Newmark gamma parameter
   * @param v_old Velocity state from the previous time step
   * @param a_old Acceleration state from the previous time step
   */
  void update_newmark_state_from_displacement(real_t dt,
                                              real_t beta,
                                              real_t gamma,
                                              const VectorXd& v_old,
                                              const VectorXd& a_old);

  /**
   * @brief Enforces translational velocity and acceleration boundary
   * conditions on the stored Newmark state.
   */
  void apply_dynamic_state_boundary_conditions();

  /**
   * @brief Applies a lumped-mass velocity projection driven by the
   * velocity-level inextensibility defect.
   */
  void project_velocity_onto_constraint_manifold();

  /**
   * @brief Assemble the residual vector for the nonlinear system
   *
   * @tparam T Scalar type (real_t or autodiff)
   * @param u Current solution vector
   * @param load Current load on the beam
   * @return Assembled residual vector containing:
   *         - Bending energy terms
   *         - External load contributions
   *         - Inextensibility constraints
   *
   * Uses Gauss quadrature with 3 points for numerical integration.
   */
  template<typename T>
  Matrix<T, Dynamic, 1> assemble_residual_template(
    const Matrix<T, Dynamic, 1>& u,
    const std::array<real_t, 3> load) const
  {
    real_t xi_q[] = { 0.1127016654, 0.5, 0.8872983346 };
    real_t w_q[] = { 0.2777777778, 0.4444444444, 0.2777777778 };

    Matrix<T, Dynamic, 1> residual = Matrix<T, Dynamic, 1>::Zero(ndof);

    for (size_t e = 0; e < elements; ++e) {
      std::vector<size_t> elem_nodes = { e, e + 1 };
      std::vector<size_t> idx_x = { offset_x + 2 * elem_nodes[0],
                                    offset_x + 2 * elem_nodes[0] + 1,
                                    offset_x + 2 * elem_nodes[1],
                                    offset_x + 2 * elem_nodes[1] + 1 };
      std::vector<size_t> idx_y = { offset_y + 2 * elem_nodes[0],
                                    offset_y + 2 * elem_nodes[0] + 1,
                                    offset_y + 2 * elem_nodes[1],
                                    offset_y + 2 * elem_nodes[1] + 1 };
      std::vector<size_t> idx_z = { offset_z + 2 * elem_nodes[0],
                                    offset_z + 2 * elem_nodes[0] + 1,
                                    offset_z + 2 * elem_nodes[1],
                                    offset_z + 2 * elem_nodes[1] + 1 };
      std::vector<size_t> idx_l = { elem_nodes[0], elem_nodes[1] };

      std::array<T, 4> ux = {
        u[idx_x[0]], u[idx_x[1]], u[idx_x[2]], u[idx_x[3]]
      };
      std::array<T, 4> uy = {
        u[idx_y[0]], u[idx_y[1]], u[idx_y[2]], u[idx_y[3]]
      };
      std::array<T, 4> uz = {
        u[idx_z[0]], u[idx_z[1]], u[idx_z[2]], u[idx_z[3]]
      };
      std::array<T, 2> ul = { lambda[idx_l[0]], lambda[idx_l[1]] };

      std::vector<T> R_loc_x(4, 0);
      std::vector<T> R_loc_y(4, 0);
      std::vector<T> R_loc_z(4, 0);
      // std::vector<T> R_loc_l(2, 0);

      for (size_t qi = 0; qi < 3; ++qi) {
        real_t xi = xi_q[qi];
        real_t w = w_q[qi];

        auto H = ELFF::FEM::CubicHermite<real_t>::values(xi, ds);
        auto dH = ELFF::FEM::CubicHermite<real_t>::derivs(xi, ds);
        auto ddH = ELFF::FEM::CubicHermite<real_t>::second_derivs(xi, ds);
        auto M = ELFF::FEM::LinearShape<real_t>::values(xi);

        T x = 0, xp = 0, xpp = 0;
        T y = 0, yp = 0, ypp = 0;
        T z = 0, zp = 0, zpp = 0;
        for (size_t i = 0; i < 4; i++) {
          x += H[i] * ux[i];
          xp += dH[i] * ux[i];
          xpp += ddH[i] * ux[i];
          y += H[i] * uy[i];
          yp += dH[i] * uy[i];
          ypp += ddH[i] * uy[i];
          z += H[i] * uz[i];
          zp += dH[i] * uz[i];
          zpp += ddH[i] * uz[i];
        }

        T l = 0;
        for (size_t i = 0; i < 2; i++) {
          l += M[i] * ul[i];
        }

        T S = xp * xp + yp * yp + zp * zp - 1.0;

        for (size_t a = 0; a < 4; ++a) {
          // Bending Energy Contributions
          R_loc_x[a] += EI * xpp * ddH[a] * w * ds;
          R_loc_y[a] += EI * ypp * ddH[a] * w * ds;
          R_loc_z[a] += EI * zpp * ddH[a] * w * ds;
          // External load contribution in y
          R_loc_x[a] -= load[0] * H[a] * w * ds;
          R_loc_y[a] -= load[1] * H[a] * w * ds;
          R_loc_z[a] -= load[2] * H[a] * w * ds;
          // Constraint contributions
          R_loc_x[a] += 2 * (l + r_penalty * S) * xp * dH[a] * w * ds;
          R_loc_y[a] += 2 * (l + r_penalty * S) * yp * dH[a] * w * ds;
          R_loc_z[a] += 2 * (l + r_penalty * S) * zp * dH[a] * w * ds;
        }

        // for (size_t a = 0; a < 2; ++a) {
        //   // Variation w.r.t. lambda
        //   R_loc_l[a] += S * M[a] * w * ds;
        // }
      }

      // Scatter local residual contribution to global residual
      for (size_t i = 0; i < 4; ++i) {
        residual[idx_x[i]] += R_loc_x[i];
        residual[idx_y[i]] += R_loc_y[i];
        residual[idx_z[i]] += R_loc_z[i];
      }

      // for (size_t i = 0; i < 2; ++i) {
      //   residual[idx_l[i]] += R_loc_l[i];
      // }
    }
    return residual;
  }

  template<typename T>
  Matrix<T, Dynamic, 1> assemble_residual_template(
    const Matrix<T, Dynamic, 1>& u,
    const std::vector<std::array<real_t, 3>> load) const
  {
    ELFF_ASSERT(load.size() == nodes,
                "Nodes does not match load vector size.\n");

    real_t xi_q[] = { 0.1127016654, 0.5, 0.8872983346 };
    real_t w_q[] = { 0.2777777778, 0.4444444444, 0.2777777778 };

    Matrix<T, Dynamic, 1> residual = Matrix<T, Dynamic, 1>::Zero(ndof);

    for (size_t e = 0; e < elements; ++e) {
      std::vector<size_t> elem_nodes = { e, e + 1 };
      std::vector<size_t> idx_x = { offset_x + 2 * elem_nodes[0],
                                    offset_x + 2 * elem_nodes[0] + 1,
                                    offset_x + 2 * elem_nodes[1],
                                    offset_x + 2 * elem_nodes[1] + 1 };
      std::vector<size_t> idx_y = { offset_y + 2 * elem_nodes[0],
                                    offset_y + 2 * elem_nodes[0] + 1,
                                    offset_y + 2 * elem_nodes[1],
                                    offset_y + 2 * elem_nodes[1] + 1 };
      std::vector<size_t> idx_z = { offset_z + 2 * elem_nodes[0],
                                    offset_z + 2 * elem_nodes[0] + 1,
                                    offset_z + 2 * elem_nodes[1],
                                    offset_z + 2 * elem_nodes[1] + 1 };
      std::vector<size_t> idx_load = { elem_nodes[0], elem_nodes[1] };
      std::vector<size_t> idx_l = { elem_nodes[0], elem_nodes[1] };

      std::array<T, 4> ux = {
        u[idx_x[0]], u[idx_x[1]], u[idx_x[2]], u[idx_x[3]]
      };
      std::array<T, 4> uy = {
        u[idx_y[0]], u[idx_y[1]], u[idx_y[2]], u[idx_y[3]]
      };
      std::array<T, 4> uz = {
        u[idx_z[0]], u[idx_z[1]], u[idx_z[2]], u[idx_z[3]]
      };
      std::array<real_t, 2> fx = { load[idx_load[0]][0], load[idx_load[1]][0] };
      std::array<real_t, 2> fy = { load[idx_load[0]][1], load[idx_load[1]][1] };
      std::array<real_t, 2> fz = { load[idx_load[0]][2], load[idx_load[1]][2] };
      std::array<T, 2> ul = { lambda[idx_l[0]], lambda[idx_l[1]] };

      std::vector<T> R_loc_x(4, 0);
      std::vector<T> R_loc_y(4, 0);
      std::vector<T> R_loc_z(4, 0);

      for (size_t qi = 0; qi < 3; ++qi) {
        real_t xi = xi_q[qi];
        real_t w = w_q[qi];

        auto H = ELFF::FEM::CubicHermite<real_t>::values(xi, ds);
        auto dH = ELFF::FEM::CubicHermite<real_t>::derivs(xi, ds);
        auto ddH = ELFF::FEM::CubicHermite<real_t>::second_derivs(xi, ds);
        auto M = ELFF::FEM::LinearShape<real_t>::values(xi);

        T x = 0, xp = 0, xpp = 0;
        T y = 0, yp = 0, ypp = 0;
        T z = 0, zp = 0, zpp = 0;
        for (size_t i = 0; i < 4; i++) {
          x += H[i] * ux[i];
          xp += dH[i] * ux[i];
          xpp += ddH[i] * ux[i];
          y += H[i] * uy[i];
          yp += dH[i] * uy[i];
          ypp += ddH[i] * uy[i];
          z += H[i] * uz[i];
          zp += dH[i] * uz[i];
          zpp += ddH[i] * uz[i];
        }

        T l = 0;
        real_t fxp = 0;
        real_t fyp = 0;
        real_t fzp = 0;
        for (size_t i = 0; i < 2; i++) {
          l += M[i] * ul[i];
          fxp += M[i] * fx[i];
          fyp += M[i] * fy[i];
          fzp += M[i] * fz[i];
        }

        T S = xp * xp + yp * yp + zp * zp - 1.0;

        for (size_t a = 0; a < 4; ++a) {
          R_loc_x[a] += EI * xpp * ddH[a] * w * ds;
          R_loc_y[a] += EI * ypp * ddH[a] * w * ds;
          R_loc_z[a] += EI * zpp * ddH[a] * w * ds;
          R_loc_x[a] -= fxp * H[a] * w * ds;
          R_loc_y[a] -= fyp * H[a] * w * ds;
          R_loc_z[a] -= fzp * H[a] * w * ds;
          R_loc_x[a] += 2 * (l + r_penalty * S) * xp * dH[a] * w * ds;
          R_loc_y[a] += 2 * (l + r_penalty * S) * yp * dH[a] * w * ds;
          R_loc_z[a] += 2 * (l + r_penalty * S) * zp * dH[a] * w * ds;
        }
      }

      for (size_t i = 0; i < 4; ++i) {
        residual[idx_x[i]] += R_loc_x[i];
        residual[idx_y[i]] += R_loc_y[i];
        residual[idx_z[i]] += R_loc_z[i];
      }
    }
    return residual;
  }
};

} // namespace ELFF
} // namespace Models
