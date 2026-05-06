#pragma once

#include "elff/config/config.hpp"
#include "elff/models/beam/EulerBeamMesh.hpp"

#include <array>
#include <vector>

namespace ELFF {
namespace Models {

/**
 * @brief Base class for Euler-Bernoulli beam implementations.
 *
 * This abstract class provides the interface for solving Euler-Bernoulli beam
 * problems with different boundary conditions, loading scenarios, and solution
 * methods. The class supports both static and dynamic problems, as well as
 * uniform and non-uniform loading conditions.
 *
 * Each beam is characterized by:
 * - Flexural rigidity \f(EI\f)
 * - Mass per unit length \f(\mu\f) for dynamic problems
 * - Boundary conditions at both ends
 * - Spatial discretization (mesh)
 * - Dimension (2D or 3D)
 */
class EulerBeam {
public:
  /**
   * @brief Boundary-condition types for Euler–Bernoulli / elastica beams.
   *
   * Use these to specify end constraints or lumped end loads.
   *  - free_bc:         shear \f(V(s) = 0\f) and bending moment \f(M(s) = 0\f)
   *  - simple_bc:       deflection \f(r(s) = 0\f) and moment \f(M(s) = 0\f)
   *  - clamped_bc:      deflection \f(r(s) = 0\f) and rotation \f(θ(s) = 0\f)
   *  - point_force_bc:  jump in shear (\f(V\f) has an impulse)
   *  - point_torque_bc: jump in moment (\f(M\f) has an impulse)
   *
   * @see EulerBeamBC for how these are carried with magnitudes/locations.
   * @memberof EulerBeamBCs
   */
  typedef enum {
    free_bc,
    simple_bc,
    clamped_bc,
    point_force_bc,
    point_torque_bc,
  } EulerBeamBCType;

  /**
   * @brief Beam end to denote where to apply a given boundary condition.
   *
   * Note that "left" should be applied to the \f(s = 0\f), \f(x = 0\f), or
   * otherwise whatever is first coordinate in the particular E-B beam model's
   * netural axis coordinate system. Correspondingly the "right" should be
   * applied to the final coordinate \f(s = L\f) or \f(x = Lf\).
   *
   * @memberof EulerBeamBCs
   */
  typedef enum {
    left,
    right,
  } EulerBeamBCEnd;

  /**
   * @brief Values to apply for a given boundary condition.
   *
   * Note that each is a vector in \f(\mathbb{R}^3\f), extraneous coordinats in
   * particular models will be ignored if they are not used. Similarly,
   * extraneous values for a given boundary condition to be applied will also be
   * ignored.
   *
   * @memberof EulerBeamBCs
   */
  typedef struct {
    double position[3];
    double slope[3];
    double force[3];
    double torque[3];
  } EulerBeamBCVals;

  /**
   * @brief The boundary condition struct contains user-supplied boundary
   * conditions to be applied to the problem solved.
   *
   *
   * @see EulerBeamBCType for the type of boundary conditions supported.
   * @memberof EulerBeam
   */
  typedef struct {
    EulerBeamBCEnd end[2];
    EulerBeamBCType type[2];
    EulerBeamBCVals vals[2];
  } EulerBeamBCs;

  /**
   * @brief Boundary values for time-dependent boundary conditions including
   * imposed values and kinematic data
   */
  typedef struct {
    std::vector<std::array<real_t, 3>> position_history;
    std::vector<std::array<real_t, 3>> slope_history;
    std::vector<std::array<real_t, 3>> force_history;
    std::vector<std::array<real_t, 3>> torque_history;
    std::vector<std::array<real_t, 3>> velocity_history;
    std::vector<std::array<real_t, 3>> acceleration_history;
    std::vector<std::array<real_t, 3>> slope_velocity_history;
    std::vector<std::array<real_t, 3>> slope_acceleration_history;
  } EulerBeamBoundaryHistory;

  /**
   * @brief The boundary condition struct contains user-supplied time-dependent
   * boundary conditions to be applied to the problem solved.
   *
   * The normal structs are used for current time boundary conditions, and
   * historical values are kept in the history member. Derived classes will
   * impose the number of historical time-steps required at each time step,
   * compare this to time_zero_idx, and throw errors if this is not satisfied.
   * If solve() is called and runs out of prescribed values, errors will also be
   * thrown.
   *
   * @see EulerBeamBCType for the type of boundary conditions supported.
   * @memberof EulerBeam
   */
  typedef struct {
    EulerBeamBCEnd end[2];
    EulerBeamBCType type[2];
    EulerBeamBoundaryHistory history[2];
    size_t time_zero_idx;
  } EulerBeamTimeDependentBCs;

  /**
   * @brief Solves the beam problem with no external loading.
   *
   * This virtual method should be implemented by derived classes to solve
   * the beam problem under the specified boundary conditions without any
   * external loads.
   */
  virtual void solve();

  /**
   * @brief Solves the beam problem under uniform loading.
   *
   * @param uniform_load A 3D vector representing the uniform load per unit
  length applied to the beam. Components represent forces in \f(x\f), \f(y\f),
  and \f(z\f) directions.
   */
  virtual void solve(std::array<real_t, 3> uniform_load);

  /**
   * @brief Solves the beam problem under non-uniform loading.
   *
   * @param nonuniform_load Vector of 3D force vectors representing the load
   *        at each node of the beam mesh.
   */
  virtual void solve(std::vector<std::array<real_t, 3>> nonuniform_load);

  /**
   * @brief Advances the dynamic beam solution by one time step.
   *
   * @param dt Time step size for temporal discretization
   */
  virtual void solve(real_t dt);

  /**
   * @brief Advances the dynamic beam solution under uniform loading by one time
   * step.
   *
   * @param dt Time step size for temporal discretization
   * @param uniform_load A 3D vector representing the uniform load per unit
   * length
   */
  virtual void solve(real_t dt, std::array<real_t, 3> uniform_load);

  /**
   * @brief Advances the dynamic beam solution under non-uniform loading by one
   * time step.
   *
   * @param dt Time step size for temporal discretization
   * @param nonuniform_load Vector of 3D force vectors for each node
   */
  virtual void solve(real_t dt,
                     std::vector<std::array<real_t, 3>> nonuniform_load);

  /**
   * @brief Updates the beam boundary conditions.
   *
   * Derived classes may override this to refresh any cached operators or time
   * integration history tied to constrained degrees of freedom.
   *
   * @param bcs New boundary conditions to apply
   */
  virtual void
  set_time_dependent_boundary_conditions(EulerBeamTimeDependentBCs bcs);

  /**
   * @brief Applies initial conditions to the beam mesh.
   *
   * Derived classes must implement this method to set the initial
   * displacement and velocity fields for dynamic problems.
   *
   * @param mesh The beam mesh to initialize
   */
  virtual void apply_initial_condition(EulerBeamMesh &mesh) = 0;

  /**
   * @brief Gets a mutable reference to the beam mesh.
   *
   * @return Reference to the EulerBeamMesh object containing the beam geometry
   */
  virtual EulerBeamMesh &get_mesh();

  /**
   * @brief Gets a const reference to the beam mesh.
   *
   * @return Const reference to the EulerBeamMesh object containing the beam
   * geometry
   */
  virtual const EulerBeamMesh &get_mesh() const;

protected:
  bool is_time_dependent;
  size_t time_iter;
  real_t t;

  size_t dimension;

  real_t EI;
  real_t mu;

  EulerBeamMesh mesh;
  EulerBeamBCs boundary_conditions;
  EulerBeamTimeDependentBCs time_dependent_boundary_conditions;
  bool have_time_dependent_boundary_conditions;

  /**
   * @brief Default constructor for static 3D beam.
   *
   * Creates a beam with:
   * - Unit flexural rigidity (\f(EI = 1\f))
   * - Length \f(L=1\f)
   * - 20 nodes
   * - Left end clamped, right end free
   * - No mass (static problem)
   */
  EulerBeam();

  /**
   * @brief Constructor for static 3D beam with custom parameters.
   *
   * @param length The length of the beam
   * @param EI The flexural rigidity
   * @param nodes Number of nodes in the discretization
   * @param bcs Boundary conditions at both ends
   */
  EulerBeam(real_t length, real_t EI, size_t nodes, EulerBeamBCs bcs);

  /**
   * @brief Constructor for dynamic 3D beam with custom parameters.
   *
   * @param length The length of the beam
   * @param EI The flexural rigidity
   * @param mu Mass per unit length
   * @param nodes Number of nodes in the discretization
   * @param bcs Boundary conditions at both ends
   */
  EulerBeam(real_t length, real_t EI, real_t mu, size_t nodes,
            EulerBeamBCs bcs);

/**
 * @brief Constructor for dynamic 3D beam with custom parameters and
 * time-dependent boundary conditions
 *
 * @param length The length of the beam
 * @param EI The flexural rigidity
 * @param mu Mass per unit length
 * @param nodes Number of nodes in the discretization
 * @param bcs Boundary conditions at both ends
 */
 EulerBeam(real_t length, real_t EI, real_t mu, size_t nodes,
                     EulerBeamTimeDependentBCs bcs);

};

} // namespace Models
} // namespace ELFF
