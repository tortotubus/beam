#pragma once

#include "EulerBeamMesh.hpp"
#include "config/config.hpp"

#include <array>
#include <vector>

namespace beam {

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
class EulerBeam
{
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
  typedef enum
  {
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
  typedef enum
  {
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
  typedef struct
  {
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
  typedef struct
  {
    EulerBeamBCEnd end[2];
    EulerBeamBCType type[2];
    EulerBeamBCVals vals[2];
  } EulerBeamBCs;

  /**
   * @brief Denotes the type of loading applied to the beam at nodes. In
   * particular, we support uniform (spatially invariant) loads represented as
   * constants as well as nonuniform (spatially depdendent) loads at the nodes.
   *
   * For more details, see the solve() methods in EulerBeam base class.
   */
  typedef enum
  {
    INVALID = -1,
    NONE,
    UNIFORM,
    NONUNIFORM,
    // UNIFORM_TIME_INDEPENDENT,
    // NONUNIFORM_TIME_INDEPENDENT,
    // UNIFORM_TIME_DEPENDENT,
    // NONUNIFORM_TIME_DEPENDENT,
  } EulerBeamLoadType;

  /**
   * @brief Solves the beam problem with no external loading.
   *
   * This virtual method should be implemented by derived classes to solve
   * the beam problem under the specified boundary conditions without any
   * external loads.
   */
  virtual void solve()
  {
    BEAM_ABORT("EulerBeam base class does not implement a solution method.\n");
  };

  /**
   * @brief Solves the beam problem under uniform loading.
   *
   * @param uniform_load A 3D vector representing the uniform load per unit
  length applied to the beam. Components represent forces in \f(x\f), \f(y\f),
  and \f(z\f) directions.
   */
  virtual void solve(std::array<real_t, 3> uniform_load)
  {
    BEAM_ABORT("EulerBeam base class does not implement a solution method.\n");
  };

  /**
   * @brief Solves the beam problem under non-uniform loading.
   *
   * @param nonuniform_load Vector of 3D force vectors representing the load
   *        at each node of the beam mesh.
   */
  virtual void solve(std::vector<std::array<real_t, 3>> nonuniform_load)
  {
    BEAM_ABORT("EulerBeam base class does not implement a solution method.\n");
  };

  /**
   * @brief Advances the dynamic beam solution by one time step.
   *
   * @param dt Time step size for temporal discretization
   */
  virtual void solve(real_t dt)
  {
    BEAM_ABORT("EulerBeam base class does not implement a solution method.\n");
  }

  /**
   * @brief Advances the dynamic beam solution under uniform loading by one time
   * step.
   *
   * @param dt Time step size for temporal discretization
   * @param uniform_load A 3D vector representing the uniform load per unit
   * length
   */
  virtual void solve(real_t dt, std::array<real_t, 3> uniform_load)
  {
    BEAM_ABORT("EulerBeam base class does not implement a solution method.\n");
  }

  /**
   * @brief Advances the dynamic beam solution under non-uniform loading by one
   * time step.
   *
   * @param dt Time step size for temporal discretization
   * @param nonuniform_load Vector of 3D force vectors for each node
   */
  virtual void solve(real_t dt,
                     std::vector<std::array<real_t, 3>> nonuniform_load)
  {
    BEAM_ABORT("EulerBeam base class does not implement a solution method.\n");
  }

  /**
   * @brief Applies initial conditions to the beam mesh.
   *
   * Derived classes must implement this method to set the initial
   * displacement and velocity fields for dynamic problems.
   *
   * @param mesh The beam mesh to initialize
   */
  virtual void apply_initial_condition(EulerBeamMesh& mesh) = 0;

  /**
   * @brief Plots the current beam configuration using gnuplot.
   *
   * Creates a visualization of the beam's centerline with default title.
   */
  virtual void plot() { mesh.plot_gnuplot(); }

  /**
   * @brief Plots the current beam configuration using gnuplot with custom
   * title.
   *
   * @param title The title to display on the plot
   */
  virtual void plot(std::string title) { mesh.plot_gnuplot(title); }

  /**
   * @brief Gets a mutable reference to the beam mesh.
   *
   * @return Reference to the EulerBeamMesh object containing the beam geometry
   */
  virtual EulerBeamMesh& get_mesh() { return mesh; }

  /**
   * @brief Gets a const reference to the beam mesh.
   *
   * @return Const reference to the EulerBeamMesh object containing the beam
   * geometry
   */
  virtual const EulerBeamMesh& get_mesh() const { return mesh; }

protected:
  bool is_time_dependent;
  size_t time_iter;
  real_t t;

  size_t dimension;

  real_t EI;
  real_t mu;

  EulerBeamMesh mesh;
  EulerBeamBCs boundary_conditions;

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
  EulerBeam()
    : EI(1)
    , mu(0)
    , is_time_dependent(false)
    , time_iter(0)
    , t(0)
    , dimension(3)
    , mesh(20, 1.)
    , boundary_conditions(
        { .end = { left, right },
          .type = { clamped_bc, free_bc },
          .vals = {
            { .position = { 0., 0., 0. }, .slope = { 1., 0., 0. } },
            { .position = { 0., 0., 0. }, .slope = { 0., 0., 0. } } } }) {};

  /**
   * @brief Constructor for static 3D beam with custom parameters.
   *
   * @param length The length of the beam
   * @param EI The flexural rigidity
   * @param nodes Number of nodes in the discretization
   * @param bcs Boundary conditions at both ends
   */
  EulerBeam(real_t length, real_t EI, size_t nodes, EulerBeamBCs bcs)
    : EI(EI)
    , mu(0)
    , is_time_dependent(false)
    , time_iter(0)
    , t(0)
    , dimension(3)
    , mesh(nodes, length)
    , boundary_conditions(bcs) {};

  /**
   * @brief Constructor for dynamic 3D beam with custom parameters.
   *
   * @param length The length of the beam
   * @param EI The flexural rigidity
   * @param mu Mass per unit length
   * @param nodes Number of nodes in the discretization
   * @param bcs Boundary conditions at both ends
   */
  EulerBeam(real_t length, real_t EI, real_t mu, size_t nodes, EulerBeamBCs bcs)
    : EI(EI)
    , mu(mu)
    , is_time_dependent(true)
    , time_iter(0)
    , t(0)
    , dimension(3)
    , mesh(nodes, length)
    , boundary_conditions(bcs) {};
};

} // namespace beam
