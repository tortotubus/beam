#pragma once

#include "elff/config/config.hpp"
#include "elff/general/error.hpp"
#include "elff/io/CXX/vtkPolyData.hpp"

#include <cstdio>
#include <cstdlib>

#include <array>
#include <vector>

namespace ELFF {
namespace Models {
/**
 * @brief This class helps provide a standard mesh for the various methods to
 * use. All implementations of E-B methods should update the mesh at the end
 * of any solve() routine from internal representations of the solution
 * to this class.
 *
 * This mesh is equipped with a curvilinear axis s which maps \f(s\in[0,L]\f)
 * to the centerline position \f(\vec{r}(s) = (x(s),y(s),z(s))\f) as well as
 * centerline tangent \f( \vec{r}'(s) = (x'(s),y'(s),z'(s))\f).
 */

class EulerBeamMesh
{
protected:
  size_t nodes;
  real_t length, ds;
  std::vector<real_t> s;

  std::vector<std::array<real_t, 3>> centerline, slope;
  std::vector<std::array<real_t, 3>> centerline_velocity;
  std::vector<std::array<real_t, 3>> centerline_acceleration;
  std::vector<std::array<real_t, 3>> slope_velocity;
  std::vector<std::array<real_t, 3>> slope_acceleration;

public:
  /**
   * @brief Construct an EulerBeamMesh with evenly spaced nodes along the beam.
   *
   * The mesh will contain `nodes` discretization points uniformly distributed
   * along the curvilinear axis \f(s\in[0,L]\f). The spacing is computed as
   * \f(\Delta s = L/(N-1)\f) and stored in the member `ds`. After allocation
   * the constructor populates the curvilinear axis and aligns the centerline
   * along the x-axis by calling the internal helpers `set_curvilinear_axis()`
   * and `set_centerline_x_axis_aligned()`.
   *
   * @param nodes Number of discretization nodes N (must be >= 2).
   * @param length Physical length L of the centerline.
   */
  EulerBeamMesh(size_t nodes, real_t length);

  EulerBeamMesh() : EulerBeamMesh(40, 1) {};

  /**
   * @brief Returns a reference to the entire centerline \f(\{\vec{r}_0,
   * \dots, \vec{r}_n\}\f).
   *
   * @return Reference to the vector of 3-component positions
   *         \f(\{\vec{r}_0, \dots, \vec{r}_n\}\f)
   */
  inline std::vector<std::array<real_t, 3>>& get_centerline()
  {
    return centerline;
  }

  /**
   * @brief Returns one vector \f(\vec{r}_i\f) contained in the centerline.
   *
   * @param i Index of the requested centerline node (0..nodes-1)
   * @return Reference to the 3-component position \f(\vec{r}_i\f)
   */
  inline std::array<real_t, 3>& get_centerline(size_t i)
  {
    return centerline[i];
  }

  /**
   * @brief Returns a reference to all tangents of the centerline
   * \f(\{\vec{r}'_0, \dots, \vec{r}'_n\}\f).
   *
   * @return Reference to the vector of 3-component tangent vectors
   */
  inline std::vector<std::array<real_t, 3>>& get_slope() { return slope; }

  /**
   * @brief Returns a reference to one tangent of the centerline
   * \f(\vec{r}'_i\f).
   *
   * @param i Index of the requested tangent (0..nodes-1)
   * @return Reference to the 3-component tangent \f(\vec{r}'_i\f)
   */
  inline std::array<real_t, 3>& get_slope(size_t i) { return slope[i]; }

  /**
   * @brief Returns a reference to the entire centerline velocity
   * \f(\{\dot{\vec{r}}_0,
   * \dots,\dot{\vec{r}}_n\}\f).
   *
   * @return Reference to the vector of centerline velocities
   */
  inline std::vector<std::array<real_t, 3>>& get_centerline_velocity()
  {
    return centerline_velocity;
  }

  /**
   * @brief Returns a reference to the entire centerline velocity
   * \f(\{\ddot{\vec{r}}_0,
   * \dots,\ddot{\vec{r}}_n\}\f).
   *
   * @return Reference to the vector of centerline accelerations
   */
  inline std::vector<std::array<real_t, 3>>& get_centerline_acceleration()
  {
    return centerline_acceleration;
  }

  /**
   * @brief Returns a reference to the entire slope velocity
   * \f(\{\partial_t \vec{r}'_0, \dots, \partial_t \vec{r}'_n\}\f).
   *
   * @return Reference to the vector of slope velocities
   */
  inline std::vector<std::array<real_t, 3>>& get_slope_velocity()
  {
    return slope_velocity;
  }

  /**
   * @brief Returns a reference to the entire slope acceleration
   * \f(\{\partial_{tt} \vec{r}'_0, \dots, \partial_{tt} \vec{r}'_n\}\f).
   *
   * @return Reference to the vector of slope accelerations
   */
  inline std::vector<std::array<real_t, 3>>& get_slope_acceleration()
  {
    return slope_acceleration;
  }

  /**
   * @brief Returns one vector \f(\dot{\vec{r}}_i\f) contained in the centerline
   * velocity.
   *
   * @param i Index of the requested velocity (0..nodes-1)
   * @return Reference to the 3-component velocity \f(\dot{\vec{r}}_i\f)
   */
  inline std::array<real_t, 3>& get_centerline_velocity(size_t i)
  {
    return centerline_velocity[i];
  }

  /**
   * @brief Returns one vector \f(\partial_t \vec{r}'_i\f) contained in the
   * slope velocity.
   *
   * @param i Index of the requested slope velocity (0..nodes-1)
   * @return Reference to the 3-component slope velocity
   */
  inline std::array<real_t, 3>& get_slope_velocity(size_t i)
  {
    return slope_velocity[i];
  }

  /**
   * @brief Returns one vector \f(\partial_{tt} \vec{r}'_i\f) contained in the
   * slope acceleration.
   *
   * @param i Index of the requested slope acceleration (0..nodes-1)
   * @return Reference to the 3-component slope acceleration
   */
  inline std::array<real_t, 3>& get_slope_acceleration(size_t i)
  {
    return slope_acceleration[i];
  }

  /**
   * @brief Returns a reference to the curvilinear axis \f(s\in[0,L]\f).
   *
   * @return Reference to the vector of curvilinear coordinates \f(s_i\f)
   */
  inline std::vector<real_t>& get_curvilinear_axis() { return s; }

  /**
   * @brief Returns a reference to one point \f(s_i \in [0,L]\f).
   *
   * @param i Index of the curvilinear coordinate (0..nodes-1)
   * @return The curvilinear coordinate value \f(s_i\f)
   */
  inline real_t get_curvilinear_axis(size_t i) { return s[i]; }

  /**
   * @brief Returns \f(\rm{d}s\f), the uniform spacing of the curvilinear axis.
   *
   * @return Spacing \f(\Delta s = ds = L/(N-1)\f)
   */
  inline real_t get_ds() const { return ds; }

  /**
   * @brief Returns the number of nodes.
   *
   * @return Number of discretization nodes along the centerline
   */
  inline size_t get_nodes() const { return nodes; }

  /**
   * @brief Returns the length of the curvilinear axis \f(s\f).
   *
   * @return Physical length \f(L\f) of the beam centerline
   */
  inline real_t get_length() const { return length; }

  /**
   * @brief Plots the mesh using gnuplot, if installed.
   */
  void plot_gnuplot(std::string title);

  /**
   * @brief Plot the mesh using gnuplot with a default title.
   *
   * This is a convenience wrapper equivalent to calling
   * plot_gnuplot("Beam Centerline"). If gnuplot is not available the
   * function will abort with an error.
   */
  void plot_gnuplot();

  /**
   * @brief Export the centerline geometry as an in-memory VTK PolyData object.
   *
   * The returned VTK wrapper will contain the centerline points and a line
   * connectivity linking successive nodes. The function pre-allocates point
   * and line storage for efficiency.
   *
   * @return A value (move) of IO::CXX::vtkPolyData containing the centerline
   *         points and line cells. The caller takes ownership of the returned
   *         object and its internal memory.
   */
  IO::CXX::vtkPolyData to_vtk_polydata();

private:
  void set_centerline_velocity_zero();

  void set_curvilinear_axis();

  void set_centerline_x_axis_aligned();
};

} // namespace Models
} // namespace ELFF 
