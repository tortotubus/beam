#pragma once

#include "config/config.hpp"
#include "general/error.hpp"
#include "io/VTKPolyData.hpp"

#include <cstdio>
#include <cstdlib>

#include <array>
// #include <cmath>
#include <vector>

namespace beam {
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

public:
  EulerBeamMesh(size_t nodes, real_t length)
    : nodes(nodes)
    , length(length)
    , ds(length / (nodes - 1))
    , s(nodes)
    , centerline(nodes)
    , slope(nodes)
    , centerline_velocity(nodes)
  {
    set_curvilinear_axis();
    set_centerline_x_axis_aligned();
  }

  /**
   * @brief Returns a reference to the entire centerline \f(\{\vec{r}_0,
   * \dots, \vec{r}_n\}\f).
   */
  inline std::vector<std::array<real_t, 3>>& get_centerline()
  {
    return centerline;
  }

  /**
   * @brief Returns one vector \f(\vec{r}_i\f) contained in the centerline.
   */
  inline std::array<real_t, 3>& get_centerline(size_t i)
  {
    return centerline[i];
  }

  /**
   * @brief Returns a reference to all tangents of the centerline
   * \f(\{\vec{r}'_0, \dots, \vec{r}'_n\}\f).
   */
  inline std::vector<std::array<real_t, 3>>& get_slope() { return slope; }

  /**
   * @brief Returns a reference to one tangent of the centerline
   * \f(\vec{r}'_i\f).
   */
  inline std::array<real_t, 3>& get_slope(size_t i) { return slope[i]; }

  /**
   * @brief Returns a reference to the entire centerline velocity
   * \f(\{\dot\vec{r}_0,
   * \dots, \dot\vec{r}_n\}\f).
   */
  inline std::vector<std::array<real_t, 3>>& get_centerline_velocity()
  {
    return centerline_velocity;
  }

  /**
   * @brief Returns one vector \f(\dot \vec{r}_i\f) contained in the centerline
   * velocity.
   */
  inline std::array<real_t, 3>& get_centerline_velocity(size_t i)
  {
    return centerline_velocity[i];
  }

  /**
   * @brief Returns a reference to the curvilinear axis \f(s\in[0,L]\f).
   */
  inline std::vector<real_t>& get_curvilinear_axis() { return s; }

  /**
   * @brief Returns a reference to one point \f(s_i \in [0,L]\f).
   */
  inline real_t get_curvilinear_axis(size_t i) { return s[i]; }

  /**
   * @brief Returns \f(\rm{d}s\f).
   */
  inline real_t get_ds() const { return ds; }

  /**
   * @brief Returns the number of nodes.
   */
  inline size_t get_nodes() const { return nodes; }

  /**
   * @brief Returns the length of the curvilinear axis \f(s\f).
   */
  inline real_t get_length() const { return length; }

  /**
   * @brief Plots the mesh using gnuplot, if installed.
   */
  void plot_gnuplot(std::string title)
  {
    FILE* pipe = popen("gnuplot -persist", "w");
    if (!pipe) {
      BEAM_ABORT("Failed to open pipe to gnuplot");
    }

    // Configure the plot
    fprintf(pipe, "set title '%s'\n", title.c_str());
    fprintf(pipe, "set xlabel 'x'\n");
    fprintf(pipe, "set ylabel 'y'\n");
    // fprintf(pipe, "set xrange [%f:%f]\n", -curvlinear_length,
    // curvlinear_length); fprintf(pipe, "set yrange [%f:%f]\n",
    // -curvlinear_length, curvlinear_length);
    fprintf(pipe, "set grid\n");
    fprintf(pipe, "set size square\n"); // <= disable extension to tics
                                        // :contentReference[oaicite:0]{index=0}
    fprintf(pipe,
            "plot '-' using 1:2 with lines title 'beam' smooth csplines\n");

    // Send the data points
    for (size_t i = 0; i < centerline.size(); ++i) {
      fprintf(pipe, "%f %f\n", centerline[i][0], centerline[i][1]);
    }
    fprintf(pipe, "e\n"); // End of data marker for gnuplot

    // Clean up
    pclose(pipe);
  }

  void plot_gnuplot() { plot_gnuplot("Beam Centerline"); }

  VTKPolyData to_vtk_polydata()
  {
    VTKPolyData pd;

    pd.reserve_points(nodes);
    pd.reserve_lines(nodes - 1);

    for (size_t i = 0; i < nodes; i++)
      pd.add_point(centerline[i][0], centerline[i][1], centerline[i][2]);

    for (size_t i = 0; i < nodes - 1; i++)
      pd.add_line(i, i + 1);

    return pd;
  }

private:
  void set_centerline_velocity_zero()
  {
    for (size_t i = 0; i < nodes; i++) {
      centerline_velocity[i] = { 0., 0., 0. };
    }
  }

  void set_curvilinear_axis()
  {
    for (size_t i = 0; i < nodes; ++i) {
      s[i] = ds * i;
    }
  }

  void set_centerline_x_axis_aligned()
  {
    for (size_t i = 0; i < nodes; ++i) {
      centerline[i] = { ds * i, 0., 0. };
      slope[i] = { 1., 0., 0. };
    }
  }
};

}