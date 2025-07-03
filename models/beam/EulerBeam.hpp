#pragma once

#include <cstdlib>
#include <cstdio>

#include <array>
#include <vector>

#include "config/config.hpp"
#include "general/error.hpp"

namespace beam {

typedef enum
{
  free_bc,
  simple_bc,
  clamped_bc,
} EulerBeamBCType;

typedef enum
{
  left,
  right,
} EulerBeamBCEnd;

typedef struct
{
  double position[3];
  double slope[3];
} EulerBeamBCVals;

typedef struct
{
  EulerBeamBCEnd end[2];
  EulerBeamBCType type[2];
  EulerBeamBCVals vals[2];
} EulerBeamBCs;

class EulerBeamMesh
{
protected:
  size_t n_nodes;
  real_t curvlinear_length, ds;
  std::vector<real_t> s;
  std::vector<std::array<real_t, 3>> centerline;

public:
  EulerBeamMesh(size_t n_nodes, real_t length)
    : n_nodes(n_nodes)
    , curvlinear_length(length)
    , ds(curvlinear_length / (n_nodes-1))
    , s(n_nodes)
    , centerline(n_nodes)
  {
    set_curvilinear_axis();
    set_centerline_x_axis_aligned();
  }

  inline std::vector<std::array<real_t, 3>>& get_centerline() { return centerline; }
  inline std::array<real_t, 3>& get_centerline(size_t i) { return centerline[i]; }

  inline std::vector<real_t>& get_curvilinear_axis() { return s; }
  inline real_t get_curvilinear_axis(size_t i) { return s[i]; }

  inline real_t get_ds() { return ds; }
  inline size_t get_n_nodes() { return n_nodes; }

  void plot_gnuplot() {
    FILE *pipe = popen("gnuplot -persist", "w");
    if (!pipe) {
      BEAM_ABORT("Failed to open pipe to gnuplot");
    }


    // Configure the plot
    fprintf(pipe, "set title 'Beam Centerline'\n");
    fprintf(pipe, "set xlabel 'x'\n");
    fprintf(pipe, "set ylabel 'y'\n");
    fprintf(pipe, "set grid\n");
    fprintf(pipe, "set autoscale fix\n");      // <= disable extension to tics :contentReference[oaicite:0]{index=0}  
    fprintf(pipe, "plot '-' using 1:2 with lines title 'beam' smooth csplines\n");

    // Send the data points
    for (size_t i = 0; i < centerline.size(); ++i) {
      fprintf(pipe, "%f %f\n", centerline[i][0], centerline[i][1]);
    }
    fprintf(pipe, "e\n"); // End of data marker for gnuplot

    // Clean up
    pclose(pipe);
  }

private:
  void set_curvilinear_axis()
  {
    for (size_t i = 0; i < n_nodes; ++i) {
      s[i] = ds * i;
    }
  }

  void set_centerline_x_axis_aligned()
  {
    for (size_t i = 0; i < n_nodes; ++i) {
      centerline[i] = { ds * i, 0., 0. };
    }
  }
};

typedef enum { INVALID = -1 , NONE, UNIFORM, NONUNIFORM } EulerBeamLoadType;

class EulerBeam
{
public:
  EulerBeam() : EI(1), area(1), dimension(2), mesh(20, 1.), load_type(NONE), uniform_load({0,0,0}), nonuniform_load(20), 
    boundary_conditions(
      { .end = { left, right },
        .type = { clamped_bc, free_bc },
        .vals = { { .position = { 0., 0., 0. }, .slope = { 0., 0., 0. } },
                  { .position = { 0., 0., 0. }, .slope = { 0., 0., 0. } } } }) {
    };

  EulerBeam(real_t length,
            real_t EI,
            real_t area,
            size_t nodes,
            EulerBeamBCs bcs) :
    EI(EI), area(area), dimension(2), mesh(nodes, length), 
    boundary_conditions(bcs) {};

  virtual void solve() = 0;

  virtual void plot() { mesh.plot_gnuplot(); }

  virtual void set_load(std::array<real_t,3> load) {
    load_type = UNIFORM;
    uniform_load = load;
  }

  virtual void set_load(std::vector<std::array<real_t,3>> load) {
    BEAM_ASSERT(load.size() == mesh.get_n_nodes(), "Nonuniform load vector must have the same size as the mesh nodes.");
    load_type = NONUNIFORM;
    nonuniform_load = load;
  }

protected:
  real_t EI;
  real_t area;
  size_t dimension;
  EulerBeamMesh mesh;
  EulerBeamBCs boundary_conditions;

  EulerBeamLoadType load_type;
  std::array<real_t,3> uniform_load;
  std::vector<std::array<real_t,3>> nonuniform_load;

};

} // namespace beam
