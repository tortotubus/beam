#include "elff/models/beam/EulerBeamMesh.hpp"
#include "elff/general/error.hpp"

namespace ELFF {
namespace Models {

EulerBeamMesh::EulerBeamMesh(size_t nodes, real_t length)
  : nodes(nodes)
  , length(length)
  , ds(length / (nodes - 1))
  , s(nodes)
  , centerline(nodes)
  , slope(nodes)
  , centerline_velocity(nodes)
  , centerline_acceleration(nodes)
  , slope_velocity(nodes)
  , slope_acceleration(nodes)
{
  set_curvilinear_axis();
  set_centerline_x_axis_aligned();
}

void
EulerBeamMesh::plot_gnuplot(std::string title)
{
  FILE* pipe = popen("gnuplot -persist", "w");
  if (!pipe) {
    ELFF_WARNING("Failed to open pipe to gnuplot");
    return;
  }

  // Configure the plot
  fprintf(pipe, "set title '%s'\n", title.c_str());
  fprintf(pipe, "set xlabel 'x'\n");
  fprintf(pipe, "set ylabel 'y'\n");
  fprintf(pipe, "set grid\n");
  fprintf(pipe, "set size square\n");
  fprintf(pipe,
          "plot '-' using 1:2 with lines title 'beam'\n");

  // Send the data points
  for (size_t i = 0; i < centerline.size(); ++i) {
    fprintf(pipe, "%f %f\n", centerline[i][0], centerline[i][1]);
  }
  fprintf(pipe, "e\n"); // End of data marker for gnuplot

  // Clean up
  pclose(pipe);
}

void
EulerBeamMesh::plot_gnuplot()
{
  plot_gnuplot("Beam Centerline");
}

IO::CXX::vtkPolyData
EulerBeamMesh::to_vtk_polydata()
{
  IO::CXX::vtkPolyData pd;

  pd.reserve_points(nodes);
  pd.reserve_lines(nodes - 1);

  for (size_t i = 0; i < nodes; i++)
    pd.add_point(centerline[i][0], centerline[i][1], centerline[i][2]);

  for (size_t i = 0; i < nodes - 1; i++)
    pd.add_line(i, i + 1);

  const int64_t centerline_velocity_field =
    pd.add_pointdata_vector("centerline_velocity", 3);

  for (size_t i = 0; i < nodes; ++i) {
    pd.set_pointdata_vector3(centerline_velocity_field,
                             i,
                             { static_cast<double>(centerline_velocity[i][0]),
                               static_cast<double>(centerline_velocity[i][1]),
                               static_cast<double>(centerline_velocity[i][2]) });
  }

  return pd;
}

void
EulerBeamMesh::set_centerline_velocity_zero()
{
  for (size_t i = 0; i < nodes; i++) {
    centerline_velocity[i] = { 0., 0., 0. };
  }
}

void
EulerBeamMesh::set_curvilinear_axis()
{
  for (size_t i = 0; i < nodes; ++i) {
    s[i] = ds * i;
  }
}

void
EulerBeamMesh::set_centerline_x_axis_aligned()
{
  for (size_t i = 0; i < nodes; ++i) {
    centerline[i] = { ds * i, 0., 0. };
    slope[i] = { 1., 0., 0. };
  }
}

} // namespace Models
} // namespace ELFF
