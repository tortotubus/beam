#include "elff/models/rope/RopeMesh.hpp"

namespace ELFF {
namespace Models {

RopeMesh::RopeMesh(size_t nodes, real_t length)
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
 
IO::CXX::vtkPolyData
RopeMesh::to_vtk_polydata()
{
  IO::CXX::vtkPolyData pd;

  pd.reserve_points(nodes);
  pd.reserve_lines(nodes - 1);

  for (size_t i = 0; i < nodes; i++)
    pd.add_point(centerline[i][0], centerline[i][1], centerline[i][2]);

  for (size_t i = 0; i < nodes - 1; i++)
    pd.add_line(i, i + 1);

  return pd;
}

void
RopeMesh::set_centerline_velocity_zero()
{
  for (size_t i = 0; i < nodes; i++) {
    centerline_velocity[i] = { 0., 0., 0. };
  }
}

void
RopeMesh::set_curvilinear_axis()
{
  for (size_t i = 0; i < nodes; ++i) {
    s[i] = ds * i;
  }
}

void
RopeMesh::set_centerline_x_axis_aligned()
{
  for (size_t i = 0; i < nodes; ++i) {
    centerline[i] = { ds * i, 0., 0. };
    slope[i] = { 1., 0., 0. };
  }
}

} // namespace Models
} // namespace ELFF
