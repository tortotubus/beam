// vtkPolyDataTest.cpp
#include <cmath>
#include <cstdint>
#include <gtest/gtest.h>

// Adjust include path / name as needed:
#include "io/C/vtkHDFPolyData.hpp"

namespace ELFF {
namespace io {
namespace C {

TEST(VtkHDFPolyDataStaticTest, WriteLine)
{
  size_t n_points = 5;
  size_t n_points_vertices = 5;
  size_t n_points_lines = 5;
  size_t n_points_strips = 0;
  size_t n_points_polygons = 0;

  double ds = 0.1;

  vtkPolyData pd = vtk_polydata_init(n_points, 100, 100, 100, 100);

  for (size_t i = 0; i < n_points; i++) {
    vtk_polydata_add_point(&pd, 0 + ds * i, 0, 0);
  }

  for (size_t i = 0; i < n_points_lines - 1; i++) {
    vtk_polydata_add_line(&pd, i, i + 1);
  }

  for (size_t i = 0; i < n_points_vertices; i++) {
    vtk_polydata_add_vertex(&pd, i);
  }

  vtkHDFPolyData hdf_pd =
    vtk_HDF_polydata_init_static("polydata_static.vtkhdf", true, &pd);
  vtk_HDF_polydata_close(&hdf_pd);

  vtk_polydata_free(&pd);
}

TEST(VtkHDFPolyDataTransientTest, WriteLine)
{
  size_t n_points = 10;

  double ds = 0.1, dt = 0.1, dtheta = M_PI * 2 * 0.01;
  const char* fname = "polydata_transient.vtkhdf";

  vtkPolyData pd = vtk_polydata_init(n_points, 100, 100, 100, 100);

  for (size_t i = 0; i < n_points; i++) {
    double x = (ds * i);
    double y = 0;
    double z = 0;
    vtk_polydata_add_point(&pd, x, y, z);
  }

  for (size_t i = 0; i < n_points - 1; i++) {
    vtk_polydata_add_line(&pd, i, i + 1);
  }

  for (size_t i = 0; i < n_points; i++) {
    vtk_polydata_add_vertex(&pd, i);
  }

  // Write the initial file
  vtkHDFPolyData hdf_pd = vtk_HDF_polydata_init_transient(fname, true, &pd, 0.0);
  vtk_HDF_polydata_close(&hdf_pd);

  for (size_t ti = 1; ti < 100; ti++) {
    // Modify the points but not connectivity
    for (size_t i = 0; i < n_points; i++) {
      pd.points[3 * i + 0] = (ds * i) * cos(ti * dtheta); // x
      pd.points[3 * i + 1] = (ds * i) * sin(ti * dtheta); // y
      pd.points[3 * i + 2] = 0;                           // z
    }

    // Write the second time step
    hdf_pd = vtk_HDF_polydata_append_transient(fname, &pd, ti*dt);
    vtk_HDF_polydata_close(&hdf_pd);
  }

  vtk_polydata_free(&pd);
}

}
}
}