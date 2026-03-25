#include <gtest/gtest.h>

#include <elff/io/C/vtkHDFPolyData.hpp>
#include <elff/io/C/vtkPolyData.hpp>

#include <filesystem>

namespace ELFF {
namespace IO {
namespace C {
namespace {

std::filesystem::path
test_output_path(const char* name)
{
  return std::filesystem::path(name);
}

void
remove_file_if_present(const std::filesystem::path& path)
{
  std::error_code ec;
  std::filesystem::remove(path, ec);
}

vtkPolyData
make_polydata(float x_offset, double cell_value_offset)
{
  vtkPolyData pd = vtk_polydata_init(3, 1, 2, 0, 0, 1, 1);

  vtk_polydata_add_point(&pd, x_offset + 0.f, 0.f, 0.f);
  vtk_polydata_add_point(&pd, x_offset + 1.f, 0.f, 0.f);
  vtk_polydata_add_point(&pd, x_offset + 1.f, 1.f, 0.f);

  vtk_polydata_add_vertex(&pd, 0);
  vtk_polydata_add_line(&pd, 1, 2);

  const int64_t point_field = vtk_polydata_add_pointdata_vector(&pd, "velocity", 3);
  double* point_data = vtk_polydata_get_pointdata_data(&pd, point_field);
  point_data[0] = x_offset + 1.;
  point_data[1] = x_offset + 2.;
  point_data[2] = x_offset + 3.;
  point_data[3] = x_offset + 4.;
  point_data[4] = x_offset + 5.;
  point_data[5] = x_offset + 6.;
  point_data[6] = x_offset + 7.;
  point_data[7] = x_offset + 8.;
  point_data[8] = x_offset + 9.;

  const int64_t cell_field = vtk_polydata_add_celldata_scalar(&pd, "cell_id");
  double* cell_data = vtk_polydata_get_celldata_data(&pd, cell_field);
  cell_data[0] = cell_value_offset + 1.;
  cell_data[1] = cell_value_offset + 2.;

  return pd;
}

TEST(vtkHDFPolyDataTransientWriterTest, WritesNewTransientFile)
{
  const std::filesystem::path output_path =
    test_output_path("elff_vtk_polydata_transient_geometry.vtkhdf");
  remove_file_if_present(output_path);

  vtkPolyData pd = make_polydata(0.f, 0.);

  vtkHDFPolyData writer =
    vtk_HDF_polydata_init_transient(output_path.string().c_str(), true, &pd, 0.f);
  vtk_HDF_polydata_close(&writer);

  EXPECT_TRUE(std::filesystem::exists(output_path));
  EXPECT_GT(std::filesystem::file_size(output_path), 0);

  vtk_polydata_free(&pd);
}

TEST(vtkHDFPolyDataTransientWriterTest, AppendsTransientFile)
{
  const std::filesystem::path output_path =
    test_output_path("elff_vtk_polydata_transient_append.vtkhdf");
  remove_file_if_present(output_path);

  vtkPolyData pd0 = make_polydata(0.f, 0.);
  vtkHDFPolyData writer0 =
    vtk_HDF_polydata_init_transient(output_path.string().c_str(), true, &pd0, 0.f);
  vtk_HDF_polydata_close(&writer0);

  ASSERT_TRUE(std::filesystem::exists(output_path));
  ASSERT_GT(std::filesystem::file_size(output_path), 0);

  vtkPolyData pd1 = make_polydata(10.f, 10.);
  vtkHDFPolyData writer1 =
    vtk_HDF_polydata_append_transient(output_path.string().c_str(), &pd1, 1.f);
  vtk_HDF_polydata_close(&writer1);

  EXPECT_TRUE(std::filesystem::exists(output_path));
  EXPECT_GT(std::filesystem::file_size(output_path), 0);

  vtk_polydata_free(&pd0);
  vtk_polydata_free(&pd1);
}

} // namespace
} // namespace C
} // namespace IO
} // namespace ELFF
