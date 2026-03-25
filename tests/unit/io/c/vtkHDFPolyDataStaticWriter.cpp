#include <gtest/gtest.h>

#include <elff/io/C/vtkHDFPolyData.hpp>
#include <elff/io/C/vtkPolyData.hpp>

#include <cstdio>
#include <filesystem>
#include <string>

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

TEST(vtkHDFPolyDataStaticWriterTest, WritesGeometryOnlyFile)
{
  const std::filesystem::path output_path =
    test_output_path("elff_vtk_polydata_static_geometry.vtkhdf");
  remove_file_if_present(output_path);

  vtkPolyData pd = vtk_polydata_init(3, 1, 2, 0, 0);

  vtk_polydata_add_point(&pd, 0.f, 0.f, 0.f);
  vtk_polydata_add_point(&pd, 1.f, 0.f, 0.f);
  vtk_polydata_add_point(&pd, 1.f, 1.f, 0.f);

  vtk_polydata_add_vertex(&pd, 0);
  vtk_polydata_add_line(&pd, 1, 2);

  vtkHDFPolyData writer =
    vtk_HDF_polydata_init_static(output_path.string().c_str(), true, &pd);
  vtk_HDF_polydata_close(&writer);

  EXPECT_TRUE(std::filesystem::exists(output_path));
  EXPECT_GT(std::filesystem::file_size(output_path), 0);

  vtk_polydata_free(&pd);
}

TEST(vtkHDFPolyDataStaticWriterTest, WritesFileWithPointAndCellData)
{
  const std::filesystem::path output_path =
    test_output_path("elff_vtk_polydata_static_fields.vtkhdf");
  remove_file_if_present(output_path);

  vtkPolyData pd = vtk_polydata_init(3, 1, 2, 0, 0, 1, 1);

  vtk_polydata_add_point(&pd, 0.f, 0.f, 0.f);
  vtk_polydata_add_point(&pd, 1.f, 0.f, 0.f);
  vtk_polydata_add_point(&pd, 1.f, 1.f, 0.f);

  vtk_polydata_add_vertex(&pd, 0);
  vtk_polydata_add_line(&pd, 1, 2);

  const int64_t point_field = vtk_polydata_add_pointdata_scalar(&pd, "point_id");
  const int64_t cell_field = vtk_polydata_add_celldata_vector(&pd, "cell_vec", 3);

  double* point_data = vtk_polydata_get_pointdata_data(&pd, point_field);
  point_data[0] = 1.;
  point_data[1] = 2.;
  point_data[2] = 3.;

  double* cell_data = vtk_polydata_get_celldata_data(&pd, cell_field);
  cell_data[0] = 1.;
  cell_data[1] = 0.;
  cell_data[2] = 0.;
  cell_data[3] = 0.;
  cell_data[4] = 1.;
  cell_data[5] = 0.;

  vtkHDFPolyData writer =
    vtk_HDF_polydata_init_static(output_path.string().c_str(), true, &pd);
  vtk_HDF_polydata_close(&writer);

  EXPECT_TRUE(std::filesystem::exists(output_path));
  EXPECT_GT(std::filesystem::file_size(output_path), 0);

  vtk_polydata_free(&pd);
}

} // namespace
} // namespace C
} // namespace IO
} // namespace ELFF
