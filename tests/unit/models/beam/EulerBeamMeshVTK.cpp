#include <gtest/gtest.h>

#include <elff/io/CXX/vtkHDFPolyData.hpp>
#include <elff/io/CXX/vtkPolyData.hpp>
#include <elff/models/beam/EulerBeamMesh.hpp>

#include <filesystem>

namespace ELFF {
namespace Models {
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

TEST(EulerBeamMeshVTKTest, ExportsCenterlineVelocityAsPointData)
{
  EulerBeamMesh mesh(3, 2.0);

  mesh.get_centerline_velocity(0) = { 1.0, 2.0, 3.0 };
  mesh.get_centerline_velocity(1) = { 4.0, 5.0, 6.0 };
  mesh.get_centerline_velocity(2) = { 7.0, 8.0, 9.0 };

  IO::CXX::vtkPolyData pd = mesh.to_vtk_polydata();

  ASSERT_EQ(pd.number_of_points(), 3);
  ASSERT_EQ(pd.number_of_lines(), 2);
  ASSERT_EQ(pd.number_of_pointdata(), 1);

  const std::vector<double>& velocity_data = pd.get_pointdata(0);
  const std::vector<double> expected = { 1.0, 2.0, 3.0,
                                         4.0, 5.0, 6.0,
                                         7.0, 8.0, 9.0 };

  EXPECT_EQ(velocity_data, expected);
}

TEST(EulerBeamMeshVTKTest, WritesStaticPolyDataWithCenterlineVelocity)
{
  const std::filesystem::path output_path =
    test_output_path("elff_euler_beam_mesh_static.vtkhdf");
  remove_file_if_present(output_path);

  EulerBeamMesh mesh(4, 3.0);
  mesh.get_centerline_velocity(0) = { 0.0, 1.0, 0.0 };
  mesh.get_centerline_velocity(1) = { 1.0, 0.0, 0.0 };
  mesh.get_centerline_velocity(2) = { 0.0, -1.0, 0.0 };
  mesh.get_centerline_velocity(3) = { -1.0, 0.0, 0.0 };

  IO::CXX::vtkPolyData pd = mesh.to_vtk_polydata();
  IO::CXX::vtkHDFPolyData writer(output_path.string(), pd);
  writer.write_new_static(true);

  EXPECT_TRUE(std::filesystem::exists(output_path));
  EXPECT_GT(std::filesystem::file_size(output_path), 0);
}

TEST(EulerBeamMeshVTKTest, WritesTransientPolyDataWithCenterlineVelocity)
{
  const std::filesystem::path output_path =
    test_output_path("elff_euler_beam_mesh_transient.vtkhdf");
  remove_file_if_present(output_path);

  EulerBeamMesh mesh(4, 3.0);
  mesh.get_centerline_velocity(0) = { 0.0, 1.0, 0.0 };
  mesh.get_centerline_velocity(1) = { 1.0, 0.0, 0.0 };
  mesh.get_centerline_velocity(2) = { 0.0, -1.0, 0.0 };
  mesh.get_centerline_velocity(3) = { -1.0, 0.0, 0.0 };

  IO::CXX::vtkPolyData pd = mesh.to_vtk_polydata();
  IO::CXX::vtkHDFPolyData writer(output_path.string(), pd);
  writer.write_new_transient(true, 0.f);

  EXPECT_TRUE(std::filesystem::exists(output_path));
  EXPECT_GT(std::filesystem::file_size(output_path), 0);
}

} // namespace
} // namespace Models
} // namespace ELFF
