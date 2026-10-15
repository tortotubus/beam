#include <gtest/gtest.h>

#include "elff/models/capsule/CapsuleMeshBuilder.hpp"
#include "elff/models/capsule/CapsuleMeshVTK.hpp"

#include <elff/io/CXX/vtkHDFPolyData.hpp>
#include <elff/io/CXX/vtkPolyData.hpp>

#include <hdf5.h>

#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

namespace ELFF {
namespace Models {
namespace {

std::filesystem::path test_output_path(const char *name) {
  return name;
}

void remove_file_if_present(const std::filesystem::path &path) {
  std::error_code ec;
  std::filesystem::remove(path, ec);
}

bool hdf5_link_exists(const std::filesystem::path &path,
                      const char *linkPath) {
  hid_t fileId = H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (fileId < 0)
    return false;

  const htri_t exists = H5Lexists(fileId, linkPath, H5P_DEFAULT);
  H5Fclose(fileId);
  return exists > 0;
}

void expect_static_capsule_plot_file(const std::string &name,
                                     const CapsuleMesh &mesh) {
  const std::filesystem::path outputPath = test_output_path(name.c_str());
  remove_file_if_present(outputPath);

  {
    IO::CXX::vtkPolyData pd = CapsuleMeshVTK::to_vtk_polydata(mesh);
    IO::CXX::vtkHDFPolyData writer(outputPath.string(), pd);
    writer.write_new_static(true);
  }

  const bool fileExists = std::filesystem::exists(outputPath);
  EXPECT_TRUE(fileExists);
  if (fileExists)
    EXPECT_GT(std::filesystem::file_size(outputPath), 0);
  EXPECT_TRUE(hdf5_link_exists(outputPath, "/VTKHDF/CellData/triangle_area"));
  EXPECT_TRUE(hdf5_link_exists(outputPath, "/VTKHDF/CellData/triangle_normal"));
  EXPECT_TRUE(hdf5_link_exists(outputPath, "/VTKHDF/CellData/stretch"));
  EXPECT_TRUE(hdf5_link_exists(outputPath, "/VTKHDF/CellData/tension"));

  // remove_file_if_present(outputPath);
}

TEST(CapsuleMeshVTKTest, ExportsSphereBuilderAsSurfacePolyDataWithNodeFields) {
  CapsuleMesh mesh = CapsuleMeshBuilder::sphere({ 1.0, Vec3::Zero(), 1 });

  for (int i = 0; i < mesh.numNodes(); ++i) {
    mesh.state.v.col(i) = Vec3(i, i + 1, i + 2);
    mesh.state.f.col(i) = Vec3(2 * i, 2 * i + 1, 2 * i + 2);
    mesh.state.meanCurv(i) = 0.1 * i;
    mesh.state.refCurv(i) = 0.2 * i;
    mesh.state.gaussCurv(i) = 0.3 * i;
  }
  for (int tid = 0; tid < mesh.numTriangles(); ++tid) {
    mesh.state.triGeom[static_cast<size_t>(tid)].stretch =
        Vec2(1.0 + 0.01 * tid, 1.5 + 0.01 * tid);
    mesh.state.triGeom[static_cast<size_t>(tid)].tension =
        Vec2(2.0 + 0.02 * tid, 2.5 + 0.02 * tid);
  }

  IO::CXX::vtkPolyData pd = CapsuleMeshVTK::to_vtk_polydata(mesh);

  ASSERT_EQ(pd.number_of_points(), 42);
  ASSERT_EQ(pd.number_of_polygons(), 80);
  ASSERT_EQ(pd.number_of_cells(), 80);
  ASSERT_EQ(pd.number_of_pointdata(), 6);
  ASSERT_EQ(pd.number_of_celldata(), 4);

  const std::vector<double> &velocity = pd.get_pointdata(0);
  const std::vector<double> &force = pd.get_pointdata(1);
  const std::vector<double> &normal = pd.get_pointdata(2);
  const std::vector<double> &meanCurvature = pd.get_pointdata(3);
  const std::vector<double> &referenceCurvature = pd.get_pointdata(4);
  const std::vector<double> &gaussianCurvature = pd.get_pointdata(5);
  const std::vector<double> &triangleArea = pd.get_celldata(0);
  const std::vector<double> &triangleNormal = pd.get_celldata(1);
  const std::vector<double> &stretch = pd.get_celldata(2);
  const std::vector<double> &tension = pd.get_celldata(3);

  EXPECT_EQ(velocity.size(), 3 * pd.number_of_points());
  EXPECT_EQ(force.size(), 3 * pd.number_of_points());
  EXPECT_EQ(normal.size(), 3 * pd.number_of_points());
  EXPECT_EQ(meanCurvature.size(), pd.number_of_points());
  EXPECT_EQ(referenceCurvature.size(), pd.number_of_points());
  EXPECT_EQ(gaussianCurvature.size(), pd.number_of_points());
  EXPECT_EQ(triangleArea.size(), pd.number_of_cells());
  EXPECT_EQ(triangleNormal.size(), 3 * pd.number_of_cells());
  EXPECT_EQ(stretch.size(), 2 * pd.number_of_cells());
  EXPECT_EQ(tension.size(), 2 * pd.number_of_cells());

  EXPECT_DOUBLE_EQ(velocity[0], 0.0);
  EXPECT_DOUBLE_EQ(velocity[1], 1.0);
  EXPECT_DOUBLE_EQ(velocity[2], 2.0);
  EXPECT_DOUBLE_EQ(force[3], 2.0);
  EXPECT_DOUBLE_EQ(force[4], 3.0);
  EXPECT_DOUBLE_EQ(force[5], 4.0);
  EXPECT_NEAR(normal[0] * normal[0] + normal[1] * normal[1] +
                normal[2] * normal[2],
              1.0,
              1e-12);
  EXPECT_DOUBLE_EQ(meanCurvature[10], 1.0);
  EXPECT_DOUBLE_EQ(referenceCurvature[10], 2.0);
  EXPECT_DOUBLE_EQ(gaussianCurvature[10], 3.0);

  EXPECT_GT(triangleArea[0], 0.0);
  EXPECT_TRUE(std::isfinite(triangleNormal[0]));
  EXPECT_TRUE(std::isfinite(triangleNormal[1]));
  EXPECT_TRUE(std::isfinite(triangleNormal[2]));
  EXPECT_NEAR(triangleNormal[0] * triangleNormal[0] +
                triangleNormal[1] * triangleNormal[1] +
                triangleNormal[2] * triangleNormal[2],
              1.0,
              1e-12);
  EXPECT_DOUBLE_EQ(stretch[0], 1.0);
  EXPECT_DOUBLE_EQ(stretch[1], 1.5);
  EXPECT_DOUBLE_EQ(tension[0], 2.0);
  EXPECT_DOUBLE_EQ(tension[1], 2.5);
}

TEST(CapsuleMeshVTKTest, WritesStaticPolyDataForCapsuleBuilders) {
  expect_static_capsule_plot_file(
      "elff_capsule_sphere_builder_static.vtkhdf",
      CapsuleMeshBuilder::sphere({ 1.0, Vec3(0.0, 0.0, 0.0), 4 }));

  expect_static_capsule_plot_file(
      "elff_capsule_ellipsoid_builder_static.vtkhdf",
      CapsuleMeshBuilder::ellipsoid(
          { Vec3(1.5, 0.75, 0.5), Vec3(0.25, -0.5, 0.75), 4 }));

  expect_static_capsule_plot_file(
      "elff_capsule_biconcave_builder_static.vtkhdf",
      CapsuleMeshBuilder::biconcave({ 1.0, Vec3(0.0, 0.0, 0.0), 4 }));
}

} // namespace
} // namespace Models
} // namespace ELFF
