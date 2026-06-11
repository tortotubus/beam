#include <gtest/gtest.h>

#include "elff/models/capsule/CapsuleMeshBuilder.hpp"
#include "elff/models/capsule/IBMCapsule.hpp"
#include "elff/models/ibm/IBRuntime.hpp"

#include <hdf5.h>

#include <filesystem>

namespace ELFF {
namespace Models {
namespace {

class NoPolyDataModel final : public IBModel {
public:
  void pack_state(IBModelState& state) const override
  {
    state.ints.clear();
    state.reals.clear();
    state.bytes.clear();
  }

  void unpack_state(const IBModelState& state) override
  {
    (void)state;
  }
};

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

bool
hdf5_link_exists(const std::filesystem::path& path, const char* link_path)
{
  hid_t file_id = H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0)
    return false;

  const htri_t exists = H5Lexists(file_id, link_path, H5P_DEFAULT);
  H5Fclose(file_id);
  return exists > 0;
}

TEST(IBRuntimeTest, DefaultModelPolyDataHookIsNoOp)
{
  IBRuntime runtime;
  NoPolyDataModel model;

  runtime.register_model(model);

  IO::CXX::vtkPolyData pd = runtime.to_vtk_polydata();

  EXPECT_EQ(pd.number_of_points(), 0);
  EXPECT_EQ(pd.number_of_cells(), 0);
  EXPECT_EQ(pd.number_of_pointdata(), 0);
  EXPECT_EQ(pd.number_of_celldata(), 0);
}

TEST(IBRuntimeTest, WritesRegisteredCapsulePolyData)
{
  const std::filesystem::path output_path =
    test_output_path("elff_ib_runtime_capsule_polydata.vtkhdf");
  remove_file_if_present(output_path);

  IBRuntime runtime;
  IBMCapsule capsule(CapsuleMeshBuilder::sphere({ 1.0, Vec3::Zero(), 1 }));

  runtime.register_model(capsule);

  IO::CXX::vtkPolyData pd = runtime.to_vtk_polydata();
  ASSERT_EQ(pd.number_of_points(), 42);
  ASSERT_EQ(pd.number_of_polygons(), 80);
  ASSERT_EQ(pd.number_of_celldata(), 4);

  EXPECT_EQ(runtime.write_polydata(output_path.string().c_str(), true), 0);
  ASSERT_TRUE(std::filesystem::exists(output_path));
  EXPECT_GT(std::filesystem::file_size(output_path), 0);
  EXPECT_TRUE(hdf5_link_exists(output_path, "/VTKHDF/CellData/triangle_area"));
  EXPECT_TRUE(hdf5_link_exists(output_path, "/VTKHDF/CellData/triangle_normal"));
}

} // namespace
} // namespace Models
} // namespace ELFF
