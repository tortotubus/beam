#include "elff/c/models/ibm/IBCoupling.h"
#include <gtest/gtest.h>

#include "elff/models/capsule/CapsuleMeshBuilder.hpp"
#include "elff/models/capsule/IBMCapsule.hpp"
#include "elff/models/runtime/ModelRuntime.hpp"
#include "elff/models/ibm/IBCoupling.hpp"

#include <hdf5.h>

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <vector>

namespace ELFF {
namespace Models {
namespace {

class MockForceModel final : public IBForceCoupled
{
public:
  explicit MockForceModel(int* destruction_count = nullptr, int64_t value = 0)
    : IBForceCoupled(1)
    , destruction_count_(destruction_count)
    , value_(value)
  {}

  ~MockForceModel() override
  {
    if (destruction_count_)
      ++*destruction_count_;
  }

  void pack_state(ModelState& state) const override
  {
    state.ints = {value_};
    state.reals.clear();
    state.bytes.clear();
  }

  void unpack_state(const ModelState& state) override
  {
    if (state.ints.size() != 1)
      throw std::invalid_argument("invalid mock state");
    value_ = state.ints[0];
  }

  int64_t value() const { return value_; }
  void set_value(int64_t value) { value_ = value; }

protected:
  void ComputeNextPoints(std::vector<IBMesh::IBVertex>, real_t) override
  {
    CopyCurrentToNext();
  }

private:
  int* destruction_count_;
  int64_t value_;
};

std::filesystem::path test_output_path(const char* name)
{
  return std::filesystem::path(name);
}

void remove_file_if_present(const std::filesystem::path& path)
{
  std::error_code ec;
  std::filesystem::remove(path, ec);
}

bool hdf5_link_exists(const std::filesystem::path& path,
                      const char* link_path)
{
  hid_t file_id = H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0)
    return false;
  const htri_t exists = H5Lexists(file_id, link_path, H5P_DEFAULT);
  H5Fclose(file_id);
  return exists > 0;
}

TEST(ModelRuntimeTest, OwnsModelsExactlyOnceAndKeepsStableBorrowedViews)
{
  int destructions = 0;
  void* first_handle = nullptr;
  {
    ModelRuntime runtime;
    auto* coupling = runtime.attach_component<IBCoupling>();
    auto first = std::make_unique<MockForceModel>(&destructions, 1);
    first_handle = first.get();
    ASSERT_EQ(coupling->adopt_owned_model(ModelType::EulerBeam,
                                        std::move(first)),
              first_handle);
    auto second = std::make_unique<MockForceModel>(&destructions, 2);
    void* second_handle = second.get();
    ASSERT_EQ(coupling->adopt_owned_model(ModelType::EulerBeamHuang,
                                        std::move(second)),
              second_handle);

    ib_model_view_array_t views = {};
    ASSERT_EQ(coupling->get_model_views(views), 0);
    ASSERT_EQ(views.count, 2U);
    ASSERT_NE(views.data, nullptr);
    const ib_model_view_t* stable_data = views.data;
    EXPECT_EQ(views.data[0].handle, first_handle);
    EXPECT_EQ(views.data[1].handle, second_handle);
    EXPECT_EQ(views.data[0].kind, IB_MODEL_VIEW_FORCE_COUPLED);

    ib_model_view_array_t again = {};
    EXPECT_EQ(coupling->get_model_views(again), 0);
    EXPECT_EQ(again.data, stable_data);
    EXPECT_EQ(again.data[0].handle, first_handle);

    auto rejected = std::make_unique<MockForceModel>(&destructions, 3);
    EXPECT_EQ(coupling->adopt_owned_model(ModelType::EulerBeam,
                                        std::move(rejected)),
              nullptr);
    EXPECT_EQ(destructions, 1);
  }
  EXPECT_EQ(destructions, 3);
}

TEST(ModelRuntimeTest, ReturnsOneThousandViewsInAdoptionOrder)
{
  ModelRuntime runtime;
    auto* coupling = runtime.attach_component<IBCoupling>();
  std::vector<void*> expected;
  expected.reserve(1000);
  for (size_t i = 0; i < 1000; ++i) {
    auto model = std::make_unique<MockForceModel>(nullptr, i);
    expected.push_back(model.get());
    ASSERT_NE(coupling->adopt_owned_model(ModelType::EulerBeam,
                                        std::move(model)),
              nullptr);
  }
  ib_model_view_array_t views = {};
  ASSERT_EQ(coupling->get_model_views(views), 0);
  ASSERT_EQ(views.count, expected.size());
  for (size_t i = 0; i < views.count; ++i)
    EXPECT_EQ(views.data[i].handle, expected[i]);
}

TEST(ModelRuntimeTest, CViewAPIHandlesNullAndSealsManualRuntime)
{
  ib_model_view_array_t views = {17, reinterpret_cast<ib_model_view_t*>(1)};
  EXPECT_NE(ib_coupling_get_model_views(ib_coupling_attach(nullptr), &views), 0);
  EXPECT_EQ(views.count, 0U);
  EXPECT_EQ(views.data, nullptr);

  model_runtime_t runtime = model_runtime_new();
  ASSERT_NE(runtime, nullptr);
  EXPECT_EQ(ib_coupling_get_model_views(ib_coupling_attach(runtime), &views), 0);
  EXPECT_EQ(views.count, 0U);
  EXPECT_EQ(views.data, nullptr);
  model_runtime_delete(runtime);
}

TEST(ModelRuntimeTest, DefaultModelPolyDataHookIsNoOp)
{
  ModelRuntime runtime;
  ASSERT_NE(runtime.adopt_owned_model(
              ModelType::EulerBeam,
              std::make_unique<MockForceModel>()),
            nullptr);
  IO::CXX::vtkPolyData pd = runtime.to_vtk_polydata();
  EXPECT_EQ(pd.number_of_points(), 0);
  EXPECT_EQ(pd.number_of_cells(), 0);
  EXPECT_EQ(pd.number_of_pointdata(), 0);
  EXPECT_EQ(pd.number_of_celldata(), 0);
}

TEST(ModelRuntimeTest, WritesOwnedCapsulePolyData)
{
  const std::filesystem::path output_path =
    test_output_path("elff_model_runtime_capsule_polydata.vtkhdf");
  remove_file_if_present(output_path);

  ModelRuntime runtime;
  auto capsule = std::make_unique<IBMCapsule>(
    CapsuleMeshBuilder::sphere({1.0, Vec3::Zero(), 1}));
  ASSERT_NE(runtime.adopt_owned_model(ModelType::Capsule,
                                      std::move(capsule)),
            nullptr);

  IO::CXX::vtkPolyData pd = runtime.to_vtk_polydata();
  ASSERT_EQ(pd.number_of_points(), 42);
  ASSERT_EQ(pd.number_of_polygons(), 80);
  ASSERT_EQ(pd.number_of_celldata(), 4);

  EXPECT_EQ(runtime.write_polydata(output_path.string().c_str(), true), 0);
  ASSERT_TRUE(std::filesystem::exists(output_path));
  EXPECT_GT(std::filesystem::file_size(output_path), 0);
  EXPECT_TRUE(hdf5_link_exists(output_path, "/VTKHDF/CellData/triangle_area"));
  remove_file_if_present(output_path);
}

TEST(ModelRuntimeTest, CheckpointV2RoundTripsOrderedMultipleModels)
{
  const auto path = test_output_path("elff_runtime_v2.checkpoint");
  remove_file_if_present(path);
  ModelRuntime runtime;
  auto first = std::make_unique<MockForceModel>(nullptr, 11);
  auto* first_view = first.get();
  auto second = std::make_unique<MockForceModel>(nullptr, 29);
  auto* second_view = second.get();
  ASSERT_NE(runtime.adopt_owned_model(ModelType::EulerBeam,
                                      std::move(first)),
            nullptr);
  ASSERT_NE(runtime.adopt_owned_model(ModelType::EulerBeamHuang,
                                      std::move(second)),
            nullptr);
  ASSERT_EQ(runtime.write_checkpoint(path.string().c_str()), 0);
  first_view->set_value(-1);
  second_view->set_value(-1);
  ASSERT_EQ(runtime.read_checkpoint(path.string().c_str()), 0);
  EXPECT_EQ(first_view->value(), 11);
  EXPECT_EQ(second_view->value(), 29);
  remove_file_if_present(path);
}

TEST(ModelRuntimeTest, RejectsCountAndTypeSequenceBeforeMutation)
{
  ModelRuntime runtime;
  auto first = std::make_unique<MockForceModel>(nullptr, 1);
  auto* first_view = first.get();
  auto second = std::make_unique<MockForceModel>(nullptr, 2);
  auto* second_view = second.get();
  ASSERT_NE(runtime.adopt_owned_model(ModelType::EulerBeam,
                                      std::move(first)),
            nullptr);
  ASSERT_NE(runtime.adopt_owned_model(ModelType::EulerBeamHuang,
                                      std::move(second)),
            nullptr);

  RuntimeState state;
  ASSERT_EQ(runtime.export_state(state), 0);
  state.models[0].state.ints[0] = 101;
  state.models[1].state.ints[0] = 202;
  std::swap(state.models[0].type, state.models[1].type);
  EXPECT_NE(runtime.import_state(state), 0);
  EXPECT_EQ(first_view->value(), 1);
  EXPECT_EQ(second_view->value(), 2);
  state.models.pop_back();
  EXPECT_NE(runtime.import_state(state), 0);
}

TEST(ModelRuntimeTest, RejectsV1AndTruncatedCheckpointFiles)
{
  const auto v1_path = test_output_path("elff_runtime_v1.checkpoint");
  const auto truncated_path =
    test_output_path("elff_runtime_truncated.checkpoint");
  remove_file_if_present(v1_path);
  remove_file_if_present(truncated_path);

  struct LegacyHeader
  {
    uint64_t magic;
    int32_t version;
    int64_t nmodels;
  };
  const LegacyHeader legacy = {0x494252554E54494DULL, 1, 0};
  {
    std::ofstream out(v1_path, std::ios::binary);
    out.write(reinterpret_cast<const char*>(&legacy), sizeof(legacy));
  }
  {
    std::ofstream out(truncated_path, std::ios::binary);
    out.write(reinterpret_cast<const char*>(&legacy), 5);
  }

  ModelRuntime runtime;
  EXPECT_NE(runtime.read_checkpoint(v1_path.string().c_str()), 0);
  EXPECT_NE(runtime.read_checkpoint(truncated_path.string().c_str()), 0);
  remove_file_if_present(v1_path);
  remove_file_if_present(truncated_path);
}

} // namespace
} // namespace Models
} // namespace ELFF
