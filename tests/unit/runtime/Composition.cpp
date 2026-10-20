#include "elff/models/ModelRuntime.hpp"
#include "elff/models/ibm/IBCoupling.hpp"
#include "elff/models/articulated/IBBendingRodImplicit.hpp"
#include <gtest/gtest.h>
#include <type_traits>

using namespace ELFF::Models;
namespace {
class PlainModel : public Model {
public:
  int value = 7;
  void pack_state(ModelState& s) const override { s = {}; s.ints = {value}; }
  void unpack_state(const ModelState& s) override { value = static_cast<int>(s.ints.at(0)); }
};
class ForceModel : public IBForceCoupled {
public:
  ForceModel() : IBForceCoupled(1) {}
  void pack_state(ModelState& s) const override { s = {}; }
  void unpack_state(const ModelState&) override {}
protected:
  void ComputeNextPoints(std::vector<IBMesh::IBVertex>, ELFF::real_t) override {}
};
class UncoupledIB : public IBModel {
public:
  void pack_state(ModelState&) const override {}
  void unpack_state(const ModelState&) override {}
};
class OtherComponent : public RuntimeComponent {
public:
  using RuntimeComponent::RuntimeComponent;
  int seal() override { return 0; }
};
class BothCouplings : public IBForceCoupled, public IBVelocityCoupled {
public:
  BothCouplings() : IBForceCoupled(1), IBVelocityCoupled(1) {}
  void pack_state(ModelState& s) const override { s = {}; }
  void unpack_state(const ModelState&) override {}
protected:
  void ComputeNextPoints(std::vector<IBMesh::IBVertex>, ELFF::real_t) override {}
  void ComputeMidpointForces() override {}
};

TEST(RuntimeCompositionTest, PlainModelsNeedNoIBComponent)
{
  ModelRuntime runtime;
  EXPECT_EQ(runtime.component<IBCoupling>(), nullptr);
  auto* model = runtime.adopt_owned_model(static_cast<ModelType>(101), std::make_unique<PlainModel>());
  ASSERT_NE(model, nullptr);
  EXPECT_EQ(runtime.model(0), model);
  EXPECT_EQ(runtime.entry(0)->owner_pid, 0);
  EXPECT_EQ(runtime.entry(99), nullptr);
  RuntimeState snapshot;
  ASSERT_EQ(runtime.export_state(snapshot), 0);
  model->value = 12;
  ASSERT_EQ(runtime.import_state(snapshot), 0);
  EXPECT_EQ(model->value, 7);
  EXPECT_EQ(runtime.component<IBCoupling>(), nullptr);
  EXPECT_EQ(runtime.attach_component<IBCoupling>(), nullptr);
}
TEST(RuntimeCompositionTest, ExplicitSubsetAndOrdering)
{
  ModelRuntime runtime;
  auto* coupling = runtime.attach_component<IBCoupling>();
  ASSERT_NE(coupling, nullptr);
  EXPECT_EQ(coupling, runtime.attach_component<IBCoupling>());
  runtime.adopt_owned_model(static_cast<ModelType>(101), std::make_unique<PlainModel>());
  auto* first = runtime.adopt_owned_model(ModelType::EulerBeam, std::make_unique<ForceModel>());
  runtime.adopt_owned_model(ModelType::EulerBeam, std::make_unique<ForceModel>()); // unregistered
  auto* last = runtime.adopt_owned_model(ModelType::EulerBeam, std::make_unique<ForceModel>());
  ASSERT_EQ(coupling->register_model(3), 0);
  ASSERT_EQ(coupling->register_model(1), 0);
  EXPECT_NE(coupling->register_model(1), 0);
  EXPECT_NE(coupling->register_model(0), 0);
  EXPECT_NE(coupling->register_model(99), 0);
  EXPECT_EQ(coupling->adopt_owned_model(ModelType::EulerBeam, std::make_unique<UncoupledIB>()), nullptr);
  EXPECT_EQ(coupling->adopt_owned_model(ModelType::Invalid, std::make_unique<ForceModel>()), nullptr);
  EXPECT_EQ(coupling->adopt_owned_model(ModelType::EulerBeam,
      std::unique_ptr<IBForceCoupled>(new BothCouplings())), nullptr);
  EXPECT_EQ(runtime.size(), 4u);
  ib_model_view_array_t views{};
  ASSERT_EQ(coupling->get_model_views(views), 0);
  ASSERT_EQ(views.count, 2u);
  EXPECT_EQ(views.data[0].handle, first);
  EXPECT_EQ(views.data[1].handle, last);
  const auto* stable = views.data;
  EXPECT_EQ(runtime.seal(), 0);
  ASSERT_EQ(coupling->get_model_views(views), 0);
  EXPECT_EQ(views.data, stable);
  EXPECT_NE(coupling->register_model(2), 0);
  EXPECT_EQ(runtime.attach_component<OtherComponent>(), nullptr);
  EXPECT_EQ(runtime.attach_component<IBCoupling>(), coupling);
  EXPECT_EQ(coupling->adopt_owned_model(ModelType::EulerBeam, std::make_unique<ForceModel>()), nullptr);
  EXPECT_EQ(runtime.size(), 4u);
}
class LifetimeModel : public PlainModel {
public:
  explicit LifetimeModel(std::vector<int>& order) : order_(order) {}
  ~LifetimeModel() override { order_.push_back(2); }
private:
  std::vector<int>& order_;
};
class LifetimeComponent : public RuntimeComponent {
public:
  using RuntimeComponent::RuntimeComponent;
  std::vector<int>* order = nullptr;
  ~LifetimeComponent() override {
    EXPECT_NE(runtime_.model(0), nullptr);
    if (order) order->push_back(1);
  }
  int seal() override { return 0; }
};
TEST(RuntimeCompositionTest, ComponentsDieBeforeModels)
{
  std::vector<int> order;
  {
    auto runtime = std::make_unique<ModelRuntime>();
    runtime->adopt_owned_model(static_cast<ModelType>(101),std::make_unique<LifetimeModel>(order));
    runtime->attach_component<LifetimeComponent>()->order = &order;
  }
  EXPECT_EQ(order,(std::vector<int>{1,2}));
}
TEST(RuntimeCompositionTest, CAttachmentAndLifetime)
{
  auto runtime = model_runtime_new();
  ASSERT_NE(runtime, nullptr);
  EXPECT_EQ(ib_coupling_get(runtime), nullptr);
  auto coupling = ib_coupling_attach(runtime);
  ASSERT_NE(coupling, nullptr);
  EXPECT_EQ(ib_coupling_get(runtime), coupling);
  EXPECT_EQ(ib_coupling_attach(runtime), coupling);
  ib_model_view_array_t views{};
  ASSERT_EQ(ib_coupling_get_model_views(coupling, &views), 0);
  EXPECT_EQ(views.count, 0u);
  model_runtime_delete(runtime);
  EXPECT_EQ(ib_coupling_attach(nullptr), nullptr);
}
TEST(RuntimeCompositionTest, LegacyRodCheckpointContinuesWithImplicitHistory)
{
  if (sizeof(ELFF::real_t) != 8) GTEST_SKIP() << "Legacy fixture uses double precision.";
  BendingRodImplicitConfig config;
  config.links = 2;
  ModelRuntime runtime;
  auto* coupling = runtime.attach_component<IBCoupling>();
  auto* rod = coupling->adopt_owned_model(ModelType::BendingRodImplicit,
                                        std::make_unique<IBBendingRodImplicit>(config));
  IBBendingRodImplicit direct(config);
  direct.solve(0.01);
  ASSERT_EQ(runtime.read_checkpoint(ELFF_LEGACY_RUNTIME_CHECKPOINT), 0);
  ModelState expected, restored;
  direct.pack_state(expected); rod->pack_state(restored);
  EXPECT_EQ(restored.ints,expected.ints); EXPECT_EQ(restored.reals,expected.reals);
  direct.solve(0.01); rod->solve(0.01);
  direct.pack_state(expected); rod->pack_state(restored);
  EXPECT_EQ(restored.reals,expected.reals);
}
} // namespace
