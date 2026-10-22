#include "elff/models/runtime/ModelInsertion.hpp"
#include "elff/models/ibm/IBCoupling.hpp"
#include "elff/models/articulated/IBBendingRodImplicit.hpp"
#include "elff/models/beam/IBEulerBeamHuang.hpp"
#include "elff/contact/NarrowPhaseSegmentDistance.hpp"
#include <gtest/gtest.h>
#include <random>

using namespace ELFF;
using namespace ELFF::Models;
using namespace ELFF::Contact;
namespace {
using Point = std::array<real_t,3>;
class PolyModel : public Model, public ContactGeometryProvider {
public:
  std::vector<Point> points;
  explicit PolyModel(std::vector<Point> p) : points(std::move(p)) {}
  PolylineView contact_geometry() const override { return {points.data(), points.size()}; }
  void pack_state(ModelState&) const override {}
  void unpack_state(const ModelState&) override {}
};
class PlainModel : public Model {
public:
  void pack_state(ModelState&) const override {}
  void unpack_state(const ModelState&) override {}
};
const InsertionGridConfig grid{{-2,-2,-2}, .5, {8,8,8}};
const ModelType plain_type = static_cast<ModelType>(100);
std::unique_ptr<PolyModel> line(real_t y) {
  return std::make_unique<PolyModel>(std::vector<Point>{{-1,y,0},{1,y,0}});
}
TEST(ModelInsertion, ConfigurationAndValidation) {
  ModelRuntime runtime;
  auto* insertion = runtime.attach_component<ModelInsertion>();
  EXPECT_EQ(insertion->check_clearance(*line(0),0).status, InsertionStatus::InvalidLifecycle);
  auto bad = grid; bad.cell_width = 0;
  EXPECT_EQ(insertion->configure(bad), InsertionStatus::InvalidInput);
  bad = grid; bad.dimensions[1] = 0;
  EXPECT_EQ(insertion->configure(bad), InsertionStatus::InvalidInput);
  ASSERT_EQ(insertion->configure(grid), InsertionStatus::Accepted);
  EXPECT_EQ(insertion->configure(grid), InsertionStatus::InvalidLifecycle);
  PolyModel empty({});
  EXPECT_EQ(insertion->check_clearance(empty,0).status, InsertionStatus::InvalidInput);
  EXPECT_EQ(insertion->check_clearance(*line(0),-1).status, InsertionStatus::InvalidInput);
  EXPECT_EQ(insertion->check_clearance(*line(0),std::numeric_limits<real_t>::infinity()).status, InsertionStatus::InvalidInput);
  auto nan = line(0); nan->points[0][0] = std::numeric_limits<real_t>::quiet_NaN();
  EXPECT_EQ(insertion->check_clearance(*nan,0).status, InsertionStatus::InvalidInput);
  EXPECT_EQ(insertion->check_clearance(*line(3),0).status, InsertionStatus::OutOfDomain);
  EXPECT_EQ(insertion->try_adopt(plain_type,std::make_unique<PlainModel>(),0).status, InsertionStatus::UnsupportedGeometry);
  EXPECT_EQ(insertion->try_adopt(ModelType::Invalid,line(0),0).status, InsertionStatus::InvalidInput);
  EXPECT_EQ(runtime.size(),0u);
}
TEST(ModelInsertion, RejectionSnapshotsAndStableIds) {
  ModelRuntime runtime;
  auto* insertion = runtime.attach_component<ModelInsertion>();
  ASSERT_EQ(insertion->configure(grid),InsertionStatus::Accepted);
  auto obstacle = line(0); auto* ptr = obstacle.get();
  auto result = insertion->try_adopt(plain_type,std::move(obstacle),.25);
  ASSERT_TRUE(result.accepted()); EXPECT_EQ(result.model_id,0u);
  EXPECT_EQ(insertion->register_model(0),InsertionStatus::InvalidInput);
  EXPECT_EQ(insertion->register_model(99),InsertionStatus::InvalidInput);
  ptr->points = {{-1,1,0},{1,1,0}}; // Snapshot deliberately stays at y=0.
  result = insertion->try_adopt(plain_type,line(.25),.25);
  EXPECT_EQ(result.status,InsertionStatus::ClearanceViolation);
  EXPECT_EQ(result.conflicting_model,0u);
  EXPECT_EQ(result.candidate_segment,0u); EXPECT_EQ(result.existing_segment,0u);
  EXPECT_EQ(runtime.size(),1u); EXPECT_EQ(insertion->size(),1u);
  result = insertion->try_adopt(plain_type,line(.5),.25);
  ASSERT_TRUE(result.accepted()); EXPECT_EQ(result.model_id,1u);
  runtime.adopt_owned_model(plain_type,line(1)); // Explicitly unregistered.
  EXPECT_TRUE(insertion->check_clearance(*line(1),0).accepted());
  EXPECT_EQ(insertion->register_model(2),InsertionStatus::Accepted);
  EXPECT_FALSE(insertion->check_clearance(*line(1),0).accepted());
  ASSERT_EQ(runtime.seal(),0);
  EXPECT_EQ(insertion->check_clearance(*line(1.5),0).status,InsertionStatus::InvalidLifecycle);
  EXPECT_EQ(insertion->register_model(2),InsertionStatus::InvalidLifecycle);
  EXPECT_EQ(insertion->try_adopt(plain_type,line(1.5),0).status,InsertionStatus::InvalidLifecycle);
}
TEST(ModelInsertion, ExactGeometryAndDomainEdges) {
  ModelRuntime runtime;
  auto* insertion = runtime.attach_component<ModelInsertion>();
  insertion->configure(grid);
  ASSERT_TRUE(insertion->try_adopt(plain_type,line(0),0).accepted());
  PolyModel crossing({{0,-1,0},{0,1,0}}), skew({{0,-1,.25},{0,1,.25}});
  EXPECT_FALSE(insertion->check_clearance(crossing,0).accepted());
  EXPECT_TRUE(insertion->check_clearance(skew,.125).accepted());
  EXPECT_FALSE(insertion->check_clearance(skew,.25).accepted());
  PolyModel endpoint({{1,0,0},{2,0,0}}), repeated({{1,0,0},{1,0,0}});
  EXPECT_FALSE(insertion->check_clearance(endpoint,0).accepted());
  EXPECT_FALSE(insertion->check_clearance(repeated,0).accepted());
  EXPECT_TRUE(insertion->check_clearance(*line(2),.5).accepted());
  EXPECT_FALSE(insertion->check_clearance(*line(2),2).accepted());
}
TEST(ModelInsertion, MatchesExhaustiveSegmentChecks) {
  ModelRuntime runtime;
  auto* insertion = runtime.attach_component<ModelInsertion>(); insertion->configure(grid);
  std::mt19937 generator(917);
  std::uniform_real_distribution<real_t> value(-2,2);
  std::vector<std::vector<Point>> obstacles;
  auto points = [&] { std::vector<Point> p(4); for(auto& x:p) for(auto& v:x) v=value(generator); return p; };
  for(size_t i=0;i<8;++i) {
    obstacles.push_back(points());
    runtime.adopt_owned_model(plain_type,std::make_unique<PolyModel>(obstacles.back()));
    ASSERT_EQ(insertion->register_model(i),InsertionStatus::Accepted);
  }
  NarrowPhaseSegmentDistance narrow;
  for(int trial=0;trial<60;++trial) {
    PolyModel candidate(points()); real_t clearance = (trial%4)*.125;
    bool violation=false;
    for(const auto& obstacle:obstacles)
      for(size_t a=0;a+1<candidate.points.size();++a)
        for(size_t b=0;b+1<obstacle.size();++b)
          violation |= narrow.within_clearance({candidate.points[a],candidate.points[a+1]}, {obstacle[b],obstacle[b+1]},clearance);
    EXPECT_EQ(insertion->check_clearance(candidate,clearance).accepted(),!violation);
  }
}
TEST(ModelInsertion, IBAdoptionConsistency) {
  ModelRuntime runtime, other;
  auto* insertion=runtime.attach_component<ModelInsertion>(); insertion->configure(grid);
  auto* ib=runtime.attach_component<IBCoupling>(); auto* foreign=other.attach_component<IBCoupling>();
  BendingRodImplicitConfig config; config.links=2; config.boundary=Boundary::Free;
  auto rod=[&]{ return std::make_unique<IBBendingRodImplicit>(config); };
  EXPECT_EQ(insertion->try_adopt_ib(*foreign,ModelType::BendingRodImplicit,rod(),0).status,InsertionStatus::InvalidInput);
  EXPECT_EQ(insertion->try_adopt_ib(*ib,plain_type,line(0),0).status,InsertionStatus::InvalidInput);
  EXPECT_EQ(insertion->try_adopt_ib(*ib,ModelType::BendingRodImplicit,rod(),0,{2,0,-1,0}).status,InsertionStatus::InvalidInput);
  EXPECT_EQ(runtime.size(),0u); EXPECT_EQ(insertion->size(),0u);
  auto result=insertion->try_adopt_ib(*ib,ModelType::BendingRodImplicit,rod(),0);
  ASSERT_TRUE(result.accepted()); EXPECT_EQ(result.model_id,0u);
  EXPECT_EQ(insertion->try_adopt_ib(*ib,ModelType::BendingRodImplicit,rod(),0).status,InsertionStatus::ClearanceViolation);
  EXPECT_EQ(runtime.size(),1u); EXPECT_EQ(insertion->size(),1u);
  ib_model_view_array_t views{}; ASSERT_EQ(ib->get_model_views(views),0); EXPECT_EQ(views.count,1u);
}
TEST(ModelInsertion, ProvidersFollowStructuralInitializationAndRestore) {
  BendingRodImplicitConfig config; config.links=2; config.boundary=Boundary::Free;
  BendingRodImplicit dry(config);
  EXPECT_EQ(dry.contact_geometry().points,dry.get_mesh().get_centerline().data());
  IBBendingRodImplicit rod(config);
  ModelState state; rod.pack_state(state);
  auto mesh=rod.get_mesh(); for(auto& p:mesh.get_centerline()) p[1]+=1;
  rod.apply_initial_condition(mesh);
  EXPECT_EQ(rod.contact_geometry().points,rod.get_mesh().get_centerline().data());
  EXPECT_EQ(rod.contact_geometry().points[0][1],mesh.get_centerline()[0][1]);
  rod.unpack_state(state);
  EXPECT_EQ(rod.contact_geometry().points[0][1],dry.contact_geometry().points[0][1]);
  EulerBeam::EulerBeamBCs bcs{}; bcs.end[0]=EulerBeam::left; bcs.end[1]=EulerBeam::right;
  bcs.type[0]=EulerBeam::free_bc; bcs.type[1]=EulerBeam::free_bc;
  IBEulerBeamHuang beam(1,.01,1,5,bcs);
  auto bm=beam.get_mesh(); for(auto& p:bm.get_centerline()) p[1]=.5;
  beam.apply_initial_condition(bm); beam.pack_state(state);
  EXPECT_EQ(beam.contact_geometry().points,beam.get_mesh().get_centerline().data());
  EXPECT_EQ(beam.contact_geometry().count,5u);
  for(auto& p:bm.get_centerline()) p[1]=1;
  beam.apply_initial_condition(bm); beam.unpack_state(state);
  EXPECT_EQ(beam.contact_geometry().points[0][1],.5);
}
TEST(ModelInsertion, ReplicatedAPIWorksSerially) {
  ModelRuntime runtime;
  auto* insertion=runtime.attach_component<ModelInsertion>(); insertion->configure(grid);
  auto factory=[](const ELFF::CounterRandom&, std::uint64_t) -> std::unique_ptr<Model> { return line(0); };
  auto result=insertion->try_adopt_replicated(plain_type,factory,0,0);
  ASSERT_TRUE(result.accepted()); EXPECT_EQ(result.model_id,0u);
  EXPECT_EQ(insertion->try_adopt_replicated(plain_type,factory,0,0).status,
            InsertionStatus::ClearanceViolation);
  EXPECT_EQ(insertion->try_adopt_replicated(plain_type,factory,0,1).status,
            InsertionStatus::InvalidInput);
  EXPECT_EQ(runtime.size(),1u);
}

TEST(ModelInsertion, OwnsOneFixedStreamAndCountsRejectedAttempts) {
  ModelRuntime runtime;
  auto* insertion=runtime.attach_component<ModelInsertion>(); insertion->configure(grid);
  std::vector<std::uint64_t> indices, draws;
  auto factory=[&](const ELFF::CounterRandom& random, std::uint64_t attempt) -> std::unique_ptr<Model> {
    EXPECT_EQ(random.seed(),1u);
    EXPECT_EQ(random.stream_id(),0x3333333333333333ULL);
    indices.push_back(attempt); draws.push_back(random.u64(attempt));
    if(attempt==2) throw std::bad_alloc();
    return line(0);
  };
  EXPECT_TRUE(insertion->try_adopt_replicated(plain_type,factory,0,0).accepted());
  EXPECT_EQ(insertion->try_adopt_replicated(plain_type,factory,0,0).status,InsertionStatus::ClearanceViolation);
  EXPECT_EQ(insertion->try_adopt_replicated(plain_type,factory,0,0).status,InsertionStatus::PreparationFailed);
  EXPECT_EQ(insertion->try_adopt_replicated(plain_type,factory,0,1).status,InsertionStatus::InvalidInput);
  EXPECT_EQ(insertion->attempts(),3u);
  EXPECT_EQ(indices,(std::vector<std::uint64_t>{0,1,2}));
  EXPECT_EQ(draws,(std::vector<std::uint64_t>{0x4f20f7ec5cf695e7ULL,0x6da6be9294ad1f2fULL,0x7aa471b9e3125c14ULL}));
  EXPECT_EQ(runtime.size(),1u);
}

}

TEST(ModelInsertion, ExactDomainAndCollectiveRegistrationSerial) {
  ModelRuntime runtime;
  auto* insertion = runtime.attach_component<ModelInsertion>();
  ASSERT_EQ(insertion->configure(grid, std::array<real_t,3>{2, 1.1, 2}), InsertionStatus::Accepted);
  EXPECT_EQ(insertion->check_clearance(*line(1.2), 0).status, InsertionStatus::OutOfDomain);
  ASSERT_NE(runtime.adopt_owned_model(plain_type, line(0)), nullptr);
  EXPECT_EQ(insertion->register_model_collective(0), InsertionStatus::Accepted);
  EXPECT_EQ(insertion->register_model_collective(0), InsertionStatus::InvalidInput);
  EXPECT_EQ(insertion->register_model_collective(99), InsertionStatus::InvalidInput);
  EXPECT_EQ(insertion->check_clearance(*line(0), 0).status, InsertionStatus::ClearanceViolation);
  EXPECT_EQ(insertion->size(), 1u);
  EXPECT_EQ(runtime.seal(), 0);
  EXPECT_EQ(insertion->register_model_collective(0), InsertionStatus::InvalidLifecycle);
}
