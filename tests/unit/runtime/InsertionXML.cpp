#include "../../InsertionXMLFixture.hpp"
#include "elff/models/runtime/ModelInsertion.hpp"
#include "elff/models/ibm/IBCoupling.hpp"
#include "elff/io/case/InitialLineXML.hpp"
#include "elff/io/case/models/EulerBeamHuangXML.hpp"
#include "elff/io/case/models/BendingRodImplicitXML.hpp"
#include "elff/c/models/ModelRuntime.h"
#include "elff/contact/NarrowPhaseSegmentDistance.hpp"
#include "embeddedELFFDtd.hpp"
#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
using namespace ELFF;
using namespace ELFF::Models;
using namespace insertion_fixture;
#ifdef ELFF_USE_LIBXML2
namespace {
using Runtime = std::unique_ptr<void, decltype(&model_runtime_delete)>;
Runtime create(const std::string& xml, model_runtime_error& error) {
  elff_xml_document_t doc = nullptr;
  elff_xml_error diagnostic{};
  if (elff_xml_document_load_memory(xml.data(), xml.size(), &doc, &diagnostic) != ELFF_XML_OK) {
    ADD_FAILURE() << diagnostic.message; return {nullptr, model_runtime_delete};
  }
  model_runtime_t runtime = nullptr;
  model_runtime_create_from_case_document(doc, &runtime, &error);
  elff_xml_document_destroy(doc);
  return {runtime, model_runtime_delete};
}
RuntimeState state(void* raw) {
  RuntimeState result;
  EXPECT_EQ(static_cast<ModelRuntime*>(raw)->export_state(result), 0);
  return result;
}
void equal(const RuntimeState& a, const RuntimeState& b) {
  ASSERT_EQ(a.models.size(), b.models.size());
  for (size_t i = 0; i < a.models.size(); ++i) {
    EXPECT_EQ(a.models[i].type, b.models[i].type);
    EXPECT_EQ(a.models[i].state.ints, b.models[i].state.ints);
    EXPECT_EQ(a.models[i].state.reals, b.models[i].state.reals);
    EXPECT_EQ(a.models[i].state.bytes, b.models[i].state.bytes);
  }
}
}
TEST(InsertionXML, DeterministicMixedModelsAndFixedFirst) {
  model_runtime_error error{};
  auto a = create(mixed(), error); ASSERT_TRUE(a) << error.message;
  auto b = create(mixed(), error); ASSERT_TRUE(b) << error.message;
  equal(state(a.get()), state(b.get()));
  auto& runtime = *static_cast<ModelRuntime*>(a.get());
  ASSERT_EQ(runtime.size(), 7u);
  EXPECT_EQ(runtime.entry(0)->type, ModelType::EulerBeamHuang);
  EXPECT_EQ(runtime.entry(1)->type, ModelType::BendingRodImplicit);
  EXPECT_EQ(runtime.entry(4)->type, ModelType::EulerBeamHuang);
  EXPECT_EQ(runtime.component<ModelInsertion>()->size(), 7u);
  EXPECT_EQ(runtime.component<ModelInsertion>()->register_model(0), InsertionStatus::InvalidLifecycle);
  Contact::NarrowPhaseSegmentDistance narrow;
  for (size_t i = 0; i < runtime.size(); ++i) {
    auto g = dynamic_cast<Contact::ContactGeometryProvider*>(runtime.model(i))->contact_geometry();
    for (size_t p = 0; p < g.count; ++p) EXPECT_EQ(g.points[p][2], 0);
    for (size_t j = 0; j < i; ++j) {
      auto other = dynamic_cast<Contact::ContactGeometryProvider*>(runtime.model(j))->contact_geometry();
      for (size_t p = 1; p < g.count; ++p) for (size_t q = 1; q < other.count; ++q)
        EXPECT_FALSE(narrow.within_clearance({g.points[p-1],g.points[p]}, {other.points[q-1],other.points[q]}, .1));
    }
  }
}
TEST(InsertionXML, AllPlacementChoicesAndThreeDimensions) {
  const auto box3 = replace(domain, "z='0'/></InsertionDomain>", "z='10'/></InsertionDomain>");
  for (bool beam : {false, true}) for (bool random_position : {false, true})
    for (bool random_orientation : {false, true}) {
      auto text = document(model(beam, random_position ? random_anchor : fixed_anchor,
          random_orientation ? "<RandomDirection/>" : fixed_direction,
          "<Insertion count='1' clearance='0' maxAttempts='20'/>"), box3, "rank0", "3");
      // Both fixed and random proposals begin inside the 3D domain.
      text = replace(text, random_position ? "<Min x='1' y='1' z='0'/>" : "<Anchor x='1' y='1'/>",
                           random_position ? "<Min x='1' y='1' z='0'/>" : "<Anchor x='1' y='1' z='5'/>");
      model_runtime_error error{};
      auto runtime = create(text, error); ASSERT_TRUE(runtime) << text << error.message;
    }
}
TEST(InsertionXML, FixedObstaclesAfterGroupsAndBudgetDiagnostics) {
  auto text = document(model(false, fixed_anchor, fixed_direction, group) + model(true, fixed_anchor, fixed_direction));
  model_runtime_error error{};
  auto runtime = create(text, error);
  EXPECT_FALSE(runtime); EXPECT_EQ(error.model_index, 0u); EXPECT_GT(error.line, 0);
  EXPECT_NE(std::string(error.message).find("accepted 0 after 100 attempts"), std::string::npos);
  EXPECT_NE(std::string(error.message).find("clearance violation"), std::string::npos);
  text = document(model(false, "<Anchor x='9.9' y='1'/>", fixed_direction, group));
  runtime = create(text, error); EXPECT_FALSE(runtime);
  EXPECT_NE(std::string(error.message).find("outside insertion domain"), std::string::npos);
}
TEST(InsertionXML, SemanticAndDtdFailures) {
  const auto valid = document(model(false, random_anchor, "<RandomDirection/>", group));
  const std::vector<std::string> bad = {
    replace(valid, domain, ""),
    document(model(false, fixed_anchor, fixed_direction)),
    replace(valid, group, ""),
    replace(valid, "count='3'", "count='0'"),
    replace(valid, "count='3'", "count='18446744073709551616'"),
    replace(valid, "maxAttempts='100'", "maxAttempts='2'"),
    replace(valid, "maxAttempts='100'", "maxAttempts='18446744073709551615'"),
    replace(valid, "clearance='0.1'", "clearance='-1'"),
    replace(valid, "clearance='0.1'", "clearance='NaN'"),
    replace(valid, "<Max x='9'", "<Max x='11'"),
    replace(valid, "<Max x='9'", "<Max x='0'"),
    replace(valid, "<Max x='10'", "<Max x='0'"),
    replace(valid, "<RandomDirection/>", "<Direction x='0' y='0'/>"),
    replace(valid, "<RandomDirection/>", "<Direction x='1' y='0' z='1'/>"),
    replace(valid, "<RandomDirection/>", "<RandomDirection/><Direction x='1' y='0'/>"),
    replace(valid, "<RandomDirection/>", ""),
    replace(valid, "dimensions='2'", "dimensions='3'"), // zero-depth domain
    replace(valid, "<Min x='0' y='0' z='0'/>", "<Min x='0' y='0'/>"), // valid in 2D: checked below separately
  };
  for (size_t i = 0; i + 1 < bad.size(); ++i) {
    model_runtime_error error{}; auto runtime = create(bad[i], error);
    EXPECT_FALSE(runtime) << bad[i]; EXPECT_NE(error.message[0], '\0'); EXPECT_GT(error.line, 0);
  }
  model_runtime_error error{}; EXPECT_TRUE(create(bad.back(), error)) << error.message;
  const auto outside = document(model(false, fixed_anchor, fixed_direction, group) +
                               model(true, "<Anchor x='20' y='1'/>", fixed_direction));
  EXPECT_FALSE(create(outside, error)); EXPECT_EQ(error.model_index, 1u);
}
TEST(InsertionXML, ExactDomainAndDegenerateSamplingBox) {
  auto box = replace(domain, "y='10'", "y='1.1'");
  model_runtime_error error{};
  auto outside = create(document(model(false, "<Anchor x='1' y='1.2'/>", fixed_direction,
    "<Insertion count='1' clearance='0' maxAttempts='2'/>"), box), error);
  EXPECT_FALSE(outside);
  EXPECT_NE(std::string(error.message).find("outside insertion domain"), std::string::npos);
  auto anchor = "<RandomAnchor><Min x='1' y='1'/><Max x='1' y='1'/></RandomAnchor>";
  auto inside = create(document(model(false, anchor, fixed_direction,
    "<Insertion count='1' clearance='0' maxAttempts='1'/>")), error);
  ASSERT_TRUE(inside) << error.message;
  auto geometry = dynamic_cast<Contact::ContactGeometryProvider*>(static_cast<ModelRuntime*>(inside.get())->model(0))->contact_geometry();
  // Free rods reconstruct endpoints from center-of-mass coordinates.
  EXPECT_NEAR(geometry.points[0][0], 1, 1e-14);
  EXPECT_NEAR(geometry.points[0][1], 1, 1e-14);
  // A 3D domain with a missing z bound is invalid, even if XML defaults are absent.
  auto box3 = replace(domain, "z='0'/></InsertionDomain>", "z='10'/></InsertionDomain>");
  box3 = replace(box3, "<Min x='0' y='0' z='0'/>", "<Min x='0' y='0'/>");
  EXPECT_FALSE(create(document(model(false, fixed_anchor, fixed_direction, group), box3, "rank0", "3"), error));
}
TEST(InsertionXML, SamplingUsesFiveReservedCounters) {
  ELFF::IO::Case::InitialLineTemplate definition;
  definition.anchors = ELFF::IO::Case::PlacementBox{{1,2,3}, {4,5,6}, 0};
  definition.random_direction = true;
  CounterRandom random(1, 0x3333333333333333ULL);
  auto full = definition.sample(random, 7, 3);
  for (unsigned d = 0; d < 3; ++d)
    EXPECT_DOUBLE_EQ(full.anchor[d], definition.anchors->min[d] + 3 * random.uniform01(35 + d));
  EXPECT_NEAR(std::hypot(full.direction[0], full.direction[1], full.direction[2]), 1, 1e-14);
  EXPECT_DOUBLE_EQ(full.direction[2], 2 * random.uniform01(39) - 1);
  definition.anchors.reset();
  EXPECT_EQ(full.direction, definition.sample(random, 7, 3).direction);
}
TEST(InsertionXML, StandaloneParsersRejectRandomAndTemplatesOwnTheirData) {
  for (bool beam : {false, true}) {
    const auto text = model(beam, random_anchor, "<RandomDirection/>", group);
    elff_xml_document_t doc = nullptr; elff_xml_error error{}; elff_xml_node_t node = nullptr;
    ASSERT_EQ(elff_xml_document_load_memory(text.data(), text.size(), &doc, &error), ELFF_XML_OK);
    ASSERT_EQ(elff_xml_document_root(doc, beam ? "EulerBeamHuang" : "BendingRodSegmentsImplicit", &node, &error), ELFF_XML_OK);
    if (beam) EXPECT_THROW({
      const IO::CXX::XMLNode model(node);
      const auto d = IO::Case::EulerBeamHuangXML::parse(model);
      IO::Case::EulerBeamHuangXML::parse_placement(model, d);
    }, IO::CXX::XMLError);
    else EXPECT_THROW(IO::Case::parse_bending_rod_implicit_xml(node), IO::CXX::XMLError);
    auto definition = IO::Case::parse_initial_line_template(IO::CXX::XMLNode(node), beam, 2, true);
    elff_xml_document_destroy(doc);
    auto sample = definition.sample(CounterRandom(1, 0x3333333333333333ULL), 5, 2);
    EXPECT_GE(sample.anchor[0], 1); EXPECT_LE(sample.anchor[0], 9);
    EXPECT_NEAR(std::hypot(sample.direction[0], sample.direction[1]), 1, 1e-14);
  }
}
TEST(InsertionXML, FileRenderedDocumentAndCheckpoint) {
  const auto file = std::filesystem::temp_directory_path() / "elff_insertion_case.xml";
  const auto checkpoint = std::filesystem::temp_directory_path() / "elff_insertion_state.bin";
  auto text = mixed();
  text = replace(text, "<ELFF", "<Parameters><Parameter name='n' value='3'/></Parameters><ELFF");
  text = replace(text, "count='3'", "count='${parameters.n}'");
  { std::ofstream out(file); out << text; }
  model_runtime_error error{}; model_runtime_t raw = nullptr;
  ASSERT_EQ(model_runtime_create_from_case_file(file.c_str(), &raw, &error), 0) << error.message;
  Runtime first(raw, model_runtime_delete);
  elff_xml_document_t source = nullptr, rendered = nullptr; elff_xml_error diagnostic{};
  ASSERT_EQ(elff_xml_document_load_memory(text.data(), text.size(), &source, &diagnostic), ELFF_XML_OK);
  ASSERT_EQ(elff_xml_case_render(source, &rendered, &diagnostic), ELFF_XML_OK);
  elff_xml_document_destroy(source);
  ASSERT_EQ(model_runtime_create_from_case_document(rendered, &raw, &error), 0) << error.message;
  elff_xml_document_destroy(rendered);
  Runtime second(raw, model_runtime_delete);
  equal(state(first.get()), state(second.get()));
  EXPECT_EQ(model_runtime_checkpoint(first.get(), checkpoint.c_str()), 0);
  EXPECT_EQ(model_runtime_restore(second.get(), checkpoint.c_str()), 0);
  equal(state(first.get()), state(second.get()));
  std::filesystem::remove(file); std::filesystem::remove(checkpoint);
  std::filesystem::remove(file.parent_path() / "elff_insertion_case.rendered.xml");
}
#else
TEST(InsertionXML, XmlDisabledReportsUnavailable) {
  model_runtime_t runtime = nullptr; model_runtime_error error{};
  EXPECT_NE(model_runtime_create_from_case_file("unused.xml", &runtime, &error), 0);
  EXPECT_EQ(runtime, nullptr); EXPECT_NE(error.message[0], '\0');
}
#endif
