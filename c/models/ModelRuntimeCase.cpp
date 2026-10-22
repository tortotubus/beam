#include "elff/c/models/ModelRuntime.h"

#include "elff/io/C/xmlCase.h"
#include "elff/io/CXX/xmlCase.hpp"
#include "elff/io/case/InitialLineXML.hpp"
#include "elff/io/case/models/EulerBeamHuangXML.hpp"
#include "elff/io/case/models/BendingRodImplicitXML.hpp"
#include "elff/models/articulated/BendingRodImplicitBuilder.hpp"
#include "elff/models/beam/EulerBeamHuangBuilder.hpp"
#include "elff/models/runtime/ModelRuntime.hpp"
#include "elff/models/runtime/ModelInsertion.hpp"
#include "elff/models/ibm/IBCoupling.hpp"
#include "elff/general/error.hpp"
#include "embeddedELFFDtd.hpp"

#include <charconv>
#include <chrono>
#include <cmath>
#include <optional>
#include <cstdio>
#include <cstring>
#include <exception>
#include <limits>
#include <memory>
#include <string>
#include <variant>
#include <vector>

using namespace ELFF::Models;
using namespace ELFF;

namespace {

using ModelRuntimeModelDefinition = std::variant<EulerBeamHuangDefinition, BendingRodImplicitDefinition>;

enum class ModelDistribution { Rank0, RoundRobin };

struct InsertionDefinition {
  size_t count = 0;
  std::uint64_t max_attempts = 0;
  real_t clearance = 0;
  long source_line = 0;
};
struct ModelTemplate {
  ModelRuntimeModelDefinition model;
  ELFF::IO::Case::InitialLineTemplate line;
  std::optional<InsertionDefinition> insertion;
  long source_line = 0;
};
struct ModelRuntimeDefinition
{
  ModelDistribution distribution = ModelDistribution::Rank0;
  std::vector<ModelTemplate> models;
  unsigned dimensions = 3;
  std::string units;
  std::optional<ELFF::IO::Case::PlacementBox> domain;
};

void clear_error(model_runtime_error* error)
{
  if (!error)
    return;
  error->model_index = std::numeric_limits<size_t>::max();
  error->line = 0;
  error->column = 0;
  error->message[0] = '\0';
}

void set_error(model_runtime_error* error,
               size_t model_index,
               long line,
               long column,
               const char* message)
{
  if (!error)
    return;
  error->model_index = model_index;
  error->line = line;
  error->column = column;
  std::snprintf(error->message,
                sizeof(error->message),
                "%s",
                message ? message : "Unknown runtime error.");
}

bool copy_xml_error(model_runtime_error* error,
                    size_t model_index,
                    const elff_xml_error& xml_error,
                    const char* fallback)
{
  set_error(error,
            model_index,
            static_cast<long>(xml_error.line),
            static_cast<long>(xml_error.column),
            xml_error.message[0] ? xml_error.message : fallback);
  return false;
}

bool unique_child(elff_xml_node_t parent,
                  const char* name,
                  elff_xml_node_t& child,
                  size_t model_index,
                  model_runtime_error* error)
{
  elff_xml_error xml_error = {};
  if (elff_xml_node_unique_child(parent, name, &child, &xml_error) !=
      ELFF_XML_OK)
    return copy_xml_error(error, model_index, xml_error, "Invalid XML child.");
  return true;
}

std::uint64_t positive_count(const IO::CXX::XMLNode& node, const char* name) {
  auto text = node.string(name);
  std::uint64_t value = 0;
  const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value);
  if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size() || value == 0)
    node.invalid("Insertion counts must be positive integers within uint64 range.");
  return value;
}

bool parse_runtime_definition(elff_xml_node_t root, elff_xml_node_t elff,
                              ModelRuntimeDefinition& definition,
                              model_runtime_error* error)
{
  elff_xml_error xml_error = {};
  elff_xml_node_t models = nullptr;
  if (!unique_child(elff,
                    "Models",
                    models,
                    std::numeric_limits<size_t>::max(),
                    error))
    return false;

  char policy[32] = {};
  const auto status = elff_xml_node_attribute_string(
    models, "distribution", policy, sizeof(policy), nullptr, &xml_error);
  if (status != ELFF_XML_NOT_FOUND && status != ELFF_XML_OK)
    return copy_xml_error(error, std::numeric_limits<size_t>::max(),
                          xml_error, "Invalid model distribution.");
  if (status == ELFF_XML_OK && std::strcmp(policy, "round-robin") == 0)
    definition.distribution = ModelDistribution::RoundRobin;
  else if (status == ELFF_XML_OK && std::strcmp(policy, "rank0") != 0) {
    set_error(error, std::numeric_limits<size_t>::max(), 0, 0,
              "Model distribution must be rank0 or round-robin.");
    return false;
  }

  size_t count = 0;
  if (elff_xml_node_element_child_count(models, &count, &xml_error) !=
      ELFF_XML_OK)
    return copy_xml_error(error,
                          std::numeric_limits<size_t>::max(),
                          xml_error,
                          "Unable to enumerate models.");

  const IO::CXX::XMLNode section(elff);
  const IO::CXX::XMLNode case_node(root);
  const bool v1 = case_node.integer("formatVersion") == 1;
  definition.dimensions = v1 ? section.integer("dimensions") : case_node.integer("dimension");
  if (definition.dimensions != 2 && definition.dimensions != 3)
    case_node.invalid("Case dimension must be 2 or 3.");
  definition.units = v1 ? section.string("units") : "";
  if (auto domain = IO::Case::optional_xml_child(section, "InsertionDomain"))
    definition.domain = IO::Case::parse_placement_box(*domain, definition.dimensions, true);
  bool has_insertion = false;
  std::uint64_t total_attempts = 0;
  size_t total_models = 0;
  definition.models.clear();
  definition.models.reserve(count);
  for (size_t i = 0; i < count; ++i) {
    elff_xml_node_t model = nullptr;
    const char* name = nullptr;
    if (elff_xml_node_element_child_at(models, i, &model, &xml_error) !=
          ELFF_XML_OK ||
        elff_xml_node_name(model, &name, &xml_error) != ELFF_XML_OK)
      return copy_xml_error(
        error, i, xml_error, "Unable to read an ordered model element.");
    if (std::strcmp(name, "EulerBeamHuang") != 0 &&
        std::strcmp(name, "BendingRodSegmentsImplicit") != 0 &&
        std::strcmp(name, "BendingRodSegments") != 0) {
      size_t line = 0;
      elff_xml_node_line(model, &line, nullptr);
      set_error(error,
                i,
                static_cast<long>(line),
                0,
                "Unsupported model in ELFF schema.");
      return false;
    }
    try {
      ModelTemplate entry;
      IO::CXX::XMLNode node(model);
      entry.source_line = node.line();
      if (auto insertion = IO::Case::optional_xml_child(node, "Insertion")) {
        has_insertion = true;
        const auto count = positive_count(*insertion, "count");
        const auto attempts = positive_count(*insertion, "maxAttempts");
        const real_t clearance = insertion->number("clearance");
        const auto limit = (std::numeric_limits<std::uint64_t>::max() - 4) / 5;
        if (count > std::numeric_limits<size_t>::max() || attempts < count ||
            attempts > limit - total_attempts || !std::isfinite(clearance) || clearance < 0)
          insertion->invalid("Invalid insertion count, attempt budget, or clearance.");
        total_attempts += attempts;
        entry.insertion = InsertionDefinition{static_cast<size_t>(count), attempts, clearance,
                                             static_cast<long>(insertion->line())};
        if (!definition.domain) insertion->invalid("Insertion requires InsertionDomain.");
      }
      const bool huang = std::strcmp(name, "EulerBeamHuang") == 0;
      const auto dimensions = definition.domain ? definition.dimensions : 3;
      if (huang) {
        auto beam = IO::Case::EulerBeamHuangXML::parse(node);
        entry.line = IO::Case::EulerBeamHuangXML::parse_placement(
            node, beam, dimensions, bool(entry.insertion));
        entry.model = std::move(beam);
      } else {
        entry.line = IO::Case::parse_initial_line_template(node, false, dimensions, bool(entry.insertion));
        entry.model = IO::Case::parse_bending_rod_implicit_xml(model, entry.line.fixed);
      }
      if (entry.line.anchors && definition.domain)
        for (size_t d = 0; d < 3; ++d)
          if (entry.line.anchors->min[d] < definition.domain->min[d] ||
              entry.line.anchors->max[d] > definition.domain->max[d])
            IO::Case::initial_line_xml(node).child("RandomAnchor").invalid("Sampling box must be inside InsertionDomain.");
      const auto added = entry.insertion ? entry.insertion->count : 1;
      if (added > std::numeric_limits<size_t>::max() - total_models)
        node.invalid("Expanded model count overflows runtime IDs.");
      total_models += added;
      definition.models.push_back(std::move(entry));
    } catch (const IO::CXX::XMLError& exception) {
      return copy_xml_error(error, i, exception.diagnostic(), exception.what());
    }
  }
  if (definition.domain && !has_insertion)
    section.child("InsertionDomain").invalid("InsertionDomain requires at least one Insertion group.");
  return true;
}

bool agree(ModelRuntime& runtime, bool success, model_runtime_error* error) {
  if (!success && error && !error->message[0])
    set_error(error, std::numeric_limits<size_t>::max(), 0, 0, "Unable to prepare runtime construction.");
#ifdef ELFF_USE_MPI
  if (runtime.ranks() > 1) {
    int all = success;
    MPI_Allreduce(MPI_IN_PLACE, &all, 1, MPI_INT, MPI_MIN, runtime.communicator());
    if (!all && success)
      set_error(error, std::numeric_limits<size_t>::max(), 0, 0,
                "Runtime construction failed on another MPI rank.");
    return all != 0;
  }
#endif
  return success;
}

// Compare explicit fields, never object padding or source-document formatting.
std::vector<char> signature(const ModelRuntimeDefinition& definition) {
  std::vector<char> data;
  auto add = [&](const auto& value) {
    const auto* p = reinterpret_cast<const char*>(&value);
    data.insert(data.end(), p, p + sizeof(value));
  };
  add(definition.distribution); add(definition.dimensions);
  add(definition.units.size());
  for (auto c : definition.units) add(c);
  add(bool(definition.domain));
  if (definition.domain) {
    for (auto x : definition.domain->min) add(x);
    for (auto x : definition.domain->max) add(x);
  }
  add(definition.models.size());
  for (const auto& entry : definition.models) {
    add(entry.model.index()); add(bool(entry.insertion));
    if (entry.insertion) {
      add(entry.insertion->count); add(entry.insertion->max_attempts); add(entry.insertion->clearance);
    }
    add(entry.line.fixed.anchor_at_right); add(entry.line.random_direction);
    for (auto x : entry.line.fixed.anchor) add(x);
    for (auto x : entry.line.fixed.direction) add(x);
    add(bool(entry.line.anchors));
    if (entry.line.anchors) {
      for (auto x : entry.line.anchors->min) add(x);
      for (auto x : entry.line.anchors->max) add(x);
    }
    if (const auto* beam = std::get_if<EulerBeamHuangDefinition>(&entry.model)) {
      add(beam->length); add(beam->nodes); add(beam->bending_rigidity); add(beam->linear_density);
      for (auto x : beam->boundary_conditions.type) add(x);
    } else {
      const auto& rod = std::get<BendingRodImplicitDefinition>(entry.model).config;
      add(rod.length); add(rod.links); add(rod.total_mass); add(rod.ei); add(rod.boundary);
      add(rod.damping_rigidity); add(rod.spectral_radius); add(rod.nonlinear_tolerance); add(rod.max_iterations);
    }
  }
  return data;
}

bool agree_definition(ModelRuntime& runtime, const ModelRuntimeDefinition& definition,
                      model_runtime_error* error) {
  std::vector<char> local;
  bool success = true;
  try { local = signature(definition); } catch (...) { success = false; }
  if (!agree(runtime, success, error)) return false;
#ifdef ELFF_USE_MPI
  if (runtime.ranks() > 1) {
    unsigned long long size = local.size(), root_size = size;
    MPI_Bcast(&root_size, 1, MPI_UNSIGNED_LONG_LONG, 0, runtime.communicator());
    bool same = root_size == size;
    char buffer[4096];
    for (unsigned long long offset = 0; offset < root_size; offset += sizeof(buffer)) {
      int count = static_cast<int>(std::min<unsigned long long>(sizeof(buffer), root_size - offset));
      if (runtime.rank() == 0) std::memcpy(buffer, local.data() + offset, count);
      MPI_Bcast(buffer, count, MPI_BYTE, 0, runtime.communicator());
      if (size == root_size && std::memcmp(buffer, local.data() + offset, count)) same = false;
    }
    if (!same) set_error(error, std::numeric_limits<size_t>::max(), 0, 0,
                         "Normalized XML runtime definitions differ across ranks.");
    return agree(runtime, same, error);
  }
#endif
  return true;
}

std::unique_ptr<IBModel> construct(const ModelTemplate& entry, const ELFF::IO::Case::InitialLineDefinition& line) {
  if (const auto* beam = std::get_if<EulerBeamHuangDefinition>(&entry.model))
    return EulerBeamHuangBuilder::build(*beam, line);
  auto rod = std::get<BendingRodImplicitDefinition>(entry.model);
  for (size_t d = 0; d < 3; ++d) { rod.config.pin[d] = line.anchor[d]; rod.direction[d] = line.direction[d]; }
  return build_bending_rod_implicit(rod);
}

const char* rejection_name(InsertionStatus status) {
  switch (status) {
    case InsertionStatus::ClearanceViolation: return "clearance violation";
    case InsertionStatus::OutOfDomain: return "outside insertion domain";
    case InsertionStatus::PreparationFailed: return "candidate preparation failed";
    case InsertionStatus::CollectiveMismatch: return "collective mismatch";
    default: return "invalid insertion setup or geometry";
  }
}

bool build_runtime(ModelRuntime& runtime, const ModelRuntimeDefinition& definition,
                   model_runtime_error* error) {
  IBCoupling* coupling = nullptr;
  ModelInsertion* insertion = nullptr;
  bool success = true;
  try {
    coupling = runtime.attach_component<IBCoupling>();
    success = coupling != nullptr;
    if (definition.domain) {
      const auto& box = *definition.domain;
      real_t longest = 0;
      for (unsigned d = 0; d < definition.dimensions; ++d)
        longest = std::max(longest, box.max[d] - box.min[d]);
      InsertionGridConfig grid{box.min, longest / 32, {1, 1, 1}};
      if (!std::isfinite(grid.cell_width) || grid.cell_width <= 0) success = false;
      else {
        for (unsigned d = 0; d < definition.dimensions; ++d) {
          grid.dimensions[d] = std::max<size_t>(1, static_cast<size_t>(std::ceil((box.max[d] - box.min[d]) / grid.cell_width)));
          if (grid.origin[d] + grid.cell_width * grid.dimensions[d] < box.max[d]) ++grid.dimensions[d];
        }
        insertion = runtime.attach_component<ModelInsertion>();
        success = success && insertion && insertion->configure(grid, box.max) == InsertionStatus::Accepted;
      }
    }
  } catch (...) { success = false; }
  if (!success) set_error(error, std::numeric_limits<size_t>::max(), definition.domain ? definition.domain->source_line : 0,
                           0, "Unable to configure XML insertion domain or coupling.");
  if (!agree(runtime, success, error)) return false;

  // Fixed obstacles first; without insertion, this is the original model order.
  for (int phase = 0; phase < 2; ++phase) {
    for (size_t i = 0; i < definition.models.size(); ++i) {
      const auto& entry = definition.models[i];
      if (bool(entry.insertion) != (phase == 1)) continue;
      const auto* beam = std::get_if<EulerBeamHuangDefinition>(&entry.model);
      const auto type = beam ? ModelType::EulerBeamHuang : ModelType::BendingRodImplicit;
      const int nodes = beam ? static_cast<int>(beam->nodes) : std::get<BendingRodImplicitDefinition>(entry.model).config.links + 1;
      auto owner = [&] { return definition.distribution == ModelDistribution::RoundRobin
          ? static_cast<int>(runtime.size() % static_cast<size_t>(runtime.ranks())) : 0; };
      if (!entry.insertion) {
        const auto id = runtime.size();
        ib_model_info_t info = {1, owner(), -1, nodes};
        success = true;
        try {
          if (runtime.rank() == info.owner_rank)
            success = coupling->adopt_owned_model(type, construct(entry, entry.line.fixed), info) != nullptr;
          else success = coupling->add_remote_model(type, IB_MODEL_VIEW_FORCE_COUPLED, info) != invalid_model_id;
        } catch (const std::exception& ex) {
          success = false; set_error(error, i, entry.source_line, 0, ex.what());
        } catch (...) { success = false; }
        if (!success && (!error || !error->message[0])) set_error(error, i, entry.source_line, 0, "Unable to construct fixed model.");
        if (!agree(runtime, success, error)) return false;
        if (insertion) {
          auto status = insertion->register_model_collective(id);
          if (status != InsertionStatus::Accepted) {
            set_error(error, i, entry.source_line, 0, rejection_name(status)); return false;
          }
        }
      } else {
        const auto& group = *entry.insertion;
        size_t accepted = 0;
        std::uint64_t attempts = 0;
        std::uint64_t clearance_rejections = 0, domain_rejections = 0;
        const auto started = std::chrono::steady_clock::now();
        InsertionStatus last = InsertionStatus::Accepted;
        ModelInsertion::CandidateFactory factory;
        success = true;
        try {
          factory = [&](const CounterRandom& random, std::uint64_t attempt) -> std::unique_ptr<Model> {
            return construct(entry, entry.line.sample(random, attempt, definition.dimensions));
          };
        } catch (...) { success = false; }
        if (!agree(runtime, success, error)) return false;
        if (runtime.rank() == 0) {
          ELFF_LOG("XML insertion group " << i << " (line " << group.source_line
                   << "): requested=" << group.count << ", clearance=" << group.clearance
                   << ", maxAttempts=" << group.max_attempts
                   << ", existing models=" << runtime.size());
        }
        while (accepted < group.count && attempts < group.max_attempts) {
          const int rank = owner();
          auto result = insertion->try_adopt_replicated(type, factory, group.clearance, rank, coupling, {1, rank, -1, nodes});
          ++attempts; last = result.status;
          if (result.accepted()) ++accepted;
          else if (last == InsertionStatus::ClearanceViolation) ++clearance_rejections;
          else if (last == InsertionStatus::OutOfDomain) ++domain_rejections;
          const bool failed = !result.accepted() &&
              last != InsertionStatus::ClearanceViolation && last != InsertionStatus::OutOfDomain;
          if (runtime.rank() == 0 &&
              (attempts <= 5 || attempts % 100 == 0 || accepted == group.count ||
               attempts == group.max_attempts || failed)) {
            ELFF_LOG("XML insertion group " << i << ": attempt=" << attempts
                     << "/" << group.max_attempts << ", accepted=" << accepted << "/" << group.count
                     << ", clearance rejections=" << clearance_rejections
                     << ", domain rejections=" << domain_rejections
                     << ", last=" << (result.accepted() ? "accepted" : rejection_name(last))
                     << ", elapsed="
                     << std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count()
                     << " s");
          }
          if (failed) break;
        }
        if (accepted != group.count) {
          char message[512];
          std::snprintf(message, sizeof(message), "Insertion group %zu: requested %zu, accepted %zu after %llu attempts; %s.",
                        i, group.count, accepted, static_cast<unsigned long long>(attempts), rejection_name(last));
          set_error(error, i, group.source_line, 0, message); return false;
        }
        if (runtime.rank() == 0) {
          ELFF_LOG("XML insertion group " << i << " completed: " << accepted
                   << " models accepted in " << attempts << " attempts.");
        }
      }
    }
  }
  success = false;
  try { success = runtime.seal() == 0; } catch (...) {}
  return agree(runtime, success, error);
}

int load_runtime_document(elff_xml_document_t document, ModelRuntimeDefinition& definition,
                          model_runtime_error* error)
{
  try {
    elff_xml_node_t root = nullptr, elff = nullptr;
    elff_xml_error xml_error = {};
    auto status = elff_xml_document_root(document, "Case", &root, &xml_error);
    if (status == ELFF_XML_OK)
      status = elff_xml_node_unique_child(root, "ELFF", &elff, &xml_error);
    if (status != ELFF_XML_OK) {
      copy_xml_error(error, std::numeric_limits<size_t>::max(), xml_error,
                     "Unable to read the <ELFF> section.");
      return -1;
    }
    const IO::CXX::XMLNode root_node(root), elff_node(elff);
    const int version = root_node.integer("formatVersion");
    if ((version != 1 && version != 2) || (version == 1 && elff_node.integer("schemaVersion") != version))
      elff_node.invalid("Case formatVersion must be 1 or 2; v1 ELFF schemaVersion must match.");
    char diagnostic[ELFF_XML_ERROR_MESSAGE_CAPACITY] = {};
    auto validation = elff_xml_dtd_validate_node(
      elff, reinterpret_cast<const char*>(version == 2 ? elff_embedded_elff_v2_dtd : elff_embedded_elff_v1_dtd),
      version == 2 ? elff_embedded_elff_v2_dtd_size : elff_embedded_elff_v1_dtd_size,
      diagnostic, sizeof(diagnostic));
    if (validation != ELFF_XML_DTD_VALID) {
      size_t line = 0;
      elff_xml_node_line(elff, &line, nullptr);
      set_error(error, std::numeric_limits<size_t>::max(), line, 0, diagnostic);
      return -1;
    }

    if (!parse_runtime_definition(root, elff, definition, error))
      return -1;
    return 0;
  } catch (const IO::CXX::XMLError& exception) {
    copy_xml_error(error, std::numeric_limits<size_t>::max(), exception.diagnostic(), exception.what());
    return -1;
  } catch (const std::exception& exception) {
    set_error(error,
              error ? error->model_index : std::numeric_limits<size_t>::max(),
              0,
              0,
              exception.what());
    return -1;
  } catch (...) {
    set_error(error,
              error ? error->model_index : std::numeric_limits<size_t>::max(),
              0,
              0,
              "Unknown exception while constructing the runtime.");
    return -1;
  }
}

int create_runtime_case(const char* path, elff_xml_document_t document,
                        model_runtime_t* out_runtime,
                        model_runtime_error* error
#ifdef ELFF_USE_MPI
                        , MPI_Comm comm
#endif
                        )
{
  clear_error(error);
  if (out_runtime)
    *out_runtime = nullptr;
  try {
    auto runtime = std::make_unique<ModelRuntime>();
#ifdef ELFF_USE_MPI
    // Establish the runtime's rank context before parsing or constructing models.
    if (comm != MPI_COMM_NULL)
      runtime->set_communicator(comm);
#endif
    int rc = -1;
    ModelRuntimeDefinition definition;
    std::unique_ptr<elff_xml_document_handle, decltype(&elff_xml_document_destroy)>
      document_owner(nullptr, elff_xml_document_destroy);
    if ((!path && !document) || !out_runtime)
      set_error(error, std::numeric_limits<size_t>::max(), 0, 0,
                "Invalid runtime case-file argument.");
    else {
      bool loaded = true;
      if (path) {
        elff_xml_error xml_error{};
        loaded = elff_xml_case_load_rendered_file(path, runtime->rank() == 0,
                                                 &document, &xml_error) == ELFF_XML_OK;
        document_owner.reset(document);
        if (!loaded)
          copy_xml_error(error, std::numeric_limits<size_t>::max(), xml_error,
                         "Unable to render case XML.");
      }
      if (loaded) rc = load_runtime_document(document, definition, error);
    }
#ifdef ELFF_USE_MPI
    if (comm != MPI_COMM_NULL) {
      // A local parse failure must be agreed before any construction collective,
      // so all participants release the duplicated communicator together.
      int local_success = rc == 0;
      int all_success = 0;
      if (MPI_Allreduce(&local_success, &all_success, 1, MPI_INT, MPI_MIN,
                        comm) != MPI_SUCCESS) {
        set_error(error, std::numeric_limits<size_t>::max(), 0, 0,
                  "Unable to agree on runtime construction status.");
        return -1;
      }
      if (!all_success) {
        if (rc == 0)
          set_error(error, std::numeric_limits<size_t>::max(), 0, 0,
                    "Runtime construction failed on another MPI rank.");
        return -1;
      }
    }
#endif
    if (rc != 0)
      return rc;
    if (!agree_definition(*runtime, definition, error) || !build_runtime(*runtime, definition, error))
      return -1;
    *out_runtime = reinterpret_cast<model_runtime_t>(runtime.release());
    return 0;
  } catch (const std::exception& exception) {
    set_error(error, std::numeric_limits<size_t>::max(), 0, 0, exception.what());
    return -1;
  } catch (...) {
    set_error(error, std::numeric_limits<size_t>::max(), 0, 0,
              "Unable to initialize runtime communicator or storage.");
    return -1;
  }
}

} // namespace

extern "C" {

int model_runtime_create_from_case_file(const char* path,
                                     model_runtime_t* out_runtime,
                                     model_runtime_error* error)
{
  return create_runtime_case(path, nullptr, out_runtime, error
#ifdef ELFF_USE_MPI
                             , MPI_COMM_NULL
#endif
                             );
}

int model_runtime_create_from_case_document(elff_xml_document_t document,
                                         model_runtime_t* out_runtime,
                                         model_runtime_error* error)
{
  return create_runtime_case(nullptr, document, out_runtime, error
#ifdef ELFF_USE_MPI
                             , MPI_COMM_NULL
#endif
                             );
}

#ifdef ELFF_USE_MPI
int model_runtime_create_from_case_file_mpi(const char* path, MPI_Comm comm,
                                         model_runtime_t* out_runtime,
                                         model_runtime_error* error)
{
  if (comm == MPI_COMM_NULL) {
    clear_error(error);
    if (out_runtime)
      *out_runtime = nullptr;
    set_error(error, std::numeric_limits<size_t>::max(), 0, 0,
              "MPI case construction requires a non-null communicator.");
    return -1;
  }
  return create_runtime_case(path, nullptr, out_runtime, error, comm);
}
int model_runtime_create_from_case_document_mpi(elff_xml_document_t document, MPI_Comm comm,
                                             model_runtime_t* out_runtime,
                                             model_runtime_error* error)
{
  if (comm == MPI_COMM_NULL) {
    clear_error(error);
    if (out_runtime) *out_runtime = nullptr;
    set_error(error, std::numeric_limits<size_t>::max(), 0, 0,
              "MPI case construction requires a non-null communicator.");
    return -1;
  }
  return create_runtime_case(nullptr, document, out_runtime, error, comm);
}

#endif

} // extern "C"
