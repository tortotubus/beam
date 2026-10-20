#include "elff/c/models/ModelRuntime.h"

#include "elff/io/C/xmlCase.h"
#include "elff/io/CXX/xmlCase.hpp"
#include "elff/models/beam/EulerBeamHuangXML.hpp"
#include "elff/models/articulated/BendingRodImplicitBuilder.hpp"
#include "elff/models/beam/EulerBeamHuangBuilder.hpp"
#include "elff/models/ModelRuntime.hpp"
#include "elff/models/ibm/IBCoupling.hpp"
#include "embeddedELFFDtd.hpp"

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

struct ModelRuntimeDefinition
{
  ModelDistribution distribution = ModelDistribution::Rank0;
  std::vector<ModelRuntimeModelDefinition> models;
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

bool parse_runtime_definition(elff_xml_node_t elff,
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
        std::strcmp(name, "BendingRodSegmentsImplicit") != 0) {
      size_t line = 0;
      elff_xml_node_line(model, &line, nullptr);
      set_error(error,
                i,
                static_cast<long>(line),
                0,
                "Unsupported model in schema v1.");
      return false;
    }
    try {
      if (std::strcmp(name, "EulerBeamHuang") == 0)
        definition.models.emplace_back(parse_euler_beam_huang_xml(model));
      else
        definition.models.emplace_back(parse_bending_rod_implicit_xml(model));
    } catch (const IO::CXX::XMLError& exception) {
      return copy_xml_error(error, i, exception.diagnostic(), exception.what());
    }
  }
  return true;
}

bool
build_runtime(ModelRuntime& runtime,
              const ModelRuntimeDefinition& definition,
              model_runtime_error* error)
{
  auto* coupling = runtime.attach_component<IBCoupling>();
  if (!coupling) return false;
  for (size_t i = 0; i < definition.models.size(); ++i) {
    const auto* huang =
      std::get_if<EulerBeamHuangDefinition>(&definition.models[i]);
    const auto* rod = std::get_if<BendingRodImplicitDefinition>(&definition.models[i]);
    const auto type = huang ? ModelType::EulerBeamHuang
                           : ModelType::BendingRodImplicit;
    const long source_line = huang ? huang->source_line : rod->source_line;
    ib_model_info_t info = {1, 0, -1, 0};
    info.node_count = huang ? static_cast<int>(huang->nodes) : rod->config.links + 1;
    info.owner_rank = definition.distribution == ModelDistribution::RoundRobin
      ? static_cast<int>(i % static_cast<size_t>(runtime.ranks())) : 0;
    if (runtime.rank() != info.owner_rank) {
      if (coupling->add_remote_model(type, IB_MODEL_VIEW_FORCE_COUPLED, info) == invalid_model_id)
        return false;
      continue;
    }
    try {
      std::unique_ptr<IBModel> model;
      if (huang) {
        EulerBeamMesh initial_mesh = make_huang_direction_initial_mesh(
          huang->anchor, huang->nodes, huang->length, huang->direction, huang->anchor_at_right);
        model = build_euler_beam_huang(std::move(initial_mesh), huang->length,
          huang->bending_rigidity, huang->linear_density, huang->nodes, huang->boundary_conditions);
      } else {
        model = build_bending_rod_implicit(*rod);
      }
      if (!coupling->adopt_owned_model(type, std::move(model), info)) {
        set_error(error, i, source_line, 0, "Unable to adopt model.");
        return false;
      }
    } catch (const std::exception& exception) {
      set_error(error, i, source_line, 0, exception.what());
      return false;
    } catch (...) {
      set_error(error,
                i,
                source_line,
                0,
                "Unknown error while constructing model.");
      return false;
    }
  }
  return runtime.seal() == 0;
}

int load_runtime_document(ModelRuntime& runtime, elff_xml_document_t document,
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
    char diagnostic[ELFF_XML_ERROR_MESSAGE_CAPACITY] = {};
    auto validation = elff_xml_dtd_validate_node(
      elff, reinterpret_cast<const char*>(elff_embedded_elff_v1_dtd),
      elff_embedded_elff_v1_dtd_size, diagnostic, sizeof(diagnostic));
    if (validation != ELFF_XML_DTD_VALID) {
      size_t line = 0;
      elff_xml_node_line(elff, &line, nullptr);
      set_error(error, std::numeric_limits<size_t>::max(), line, 0, diagnostic);
      return -1;
    }

    ModelRuntimeDefinition definition;
    if (!parse_runtime_definition(elff, definition, error))
      return -1;
    return build_runtime(runtime, definition, error) ? 0 : -1;
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
      if (loaded) rc = load_runtime_document(*runtime, document, error);
    }
#ifdef ELFF_USE_MPI
    if (comm != MPI_COMM_NULL) {
      // A local case/constructor failure must fail the factory on every rank,
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
