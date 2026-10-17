#include "elff/c/models/ibm/IBRuntime.h"

#include "elff/io/C/xmlCase.h"
#include "elff/models/beam/EulerBeamHuangBuilder.hpp"
#include "elff/models/ibm/IBRuntime.hpp"
#include "embeddedELFFDtd.hpp"

#include <array>
#include <cmath>
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

struct EulerBeamHuangDefinition
{
  long source_line;
  real_t length;
  real_t bending_rigidity;
  real_t linear_density;
  size_t nodes;
  std::array<real_t, 3> anchor;
  std::array<real_t, 3> direction;
  bool anchor_at_right;
  EulerBeam::EulerBeamBCs boundary_conditions;
};

using IBRuntimeModelDefinition = std::variant<EulerBeamHuangDefinition>;

struct IBRuntimeDefinition
{
  std::vector<IBRuntimeModelDefinition> models;
};

void clear_error(ib_runtime_error* error)
{
  if (!error)
    return;
  error->model_index = std::numeric_limits<size_t>::max();
  error->line = 0;
  error->column = 0;
  error->message[0] = '\0';
}

void set_error(ib_runtime_error* error,
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

bool copy_xml_error(ib_runtime_error* error,
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
                  ib_runtime_error* error)
{
  elff_xml_error xml_error = {};
  if (elff_xml_node_unique_child(parent, name, &child, &xml_error) !=
      ELFF_XML_OK)
    return copy_xml_error(error, model_index, xml_error, "Invalid XML child.");
  return true;
}

bool required_double(elff_xml_node_t node,
                     const char* name,
                     double& value,
                     size_t model_index,
                     ib_runtime_error* error)
{
  elff_xml_error xml_error = {};
  if (elff_xml_node_attribute_double(node, name, &value, &xml_error) !=
      ELFF_XML_OK)
    return copy_xml_error(
      error, model_index, xml_error, "Invalid floating-point attribute.");
  return true;
}

bool required_xyz(elff_xml_node_t node,
                  std::array<double, 3>& value,
                  size_t model_index,
                  ib_runtime_error* error)
{
  elff_xml_error xml_error = {};
  if (elff_xml_node_attribute_xyz(node, 0.0, value.data(), &xml_error) !=
      ELFF_XML_OK)
    return copy_xml_error(
      error, model_index, xml_error, "Invalid XYZ attributes.");
  return true;
}

bool required_string(elff_xml_node_t node,
                     const char* name,
                     std::string& value,
                     size_t model_index,
                     ib_runtime_error* error)
{
  elff_xml_error xml_error = {};
  size_t size = 0;
  elff_xml_status status = elff_xml_node_attribute_string(
    node, name, nullptr, 0, &size, &xml_error);
  if (status != ELFF_XML_BUFFER_TOO_SMALL || size == 0)
    return copy_xml_error(error, model_index, xml_error, "Invalid attribute.");
  std::vector<char> buffer(size);
  status = elff_xml_node_attribute_string(
    node, name, buffer.data(), buffer.size(), nullptr, &xml_error);
  if (status != ELFF_XML_OK)
    return copy_xml_error(error, model_index, xml_error, "Invalid attribute.");
  value.assign(buffer.data());
  return true;
}

bool parse_bc_type(const std::string& name, EulerBeam::EulerBeamBCType& type)
{
  if (name == "free")
    type = EulerBeam::free_bc;
  else if (name == "simple")
    type = EulerBeam::simple_bc;
  else if (name == "clamped")
    type = EulerBeam::clamped_bc;
  else
    return false;
  return true;
}

bool parse_huang_definition(elff_xml_node_t model,
                            size_t model_index,
                            EulerBeamHuangDefinition& definition,
                            ib_runtime_error* error)
{
  elff_xml_node_t geometry = nullptr;
  elff_xml_node_t material = nullptr;
  elff_xml_node_t initial_line = nullptr;
  elff_xml_node_t anchor = nullptr;
  elff_xml_node_t direction = nullptr;
  elff_xml_node_t boundary = nullptr;
  if (!unique_child(model, "Geometry", geometry, model_index, error) ||
      !unique_child(model, "Material", material, model_index, error) ||
      !unique_child(model, "InitialLine", initial_line, model_index, error) ||
      !unique_child(initial_line, "Anchor", anchor, model_index, error) ||
      !unique_child(initial_line,
                    "Direction",
                    direction,
                    model_index,
                    error) ||
      !unique_child(model,
                    "BoundaryConditions",
                    boundary,
                    model_index,
                    error))
    return false;

  double length = 0.0, bending = 0.0, density = 0.0;
  int nodes = 0;
  elff_xml_error xml_error = {};
  if (!required_double(
        geometry, "length", length, model_index, error) ||
      elff_xml_node_attribute_int(geometry, "nodes", &nodes, &xml_error) !=
        ELFF_XML_OK) {
    if (error && error->message[0] == '\0')
      copy_xml_error(error, model_index, xml_error, "Invalid node count.");
    return false;
  }
  if (!required_double(material,
                       "bendingRigidity",
                       bending,
                       model_index,
                       error) ||
      !required_double(material,
                       "linearDensity",
                       density,
                       model_index,
                       error))
    return false;

  std::array<double, 3> anchor_value = {};
  std::array<double, 3> direction_value = {};
  if (!required_xyz(anchor, anchor_value, model_index, error) ||
      !required_xyz(direction, direction_value, model_index, error))
    return false;

  std::string anchor_end, left_bc, right_bc;
  if (!required_string(initial_line,
                       "anchorEnd",
                       anchor_end,
                       model_index,
                       error) ||
      !required_string(
        boundary, "left", left_bc, model_index, error) ||
      !required_string(
        boundary, "right", right_bc, model_index, error))
    return false;

  size_t line = 0;
  elff_xml_node_line(model, &line, nullptr);
  const bool finite_scalars =
    std::isfinite(length) && std::isfinite(bending) &&
    std::isfinite(density);
  bool finite_vectors = true;
  for (size_t d = 0; d < 3; ++d) {
    finite_vectors = finite_vectors && std::isfinite(anchor_value[d]) &&
                     std::isfinite(direction_value[d]);
  }
  if (!finite_scalars || !(length > 0.0) || !(bending > 0.0) ||
      !(density > 0.0) || nodes < 2) {
    set_error(error,
              model_index,
              static_cast<long>(line),
              0,
              "EulerBeamHuang length, material values, and node count must "
              "be positive (nodes >= 2).");
    return false;
  }
  if (!finite_vectors) {
    set_error(error,
              model_index,
              static_cast<long>(line),
              0,
              "EulerBeamHuang anchor and direction values must be finite.");
    return false;
  }
  const double direction_norm =
    std::sqrt(direction_value[0] * direction_value[0] +
              direction_value[1] * direction_value[1] +
              direction_value[2] * direction_value[2]);
  if (!std::isfinite(direction_norm) || direction_norm == 0.0) {
    set_error(error,
              model_index,
              static_cast<long>(line),
              0,
              "EulerBeamHuang direction must be finite and nonzero.");
    return false;
  }
  if (anchor_end != "left" && anchor_end != "right") {
    set_error(error,
              model_index,
              static_cast<long>(line),
              0,
              "EulerBeamHuang anchorEnd must be 'left' or 'right'.");
    return false;
  }

  EulerBeam::EulerBeamBCs bcs = {};
  bcs.end[0] = EulerBeam::left;
  bcs.end[1] = EulerBeam::right;
  if (!parse_bc_type(left_bc, bcs.type[0]) ||
      !parse_bc_type(right_bc, bcs.type[1])) {
    set_error(error,
              model_index,
              static_cast<long>(line),
              0,
              "Unsupported EulerBeamHuang boundary condition.");
    return false;
  }

  const bool anchor_at_right = anchor_end == "right";
  const real_t sign = anchor_at_right ? -1.0 : 1.0;
  for (size_t d = 0; d < 3; ++d) {
    const real_t a = static_cast<real_t>(anchor_value[d]);
    const real_t tangent = static_cast<real_t>(direction_value[d] /
                                               direction_norm);
    if (!std::isfinite(a) || !std::isfinite(tangent)) {
      set_error(error,
                model_index,
                static_cast<long>(line),
                0,
                "EulerBeamHuang anchor or direction is outside the "
                "configured real range.");
      return false;
    }
    const real_t left_position =
      anchor_at_right ? a + static_cast<real_t>(length) * tangent : a;
    const real_t right_position =
      anchor_at_right ? a : a + static_cast<real_t>(length) * tangent;
    bcs.vals[0].position[d] = left_position;
    bcs.vals[1].position[d] = right_position;
    bcs.vals[0].slope[d] = sign * tangent;
    bcs.vals[1].slope[d] = sign * tangent;
    definition.anchor[d] = a;
    definition.direction[d] = static_cast<real_t>(direction_value[d]);
  }

  definition.source_line = static_cast<long>(line);
  definition.length = static_cast<real_t>(length);
  definition.bending_rigidity = static_cast<real_t>(bending);
  definition.linear_density = static_cast<real_t>(density);
  if (!std::isfinite(definition.length) ||
      !std::isfinite(definition.bending_rigidity) ||
      !std::isfinite(definition.linear_density)) {
    set_error(error,
              model_index,
              static_cast<long>(line),
              0,
              "EulerBeamHuang values are outside the configured real range.");
    return false;
  }
  definition.nodes = static_cast<size_t>(nodes);
  definition.anchor_at_right = anchor_at_right;
  definition.boundary_conditions = bcs;
  return true;
}

bool parse_runtime_definition(elff_xml_node_t elff,
                              IBRuntimeDefinition& definition,
                              ib_runtime_error* error)
{
  elff_xml_error xml_error = {};
  elff_xml_node_t models = nullptr;
  if (!unique_child(elff,
                    "Models",
                    models,
                    std::numeric_limits<size_t>::max(),
                    error))
    return false;

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
    if (std::strcmp(name, "EulerBeamHuang") != 0) {
      size_t line = 0;
      elff_xml_node_line(model, &line, nullptr);
      set_error(error,
                i,
                static_cast<long>(line),
                0,
                "Only <EulerBeamHuang> is supported in schema v1.");
      return false;
    }
    EulerBeamHuangDefinition huang = {};
    if (!parse_huang_definition(model, i, huang, error))
      return false;
    definition.models.emplace_back(std::move(huang));
  }
  return true;
}

std::unique_ptr<IBRuntime>
build_runtime(const IBRuntimeDefinition& definition,
              ib_runtime_error* error)
{
  auto runtime = std::make_unique<IBRuntime>();
  for (size_t i = 0; i < definition.models.size(); ++i) {
    const auto* huang =
      std::get_if<EulerBeamHuangDefinition>(&definition.models[i]);
    if (!huang) {
      set_error(error, i, 0, 0, "Unsupported runtime model definition.");
      return nullptr;
    }
    try {
      EulerBeamMesh initial_mesh = make_huang_direction_initial_mesh(
        huang->anchor,
        huang->nodes,
        huang->length,
        huang->direction,
        huang->anchor_at_right);
      auto beam = build_euler_beam_huang(std::move(initial_mesh),
                                         huang->length,
                                         huang->bending_rigidity,
                                         huang->linear_density,
                                         huang->nodes,
                                         huang->boundary_conditions);
      if (!runtime->adopt_owned_model(IBRuntimeModelType::EulerBeamHuang,
                                      std::move(beam))) {
        set_error(error,
                  i,
                  huang->source_line,
                  0,
                  "Unable to adopt EulerBeamHuang model.");
        return nullptr;
      }
    } catch (const std::exception& exception) {
      set_error(error, i, huang->source_line, 0, exception.what());
      return nullptr;
    } catch (...) {
      set_error(error,
                i,
                huang->source_line,
                0,
                "Unknown error while constructing EulerBeamHuang model.");
      return nullptr;
    }
  }
  runtime->seal();
  return runtime;
}

} // namespace

extern "C" {

ib_runtime_t ib_runtime_new()
{
  try {
    return reinterpret_cast<ib_runtime_t>(new IBRuntime());
  } catch (...) {
    return nullptr;
  }
}

int ib_runtime_create_from_case_file(const char* path,
                                     ib_runtime_t* out_runtime,
                                     ib_runtime_error* error)
{
  clear_error(error);
  if (out_runtime)
    *out_runtime = nullptr;
  if (!path || !out_runtime) {
    set_error(error,
              std::numeric_limits<size_t>::max(),
              0,
              0,
              "Invalid runtime case-file argument.");
    return -1;
  }

  try {
    elff_xml_document_t document = nullptr;
    elff_xml_node_t elff = nullptr;
    elff_xml_error xml_error = {};
    if (elff_xml_case_load_section_file(
          path,
          "Case",
          "ELFF",
          reinterpret_cast<const char*>(elff_embedded_elff_v1_dtd),
          elff_embedded_elff_v1_dtd_size,
          &document,
          &elff,
          &xml_error) != ELFF_XML_OK) {
      copy_xml_error(error,
                     std::numeric_limits<size_t>::max(),
                     xml_error,
                     "Unable to load and validate the <ELFF> section.");
      return -1;
    }
    std::unique_ptr<elff_xml_document_handle,
                    void (*)(elff_xml_document_t)>
      document_owner(document, elff_xml_document_destroy);

    IBRuntimeDefinition definition;
    if (!parse_runtime_definition(elff, definition, error))
      return -1;
    std::unique_ptr<IBRuntime> runtime = build_runtime(definition, error);
    if (!runtime)
      return -1;
    *out_runtime = reinterpret_cast<ib_runtime_t>(runtime.release());
    return 0;
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

int ib_runtime_get_model_views(ib_runtime_t runtime,
                               ib_model_view_array_t* out_views)
{
  if (!out_views)
    return -1;
  out_views->count = 0;
  out_views->data = nullptr;
  if (!runtime)
    return -1;
  try {
    return reinterpret_cast<IBRuntime*>(runtime)->get_model_views(*out_views);
  } catch (...) {
    return -1;
  }
}

#ifdef ELFF_USE_MPI
int ib_runtime_set_communicator(ib_runtime_t runtime, MPI_Comm comm)
{
  if (!runtime)
    return -1;

  try {
    reinterpret_cast<IBRuntime*>(runtime)->set_communicator(comm);
    return 0;
  } catch (...) {
    return -1;
  }
}
#endif

int ib_runtime_checkpoint(ib_runtime_t runtime, const char* fname)
{
  if (!runtime || !fname)
    return -1;

  try {
    return reinterpret_cast<IBRuntime*>(runtime)->write_checkpoint(fname);
  } catch (...) {
    return -1;
  }
}

int ib_runtime_restore(ib_runtime_t runtime, const char* fname)
{
  if (!runtime || !fname)
    return -1;

  try {
    return reinterpret_cast<IBRuntime*>(runtime)->read_checkpoint(fname);
  } catch (...) {
    return -1;
  }
}

int ib_runtime_write_polydata(ib_runtime_t runtime,
                              const char* fname,
                              int overwrite)
{
  if (!runtime || !fname)
    return -1;

  try {
    return reinterpret_cast<IBRuntime*>(runtime)->write_polydata(fname,
                                                                 overwrite != 0);
  } catch (...) {
    return -1;
  }
}

void ib_runtime_delete(ib_runtime_t runtime)
{
  delete reinterpret_cast<IBRuntime*>(runtime);
}

} // extern "C"
