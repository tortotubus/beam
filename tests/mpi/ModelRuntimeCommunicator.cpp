#include "elff/c/models/ibm/IBCoupling.h"
#include "elff/c/models/ModelRuntime.h"
#include "elff/c/models/ibm/IBForceCoupled.h"

#include "elff/models/ModelRuntime.hpp"
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

namespace {
bool mpi_active = false;

void check(bool condition, const char* message)
{
  if (condition)
    return;
  std::fprintf(stderr, "ModelRuntime communicator regression: %s\n", message);
  if (mpi_active)
    MPI_Abort(MPI_COMM_WORLD, 1);
  std::exit(1);
}
}

int main(int argc, char** argv)
{
  check(argc == 2, "expected a case-file argument");
  model_runtime_t runtime = reinterpret_cast<model_runtime_t>(1);
  model_runtime_error error = {};

  // Communicator validation must precede even opening the case file.
  check(model_runtime_create_from_case_file_mpi(
          "missing-runtime-case.xml", MPI_COMM_WORLD, &runtime, &error) != 0,
        "factory accepted uninitialized MPI");
  check(runtime == nullptr, "failed factory published a runtime");
  check(std::string(error.message).find("MPI is not initialized") != std::string::npos,
        "factory attempted case construction before communicator setup");

  MPI_Init(&argc, &argv);
  mpi_active = true;
  int world_rank = 0, size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  check(model_runtime_create_from_case_file_mpi(
          argv[1], MPI_COMM_NULL, &runtime, &error) != 0,
        "MPI factory silently accepted a null communicator");
  check(runtime == nullptr, "null-communicator failure published a runtime");

  MPI_Comm reversed = MPI_COMM_NULL;
  MPI_Comm_split(MPI_COMM_WORLD, 0, size - world_rank - 1, &reversed);
  int owner_rank = 0;
  MPI_Comm_rank(reversed, &owner_rank);
  check(model_runtime_create_from_case_file_mpi(argv[1], reversed, &runtime, &error) == 0,
        error.message);
  // The runtime must own a duplicate, independent of the caller's handle.
  MPI_Comm_free(&reversed);

  ib_model_view_array_t views = {};
  check(ib_coupling_get_model_views(ib_coupling_attach(runtime), &views) == 0 && views.count == 1,
        "MPI factory did not build the case model");
  check((views.data[0].handle != nullptr) == (owner_rank == 0),
        "default XML construction allocated a physical model outside rank zero");
  if (owner_rank == 0)
    check(ib_force_coupled_get_number_of_nodes(views.data[0].handle) == 9,
          "case geometry differs in MPI construction");

  const std::string checkpoint = "runtime_comm_rank_" + std::to_string(world_rank) + ".bin";
  std::filesystem::remove(checkpoint);
  check(model_runtime_checkpoint(runtime, checkpoint.c_str()) == 0,
        "checkpoint failed after freeing caller communicator");
  check(std::filesystem::exists(checkpoint) == (owner_rank == 0),
        "checkpoint ownership did not use the supplied communicator's ranks");
  std::filesystem::remove(checkpoint);
  model_runtime_delete(runtime);

  // All entries in a multi-model XML case have the same default owner.
  std::ifstream input(argv[1]);
  std::ostringstream contents;
  contents << input.rdbuf();
  std::string xml = contents.str();
  // Render once in memory and pass it directly to the collective factory.
  std::string parameterized = xml;
  parameterized.insert(parameterized.find('>', parameterized.find("<Case")) + 1,
    "<Parameters><Parameter name='nodes' value='9'/></Parameters>");
  auto nodes = parameterized.find("nodes=\"");
  check(nodes != std::string::npos, "fixture has no node count");
  nodes += 7;
  parameterized.replace(nodes, parameterized.find('"', nodes) - nodes, "${parameters.nodes}");
  elff_xml_document_t source = nullptr, rendered = nullptr;
  elff_xml_error xml_error{};
  check(elff_xml_document_load_memory(parameterized.data(), parameterized.size(),
                                       &source, &xml_error) == ELFF_XML_OK, xml_error.message);
  check(elff_xml_case_render(source, &rendered, &xml_error) == ELFF_XML_OK, xml_error.message);
  elff_xml_document_destroy(source);
  check(model_runtime_create_from_case_document_mpi(rendered, MPI_COMM_WORLD,
                                                 &runtime, &error) == 0, error.message);
  check(ib_coupling_get_model_views(ib_coupling_attach(runtime), &views) == 0 && views.count == 1,
        "document factory lost model entries");
  check((views.data[0].handle != nullptr) == (world_rank == 0),
        "document factory failed sparse ownership");
  if (world_rank == 0)
    check(ib_force_coupled_get_number_of_nodes(views.data[0].handle) == 9,
          "document factory did not consume rendered geometry");
  model_runtime_delete(runtime);
  check(model_runtime_create_from_case_document_mpi(world_rank == 0 ? nullptr : rendered,
                                                 MPI_COMM_WORLD, &runtime, &error) != 0,
        "document factory failed to propagate a local error");
  check(runtime == nullptr, "failed document factory published runtime");
  elff_xml_document_destroy(rendered);

  const auto begin = xml.find("<EulerBeamHuang>");
  const auto end = xml.find("</EulerBeamHuang>") + std::string("</EulerBeamHuang>").size();
  xml.insert(xml.find("</Models>"), xml.substr(begin, end - begin));
  const std::string path_multi = "runtime_multi_rank_" + std::to_string(world_rank) + ".xml";
  { std::ofstream output(path_multi); output << xml; }
  check(model_runtime_create_from_case_file_mpi(path_multi.c_str(), MPI_COMM_WORLD,
                                            &runtime, &error) == 0, error.message);
  check(ib_coupling_get_model_views(ib_coupling_attach(runtime), &views) == 0 && views.count == 2,
        "multi-model directory lost entries");
  for (size_t i = 0; i < views.count; ++i) {
    ib_model_info_t info = {};
    check(ib_model_descriptor_get_info(views.data[i].descriptor, &info) == 0 &&
          info.owner_rank == 0 && info.sparse && info.node_count == 9,
          "incorrect default ownership metadata");
    check((views.data[i].handle != nullptr) == (world_rank == 0),
          "physical model allocated outside rank zero");
  }
  model_runtime_delete(runtime);
  std::filesystem::remove(path_multi);

  // Round-robin uses communicator ranks, including a reordered communicator.
  xml.insert(xml.find("</Models>"), xml.substr(begin, end - begin));
  xml.replace(xml.find("<Models>"), 8, "<Models distribution=\"round-robin\">");
  { std::ofstream output(path_multi); output << xml; }
  MPI_Comm_split(MPI_COMM_WORLD, 0, size - world_rank - 1, &reversed);
  check(model_runtime_create_from_case_file_mpi(path_multi.c_str(), reversed,
                                            &runtime, &error) == 0, error.message);
  MPI_Comm_free(&reversed);
  check(ib_coupling_get_model_views(ib_coupling_attach(runtime), &views) == 0 && views.count == 3,
        "round-robin directory lost entries");
  for (size_t i = 0; i < views.count; ++i) {
    ib_model_info_t info = {};
    check(ib_model_descriptor_get_info(views.data[i].descriptor, &info) == 0 &&
          info.owner_rank == static_cast<int>(i % size), "wrong round-robin owner");
    check((views.data[i].handle != nullptr) == (owner_rank == info.owner_rank),
          "round-robin ignored communicator rank");
    if (views.data[i].handle) {
      std::vector<vertex_t> forces(9);
      ib_force_coupled_get_next(views.data[i].handle, forces.data(), 9, 0.001*(i+1));
    }
  }
  auto* cpp = reinterpret_cast<ELFF::Models::ModelRuntime*>(runtime);
  ELFF::Models::RuntimeState before, after;
  check(cpp->export_state(before) == 0, "distributed state collection failed");
  const char* shared = "runtime_round_robin.bin";
  check(model_runtime_checkpoint(runtime, shared) == 0, "mixed-owner checkpoint failed");
  for (size_t i = 0; i < views.count; ++i)
    if (views.data[i].handle) {
      std::vector<vertex_t> forces(9);
      ib_force_coupled_get_next(views.data[i].handle, forces.data(), 9, 0.002);
    }
  check(model_runtime_restore(runtime, shared) == 0, "mixed-owner restore failed");
  check(cpp->export_state(after) == 0, "restored state collection failed");
  for (size_t i = 0; i < before.models.size(); ++i)
    check(before.models[i].state.ints == after.models[i].state.ints &&
          before.models[i].state.reals == after.models[i].state.reals &&
          before.models[i].state.bytes == after.models[i].state.bytes,
          "checkpoint did not restore a remote owner's complete state");
  model_runtime_delete(runtime);
  xml.replace(xml.find("round-robin"), 11, "rank0");
  { std::ofstream output(path_multi); output << xml; }
  check(model_runtime_create_from_case_file_mpi(path_multi.c_str(), MPI_COMM_WORLD,
                                            &runtime, &error) == 0, error.message);
  check(model_runtime_restore(runtime, shared) == 0,
        "cannot restore distributed checkpoint with rank0 policy");
  if (world_rank == 0) {
    cpp = reinterpret_cast<ELFF::Models::ModelRuntime*>(runtime);
    check(cpp->export_state(after) == 0, "rank0 restored state export failed");
    for (size_t i = 0; i < before.models.size(); ++i)
      check(before.models[i].state.ints == after.models[i].state.ints &&
            before.models[i].state.reals == after.models[i].state.reals,
            "changing ownership policy lost checkpoint state");
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (world_rank == 0) std::filesystem::remove(shared);
  model_runtime_delete(runtime);
  std::filesystem::remove(path_multi);

  // One rank cannot read its case: all ranks must fail without leaking a
  // partially constructed runtime or hanging during communicator destruction.
  runtime = reinterpret_cast<model_runtime_t>(1);
  const char* path = world_rank == 0 ? "missing-runtime-case.xml" : argv[1];
  check(model_runtime_create_from_case_file_mpi(path, MPI_COMM_WORLD, &runtime, &error) != 0,
        "rank-local construction failure was not propagated");
  check(runtime == nullptr && error.message[0], "missing collective failure diagnostic");

  MPI_Finalize();
  mpi_active = false;
  check(model_runtime_create_from_case_file_mpi(argv[1], MPI_COMM_WORLD, &runtime, &error) != 0,
        "factory accepted finalized MPI");
  check(runtime == nullptr, "finalized-MPI failure published a runtime");
  check(std::string(error.message).find("already finalized") != std::string::npos,
        "finalized MPI was not diagnosed before constructing models");
  return 0;
}
