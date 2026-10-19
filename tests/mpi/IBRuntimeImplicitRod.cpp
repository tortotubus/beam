#include "elff/c/models/ibm/IBForceCoupled.h"
#include "elff/models/ibm/IBRuntime.hpp"
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <vector>

using namespace ELFF::Models;
namespace {
void check(bool success, const char *message) {
  if (!success) {
    std::fprintf(stderr, "Implicit rod MPI regression: %s\n", message);
    MPI_Abort(MPI_COMM_WORLD, 1);
    std::abort();
  }
}
std::string xml(const char *policy) {
  return std::string("<Case formatVersion='1'><ELFF schemaVersion='1' dimensions='3' units='SI'><Models distribution='") + policy + R"('>
<EulerBeamHuang><Geometry length='1' nodes='4'/><Material bendingRigidity='0.01' linearDensity='1'/>
<InitialLine anchorEnd='left'><Anchor x='0' y='0'/><Direction x='1' y='0'/></InitialLine>
<BoundaryConditions left='simple' right='free'/></EulerBeamHuang>
<BendingRodSegmentsImplicit><Geometry length='2' nodes='5'/><Material bendingRigidity='0.1' linearDensity='1.5'/>
<RodInitialLine><Anchor x='2' y='0'/><Direction x='1' y='0'/></RodInitialLine>
<RodBoundary type='free'/><Gravity x='0' y='-1'/></BendingRodSegmentsImplicit>
</Models></ELFF></Case>)";
}
ib_runtime_t create(const std::string &text, MPI_Comm comm) {
  elff_xml_document_t doc = nullptr;
  elff_xml_error xml_error{};
  check(elff_xml_document_load_memory(text.data(), text.size(), &doc, &xml_error) == ELFF_XML_OK, xml_error.message);
  ib_runtime_t runtime = nullptr;
  ib_runtime_error error{};
  check(ib_runtime_create_from_case_document_mpi(doc, comm, &runtime, &error) == 0, error.message);
  elff_xml_document_destroy(doc);
  return runtime;
}
void step(ib_runtime_t runtime, double dt) {
  ib_model_view_array_t views{};
  check(ib_runtime_get_model_views(runtime, &views) == 0, "get views");
  for (size_t i = 0; i < views.count; ++i) {
    if (!views.data[i].handle) continue;
    const int count = ib_force_coupled_get_number_of_nodes(views.data[i].handle);
    std::vector<vertex_t> loads(count, {0.1, 0.2, 0.3});
    ib_mesh_t mesh{};
    char error[512]{};
    check(ib_force_coupled_get_next_checked(views.data[i].handle, loads.data(), count, dt,
      &mesh, error, sizeof(error)) == 0, error);
    check(mesh.n == count, "stepping node count");
    ib_mesh_free(&mesh);
  }
}
void equal(const IBRuntimeState &a, const IBRuntimeState &b) {
  check(a.models.size() == b.models.size(), "model count mismatch");
  for (size_t i = 0; i < a.models.size(); ++i)
    check(a.models[i].type == b.models[i].type && a.models[i].state.ints == b.models[i].state.ints &&
          a.models[i].state.reals == b.models[i].state.reals && a.models[i].state.bytes == b.models[i].state.bytes,
          "checkpoint or continuation mismatch");
}
}
int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank, ranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  MPI_Comm reversed;
  MPI_Comm_split(MPI_COMM_WORLD, 0, ranks - 1 - rank, &reversed);
  int local;
  MPI_Comm_rank(reversed, &local);
  auto runtime = create(xml("round-robin"), reversed);
  MPI_Comm_free(&reversed);
  ib_model_view_array_t views{};
  check(ib_runtime_get_model_views(runtime, &views) == 0 && views.count == 2, "mixed views");
  for (size_t i = 0; i < 2; ++i) {
    ib_model_info_t info{};
    check(ib_model_descriptor_get_info(views.data[i].descriptor, &info) == 0, "descriptor");
    check(info.owner_rank == static_cast<int>(i % ranks) && info.node_count == (i == 0 ? 4 : 5) && info.sparse,
      "mixed descriptor metadata");
    check((views.data[i].handle != nullptr) == (local == info.owner_rank), "reversed communicator ownership");
  }
  step(runtime, 0.01);
  step(runtime, 0.01);
  auto *cpp = static_cast<IBRuntime*>(runtime);
  IBRuntimeState saved, expected, actual;
  check(cpp->export_state(saved) == 0, "export state");
  check(saved.models[1].type == IBRuntimeModelType::BendingRodImplicit, "persistent rod identifier");
  const char *path = "implicit_rod_mpi_checkpoint.bin";
  check(ib_runtime_checkpoint(runtime, path) == 0, "distributed checkpoint");
  MPI_Barrier(MPI_COMM_WORLD);
  step(runtime, 0.005);
  check(cpp->export_state(expected) == 0, "continuation state");
  check(ib_runtime_restore(runtime, path) == 0, "distributed restore");
  check(cpp->export_state(actual) == 0, "restored state");
  equal(saved, actual);
  step(runtime, 0.005);
  check(cpp->export_state(actual) == 0, "restored continuation");
  equal(expected, actual);
  ib_runtime_delete(runtime);

  runtime = create(xml("rank0"), MPI_COMM_WORLD);
  check(ib_runtime_get_model_views(runtime, &views) == 0, "rank0 views");
  for (size_t i = 0; i < 2; ++i)
    check((views.data[i].handle != nullptr) == (rank == 0), "rank0 ownership");
  check(ib_runtime_restore(runtime, path) == 0, "restore with new ownership");
  step(runtime, 0.005);
  if (rank == 0) {
    check(static_cast<IBRuntime*>(runtime)->export_state(actual) == 0, "rank0 continuation export");
    equal(expected, actual);
  }
  ib_runtime_delete(runtime);
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) std::filesystem::remove(path);

  auto bad = xml("round-robin");
  if (rank == 0) bad.replace(bad.find("type='free'"), 11, "type='clamped'");
  elff_xml_document_t doc = nullptr;
  elff_xml_error xml_error{};
  check(elff_xml_document_load_memory(bad.data(), bad.size(), &doc, &xml_error) == ELFF_XML_OK, xml_error.message);
  ib_runtime_error error{};
  runtime = nullptr;
  check(ib_runtime_create_from_case_document_mpi(doc, MPI_COMM_WORLD, &runtime, &error) != 0,
    "rank-local rod failure was not propagated");
  check(!runtime && error.message[0], "failed factory published runtime or lost error");
  elff_xml_document_destroy(doc);
  MPI_Finalize();
}
