#include "elff/c/models/ibm/IBForceCoupled.h"
#include "elff/models/ibm/IBForceCoupled.hpp"
#include <cmath>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <stdexcept>

using namespace ELFF::Models;
namespace {
IBForceCoupled &model_from(ib_force_coupled_t handle) {
  if (!handle) throw std::invalid_argument("null force-coupled model");
  return *static_cast<IBForceCoupled*>(handle);
}
ib_mesh_t allocate_mesh(size_t count) {
  if (count > INT_MAX) throw std::invalid_argument("node count exceeds C mesh range");
  ib_mesh_t result{};
  result.n = static_cast<int>(count);
  result.position = static_cast<vertex_t*>(calloc(count, sizeof(vertex_t)));
  result.velocity = static_cast<vertex_t*>(calloc(count, sizeof(vertex_t)));
  result.forces = static_cast<vertex_t*>(calloc(count, sizeof(vertex_t)));
  result.measure = static_cast<double*>(calloc(count, sizeof(double)));
  if (count && (!result.position || !result.velocity || !result.forces || !result.measure)) {
    ib_mesh_free(&result);
    throw std::bad_alloc();
  }
  return result;
}
void copy_mesh(const IBMesh &mesh, ib_mesh_t &result) {
  for (int i = 0; i < result.n; ++i) {
    const auto &p = mesh.GetPoints()[i];
    const auto &v = mesh.GetVelocity()[i];
    result.position[i] = {p.x, p.y, p.z};
    result.velocity[i] = {v.x, v.y, v.z};
    result.measure[i] = mesh.GetMeasures()[i];
  }
}
void report(char *error, size_t capacity, const char *message) {
  if (error && capacity) std::snprintf(error, capacity, "%s", message);
}
}

extern "C" {
int ib_force_coupled_get_number_of_nodes(ib_force_coupled_t handle) {
  try {
    const auto count = model_from(handle).GetNumberOfPoints();
    return count <= INT_MAX ? static_cast<int>(count) : -1;
  } catch (...) { return -1; }
}

ib_mesh_t ib_force_coupled_get_current(ib_force_coupled_t handle) {
  ib_mesh_t result{};
  try {
    auto &model = model_from(handle);
    auto &mesh = model.GetCurrent();
    result = allocate_mesh(mesh.GetNumberOfPoints());
    copy_mesh(mesh, result);
    return result;
  } catch (...) {
    ib_mesh_free(&result);
    return {};
  }
}

int ib_force_coupled_get_next_checked(ib_force_coupled_t handle,
    const vertex_t *force, int n, double dt, ib_mesh_t *output,
    char *error, size_t error_capacity) {
  if (output) *output = {};
  report(error, error_capacity, "");
  ib_mesh_t result{};
  try {
    if (!output) throw std::invalid_argument("null output mesh");
    auto &model = model_from(handle);
    if (n < 0 || static_cast<size_t>(n) != model.GetNumberOfPoints() ||
        (n && !force) || !std::isfinite(dt) || dt < 0)
      throw std::invalid_argument("invalid force array or coupling timestep");
    std::vector<IBMesh::IBVertex> loads(n);
    for (int i = 0; i < n; ++i) {
      if (!std::isfinite(force[i].x) || !std::isfinite(force[i].y) || !std::isfinite(force[i].z))
        throw std::invalid_argument("force densities must be finite");
      loads[i] = {static_cast<ELFF::real_t>(force[i].x),
                  static_cast<ELFF::real_t>(force[i].y), static_cast<ELFF::real_t>(force[i].z)};
      if (!std::isfinite(loads[i].x) || !std::isfinite(loads[i].y) || !std::isfinite(loads[i].z))
        throw std::invalid_argument("force density exceeds configured real range");
    }
    const auto step = static_cast<ELFF::real_t>(dt);
    if (!std::isfinite(step) || (dt > 0 && step == 0))
      throw std::invalid_argument("timestep exceeds configured real range");
    result = allocate_mesh(n); // Allocate before advancing the model.
    const auto &mesh = model.GetNext(std::move(loads), step);
    copy_mesh(mesh, result);
    *output = result;
    return 0;
  } catch (const std::exception &exception) {
    report(error, error_capacity, exception.what());
  } catch (...) {
    report(error, error_capacity, "unknown force-coupled stepping failure");
  }
  ib_mesh_free(&result);
  return -1;
}

ib_mesh_t ib_force_coupled_get_next(ib_force_coupled_t handle,
    vertex_t *force, int n, double dt) {
  ib_mesh_t output{};
  ib_force_coupled_get_next_checked(handle, force, n, dt, &output, nullptr, 0);
  return output;
}
} // extern "C"
