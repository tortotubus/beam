#include "elff/c/models/ibm/IBCoupling.h"
#include "elff/models/ibm/IBCoupling.hpp"
using namespace ELFF::Models;
extern "C" {
ib_coupling_t ib_coupling_attach(model_runtime_t runtime) {
  if (!runtime) return nullptr;
  try { return static_cast<ModelRuntime*>(runtime)->attach_component<IBCoupling>(); }
  catch (...) { return nullptr; }
}
ib_coupling_t ib_coupling_get(model_runtime_t runtime) {
  return runtime ? static_cast<ModelRuntime*>(runtime)->component<IBCoupling>() : nullptr;
}
int ib_coupling_get_model_views(ib_coupling_t coupling, ib_model_view_array_t* views) {
  if (!views) return -1;
  *views = {};
  if (!coupling) return -1;
  try { return static_cast<IBCoupling*>(coupling)->get_model_views(*views); }
  catch (...) { return -1; }
}
int ib_model_descriptor_get_info(ib_model_descriptor_t descriptor, ib_model_info_t* info) {
  if (!descriptor || !info) return -1;
  *info = *static_cast<const ib_model_info_t*>(descriptor);
  return 0;
}
}
