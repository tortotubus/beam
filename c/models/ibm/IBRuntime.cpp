#include "c/models/ibm/IBRuntime.h"

#include "models/ibm/IBModel.hpp"
#include "models/ibm/IBRuntime.hpp"

using namespace ELFF::Models;

extern "C" {

ib_runtime_t ib_runtime_new()
{
  return reinterpret_cast<ib_runtime_t>(new IBRuntime());
}

int ib_runtime_register(ib_runtime_t runtime, ib_model_t model)
{
  if (!runtime || !model)
    return -1;

  auto* rt = reinterpret_cast<IBRuntime*>(runtime);
  auto* m = reinterpret_cast<IBModel*>(model);
  return rt->register_model(*m);
}

int ib_runtime_checkpoint(ib_runtime_t runtime, const char* fname)
{
  if (!runtime || !fname)
    return -1;

  auto* rt = reinterpret_cast<IBRuntime*>(runtime);
  return rt->write_checkpoint(fname);
}

int ib_runtime_restore(ib_runtime_t runtime, const char* fname)
{
  if (!runtime || !fname)
    return -1;

  auto* rt = reinterpret_cast<IBRuntime*>(runtime);
  return rt->read_checkpoint(fname);
}

void ib_runtime_delete(ib_runtime_t runtime)
{
  delete reinterpret_cast<IBRuntime*>(runtime);
}

} // extern "C"
