#include "elff/c/models/ibm/IBRuntime.h"

#include "elff/models/ibm/IBModel.hpp"
#include "elff/models/ibm/IBRuntime.hpp"

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

int ib_runtime_deregister(ib_runtime_t runtime, ib_model_t model)
{
  if (!runtime || !model)
    return -1;

  auto* rt = reinterpret_cast<IBRuntime*>(runtime);
  auto* m = reinterpret_cast<IBModel*>(model);
  return rt->deregister_model(*m);
}

int ib_runtime_set_pid(ib_runtime_t runtime, ib_model_t model, int pid)
{
  if (!runtime || !model)
    return -1;

  auto* rt = reinterpret_cast<IBRuntime*>(runtime);
  auto* m = reinterpret_cast<IBModel*>(model);
  return rt->set_model_pid(*m, pid);
}

#ifdef ELFF_USE_MPI
int ib_runtime_set_communicator(ib_runtime_t runtime, MPI_Comm comm)
{
  if (!runtime)
    return -1;

  auto* rt = reinterpret_cast<IBRuntime*>(runtime);
  rt->set_communicator(comm);
  return 0;
}
#endif

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

int ib_runtime_write_polydata(ib_runtime_t runtime,
                              const char* fname,
                              int overwrite)
{
  if (!runtime || !fname)
    return -1;

  auto* rt = reinterpret_cast<IBRuntime*>(runtime);
  return rt->write_polydata(fname, overwrite != 0);
}

void ib_runtime_delete(ib_runtime_t runtime)
{
  delete reinterpret_cast<IBRuntime*>(runtime);
}

} // extern "C"
