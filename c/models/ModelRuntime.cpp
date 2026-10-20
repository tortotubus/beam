#include "elff/c/models/ModelRuntime.h"
#include "elff/models/ModelRuntime.hpp"
using ELFF::Models::ModelRuntime;
extern "C" {
model_runtime_t model_runtime_new()
{
  try {
    return reinterpret_cast<model_runtime_t>(new ModelRuntime());
  } catch (...) {
    return nullptr;
  }
}

#ifdef ELFF_USE_MPI
int model_runtime_set_communicator(model_runtime_t runtime, MPI_Comm comm)
{
  if (!runtime)
    return -1;

  try {
    return reinterpret_cast<ModelRuntime*>(runtime)->set_communicator(comm);
  } catch (...) {
    return -1;
  }
}
#endif

int model_runtime_checkpoint(model_runtime_t runtime, const char* fname)
{
  if (!runtime || !fname)
    return -1;

  try {
    return reinterpret_cast<ModelRuntime*>(runtime)->write_checkpoint(fname);
  } catch (...) {
    return -1;
  }
}

int model_runtime_restore(model_runtime_t runtime, const char* fname)
{
  if (!runtime || !fname)
    return -1;

  try {
    return reinterpret_cast<ModelRuntime*>(runtime)->read_checkpoint(fname);
  } catch (...) {
    return -1;
  }
}

int model_runtime_write_polydata(model_runtime_t runtime,
                              const char* fname,
                              int overwrite)
{
  if (!runtime || !fname)
    return -1;

  try {
    return reinterpret_cast<ModelRuntime*>(runtime)->write_polydata(fname,
                                                                 overwrite != 0);
  } catch (...) {
    return -1;
  }
}

void model_runtime_delete(model_runtime_t runtime)
{
  delete reinterpret_cast<ModelRuntime*>(runtime);
}

}
