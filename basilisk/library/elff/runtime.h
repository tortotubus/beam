#pragma once

@include <assert.h>

#include "models/ibm/IBRuntime.h"

static ib_runtime_t elff_runtime = NULL;

static inline ib_runtime_t elff_runtime_get () {
  if (!elff_runtime)
    elff_runtime = ib_runtime_new ();
  assert (elff_runtime);
  return elff_runtime;
}

static inline void elff_runtime_init () {
  (void) elff_runtime_get ();
}

static inline int elff_runtime_register (ib_model_t model) {
  if (!model)
    return -1;
  return ib_runtime_register (elff_runtime_get (), model);
}

static inline int elff_dump (const char * fname) {
  if (!fname)
    return -1;
  if (!elff_runtime)
    return -1;

  return ib_runtime_checkpoint (elff_runtime_get (), fname);
}

static inline int elff_restore (const char * fname) {
  if (!fname)
    return -1;
  if (!elff_runtime)
    return -1;

  return ib_runtime_restore (elff_runtime_get (), fname);
}

static inline void elff_runtime_free () {
  if (!elff_runtime)
    return;
  ib_runtime_delete (elff_runtime);
  elff_runtime = NULL;
}
