#pragma once

#ifdef __cplusplus
extern "C"
{
#endif

typedef void* ib_runtime_t;
typedef void* ib_model_t;

ib_runtime_t ib_runtime_new();
int ib_runtime_register(ib_runtime_t runtime, ib_model_t model);
int ib_runtime_checkpoint(ib_runtime_t runtime, const char *fname);
int ib_runtime_restore(ib_runtime_t runtime, const char *fname);
void ib_runtime_delete(ib_runtime_t runtime);

#ifdef __cplusplus
}
#endif