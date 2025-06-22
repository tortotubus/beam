#pragma once

#ifdef __cplusplus
extern "C" { 
#endif 

typedef void* beam_handle_t;

beam_handle_t beam_dummy_create(int init_value);
void beam_dummy_destroy(beam_handle_t handle);
int beam_dummy_compute(beam_handle_t handle, int x);

#ifdef __cplusplus
}
#endif