// include/mylib/mylib_c.h
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// Opaque handle to mylib::MyClass
typedef void* mylib_handle_t;

/// Create and destroy
mylib_handle_t mylib_create(int init_value);
void            mylib_destroy(mylib_handle_t handle);

/// Call compute()
int mylib_compute(mylib_handle_t handle, int x);

#ifdef __cplusplus
}
#endif
