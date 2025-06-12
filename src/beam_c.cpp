// src/mylib_c.cpp
#include "beam/MyClass.hpp"
#include "beam/beam_c.h"

using beam::MyClass;

extern "C" {

mylib_handle_t mylib_create(int init_value) {
    return reinterpret_cast<mylib_handle_t>(new MyClass(init_value));
}

void mylib_destroy(mylib_handle_t handle) {
    delete reinterpret_cast<MyClass*>(handle);
}

int mylib_compute(mylib_handle_t handle, int x) {
    return reinterpret_cast<MyClass*>(handle)->compute(x);
}

} // extern "C"
