#include "../../general/dummy.hpp"
#include "basilisk_interface.h"

using beam::Dummy;

extern "C" {

beam_handle_t beam_dummy_create(int init_value) {
    return reinterpret_cast<beam_handle_t>(new Dummy(init_value));
}

void beam_dummy_destroy(beam_handle_t handle) {
    delete reinterpret_cast<Dummy*>(handle);
}

int beam_dummy_compute(beam_handle_t handle, int x) {
    return reinterpret_cast<Dummy*>(handle)->compute(x);
}

} // extern "C"
