#include "../../mesh/vertex.hpp"
#include "vertex.h"
#include <stdexcept>

using beam::fem::Vertex;

extern "C"
{
  // vertex_t vertex_create() {
  //   return new Vertex();
  // }

  // vertex_t vertex_create_2d(double x, double y) {
  //   return new Vertex(x, y);
  // }

  // vertex_t vertex_create_3d(double x, double y, double z) {
  //   return new Vertex(x, y, z);
  // }

  // void vertex_destroy(vertex_t handle) {
  //   if (handle) {
  //     delete reinterpret_cast<Vertex*>(handle);
  //   }
  // }

  // double vertex_get_component(vertex_t handle, int i) {
  //   if (!handle || i < 0 || i > 2) {
  //     return 0.0; // Or some error handling strategy
  //   }
  //   return (double)reinterpret_cast<Vertex*>(handle)->operator()(i);
  // }

  // void vertex_set_component(vertex_t handle, int i, double value) {
  //   if (handle && i >= 0 && i <= 2) {
  //     reinterpret_cast<Vertex*>(handle)->operator()(i) = value;
  //   }
  // }
} // extern "C"
