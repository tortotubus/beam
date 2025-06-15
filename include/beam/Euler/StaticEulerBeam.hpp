#pragma once

#include <beam/Euler/EulerBeam.hpp>
#include <beam/FiniteElement/Mesh/Mesh.hpp>

namespace beam {
class StaticEulerBeam {
public:
  StaticEulerBeam() : StaticEulerBeam(0., 1., 1., 20) {};
  StaticEulerBeam(double load, double EI, double length, size_t N_elements)
      : load_(load), EI_(EI), length_(length), h_(length / N_elements),
        N_elements_(N_elements), N_nodes_(N_elements + 1),
        mesh(N_nodes_, N_elements_) {
    // Initialize mesh
    for (size_t i = 0; i < N_nodes_; ++i)
      mesh.add_node({i * h_});
    for (size_t i = 1; i < N_nodes_; ++i)
      mesh.add_element({i - 1, i});
  };

  const Mesh<1, 2>& get_mesh() const { return mesh; }

  ~StaticEulerBeam() = default;

private:
  double load_, EI_, area_, length_, h_;
  size_t N_elements_;
  size_t N_nodes_;
  Mesh<1, 2> mesh;
};
} // namespace beam