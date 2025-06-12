#include "Mesh.hpp"

template <std::size_t Dim, std::size_t NNodes>
class DoFHandler {
public:
  DoFHandler(const Mesh &m): mesh(m) {
    node_to_dof.resize(mesh.size_nodes());
    for (size_t i = 0; i < mesh.size_nodes(); i++) {
      node_to_dof[i] = i;
    }
    n_dofs = mesh.size_nodes();
  }
private:
  const Mesh &mesh;
  std::vector<size_t> node_to_dof;
  size_t n_dofs;
};