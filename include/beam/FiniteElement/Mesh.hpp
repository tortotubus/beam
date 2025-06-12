#pragma once

#include <utility>
#include <vector>

template <std::size_t Dim, std::size_t NNodes> class Mesh {
public:
  Mesh();
  ~Mesh();

  void reserve(std::size_t n_nodes, std::size_t n_elements) {
    nodes.reserve(n_nodes);
    elements.reserve(n_elements);
  };

  size_t add_node(std::array<double, Dim> node) {
    nodes.push_back(node);
    return nodes.size() - 1;
  };

  void add_element(std::array<size_t, NNodes> element) {
    elements.push_back(element);
  };

  size_t size_nodes() { return nodes.size(); }
  size_t size_elements() { return elements.size(); }

private:
  std::vector<std::array<double, Dim>> nodes;
  std::vector<std::array<size_t, NNodes>> elements;
};

class LineMesh2D : public Mesh<2, 2> {};
class LineMesh3D : public Mesh<3, 2> {};