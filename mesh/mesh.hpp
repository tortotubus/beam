#pragma once

#include "../general/array.hpp"
#include "element.hpp"
#include "vertex.hpp"

#include <exception>
#include <fstream>
#include <string>
#include <vector>

namespace beam {
namespace fem {
class Mesh {
public:
  void Load(const std::string &filename);


  void Loader(std::istream &input, int generate_edges, std::string parse_tag);

  /// MFEM Mesh Reader
  void ReadMFEMMesh(std::istream &input);
  Element *ReadElement(std::istream &input);
  Element *ReadElementWithoutAttr(std::istream &input);
  Element *NewElement(int geom);

  // int GetNE() const { return static_cast<int>(elements_.size()); }
  // const Element &GetElement(int e) const { return *elements_[e]; }
  // const std::vector<Vertex> &GetVertices() { return vertices_; }

protected:
  int Dim;
  int spaceDim;

  int NumOfVertices, NumOfElements, NumOfBdrElements;
  int NumOfEdges, NumOfFaces;

  Array<Element *> elements;
  Array<Vertex> vertices;
  Array<Element *> boundary;
  Array<Element *> faces;

  // std::vector<Element *> boundary_;
  // std::vector<Element *> faces_;
};
}
} // namespace beam
