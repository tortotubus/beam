#include "../general/text.hpp"

#include "mesh.hpp"
#include "point.hpp"
#include "segment.hpp"
#include "triangle.hpp"

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <map>
#include <vector>

namespace beam {
namespace fem {
void
Mesh::Load(const std::string& filename)
{
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open mesh file: " + filename);
  }
  ReadMFEMMesh(file);
}

void
Mesh::ReadMFEMMesh(std::istream& input, int version, int& curved)
{
  std::string ident;
  skip_comment_lines(input, '#');
  input >> ident; // 'dimension'

  assert(ident == "dimension");
  input >> Dim;

  skip_comment_lines(input, '#');
  input >> ident; // 'elements'

  assert(ident == "elements");
  input >> NumOfElements;
  elements.SetSize(NumOfElements);
  for (int j = 0; j < NumOfElements; j++) {
    elements[j] = ReadElement(input);
  }

  if (version == 13) {
    skip_comment_lines(input, '#');
    input >> ident; // 'attrbiute_sets'
    // attribute_sets.attr_sets.Load(input);
    // attribute_sets.attr_sets.SortAll();
    // attribute_sets.attr_sets.UniqueAll();
  }

  skip_comment_lines(input, '#');
  input >> ident; // 'boundary'
  input >> NumOfBdrElements;
  boundary.SetSize(NumOfBdrElements);

  for (int j = 0; j < NumOfBdrElements; j++) {
    boundary[j] = ReadElement(input);
  }

  if (version == 13) {
    skip_comment_lines(input, '#');
    input >> ident;

    // bdr_attribute_sets.attr_sets.Load(input);
    // bdr_attribute_sets.attr_sets.SortAll();
    // bdr_attribute_sets.attr_sets.UniqueAll();
  }

  skip_comment_lines(input, '#');
  input >> ident;

  input >> NumOfVertices;
  vertices.SetSize(NumOfVertices);

  input >> std::ws >> ident;

  if (ident != "nodes") {
    spaceDim = atoi(ident.c_str());
    for (int j = 0; j < NumOfVertices; j++) {
      for (int i = 0; i < spaceDim; i++) {
        input >> vertices[j](i);
      }
    }
  } else {
    // prepare to read the nodes
    input >> std::ws;
    curved = 1;
  }
}

Element*
Mesh::ReadElement(std::istream& input)
{
  int attr;
  Element* el;

  input >> attr;
  el = ReadElementWithoutAttr(input);
  el->SetAttribute(attr);

  return el;
}

Element*
Mesh::ReadElementWithoutAttr(std::istream& input)
{
  int geom, nv, *v;
  Element* el;

  input >> geom;
  el = NewElement(geom);
  nv = el->GetNVertices();
  v = el->GetVertices();
  for (int i = 0; i < nv; i++) {
    input >> v[i];
  }
  return el;
}

Element*
Mesh::NewElement(int geom)
{
  switch (geom) {
    case Geometry::POINT:
      return (new Point);
    case Geometry::SEGMENT:
      return (new Segment);
    case Geometry::TRIANGLE:
      return (new Triangle);
    default:
      throw std::runtime_error("Invalid Geometry::Type, geom = " + geom);
  }
}
}
} // namespace beam