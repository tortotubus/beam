#pragma once

#include "../fem/geom.hpp"
#include "../general/array.hpp"
#include "vertex.hpp"

namespace beam {
namespace fem {
class Mesh;

class Element
{
protected:
  /// Element's attribute (specifying material property, etc).
  int attribute;

  /// Element's type from the Finite Element's perspective
  Geometry::Type base_geom;

public:
  /// Constants for the classes derived from Element.
  enum Type
  {
    POINT,
    SEGMENT,
    TRIANGLE
  };

  /// Default element constructor.
  explicit Element(Geometry::Type bg = Geometry::POINT)
  {
    attribute = 1;
    base_geom = bg;
  }

  /// Returns element's type
  virtual Type GetType() const = 0;

  /// Return the Element::Type associated with the given Geometry::Type.
  static Type TypeFromGeometry(const Geometry::Type geom);

  Geometry::Type GetGeometryType() const { return base_geom; }

  /// Return element's attribute.
  inline int GetAttribute() const { return attribute; }

  /// Set element's attribute.
  inline void SetAttribute(const int attr) { attribute = attr; }

  /// Get the indices defining the vertices.
  virtual void GetVertices(Array<int>& v) const = 0;

  /// Set the indices defining the vertices.
  virtual void SetVertices(const Array<int>& v) = 0;

  /// Set the indices defining the vertices.
  virtual void SetVertices(const int* ind) = 0;

  /// @note The returned array should NOT be deleted by the caller.
  virtual int* GetVertices() = 0;

  const int* GetVertices() const
  {
    return const_cast<Element*>(this)->GetVertices();
  }

  virtual int GetNVertices() const = 0;

  virtual int GetNEdges() const = 0;

  virtual const int* GetEdgeVertices(int) const = 0;

  virtual int GetNFaces() const = 0;

  virtual int GetNFaceVertices(int fi) const = 0;

  virtual const int* GetFaceVertices(int fi) const = 0;

  /// Mark the longest edge by assuming/changing the order of the vertices.
  // virtual void MarkEdge(const DSTable &v_to_v, const int *length) {}

  /// Return 1 if the element needs refinement in order to get conforming mesh.
  // virtual int NeedRefinement(HashTable<Hashed2> &v_to_v) const { return 0; }

  /// Set current coarse-fine transformation number.
  virtual void ResetTransform(int tr) {}

  /// Add 'tr' to the current chain of coarse-fine transformations.
  virtual void PushTransform(int tr) {}

  /// Return current coarse-fine transformation.
  virtual unsigned GetTransform() const { return 0; }

  /// @note The returned object should be deleted by the caller.
  virtual Element* Duplicate(Mesh* m) const = 0;

  /// Destroys element.
  virtual ~Element() {}
};
}
} // namespace beam