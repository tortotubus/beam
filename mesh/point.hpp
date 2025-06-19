#pragma once

#include "../config/config.hpp"
#include "element.hpp"

namespace beam {
namespace fem {
/// Data type point element
class Point : public Element
{
protected:
   int indices[1];

public:
   typedef Geometry::Constants<Geometry::POINT> geom_t;

   Point() : Element(Geometry::POINT) {}

   /// Constructs point by specifying the indices and the attribute.
   Point(const int *ind, int attr = 1);

   /// Return element's type.
   Type GetType() const override { return Element::POINT; }

   /// Get the indices defining the vertices
   void GetVertices(Array<int> &v) const override;

   /// Set the indices defining the vertices
   void SetVertices(const Array<int> &v) override;

   /// @note The returned array should NOT be deleted by the caller.
   int *GetVertices() override { return indices; }

   /// Set the vertices according to the given input.
   void SetVertices(const int *ind) override;

   int GetNVertices() const override { return 1; }

   int GetNEdges() const override { return (0); }

   const int *GetEdgeVertices(int ei) const override { return NULL; }

   int GetNFaces() const override { return 0; }

   int GetNFaceVertices(int) const override { return 0; }

   const int *GetFaceVertices(int fi) const override { return NULL; }

   Element *Duplicate(Mesh *m) const override
   { return new Point(indices, attribute); }

   virtual ~Point() = default;
};
}
}