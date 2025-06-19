#include "point.hpp"

namespace beam {
namespace fem {
Point::Point( const int *ind, int attr ) : Element(Geometry::POINT)
{
   attribute = attr;
   indices[0] = ind[0];
}

void Point::GetVertices(Array<int> &v) const
{
   //v.SetSize(1);
   v[0] = indices[0];
}

void Point::SetVertices(const Array<int> &v)
{
   //MFEM_ASSERT(v.Size() == 1, "!");
   indices[0] = v[0];
}


void Point::SetVertices(const int *ind)
{
   indices[0] = ind[0];
}
}
}