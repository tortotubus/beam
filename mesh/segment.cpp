#include "segment.hpp"

namespace beam {
namespace fem {
Segment::Segment(const int *ind, int attr) : Element(Geometry::SEGMENT) {
  attribute = attr;
  for (int i = 0; i < 2; i++) {
    indices[i] = ind[i];
  }
}

Segment::Segment(int ind1, int ind2, int attr) : Element(Geometry::SEGMENT) {
  attribute = attr;
  indices[0] = ind1;
  indices[1] = ind2;
}

void Segment::SetVertices(const int *ind) {
  indices[0] = ind[0];
  indices[1] = ind[1];
}

void Segment::GetVertices(Array<int> &v) const {
  // v.SetSize(2);
  //std::copy(indices, indices + 2, v.begin());
}

void Segment::SetVertices(const Array<int> &v) {
  //MFEM_ASSERT(v.Size() == 2, "!");
  //std::copy(v.begin(), v.end(), indices);
}
}
} // namespace beam
