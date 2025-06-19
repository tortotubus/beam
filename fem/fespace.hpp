#include "fe_coll.hpp"
#include "../mesh/mesh.hpp"

namespace beam {
  namespace fem {
class FiniteElementSpace {
public:
  FiniteElementSpace(const Mesh &mesh, const FiniteElementCollection &fec)
      : mesh(mesh), fec(fec) {};

protected:
  const Mesh &mesh;
  const FiniteElementCollection &fec;
  std::vector<std::unique_ptr<FiniteElement>> elements_fe_;
  //int vdim;
  //int ndofs;
  //int nvdofs, nedofs, nbdofs, lnedofs, lnfdofs;
};
  }
} // namespace beam