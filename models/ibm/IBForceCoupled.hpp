#pragma once

#include "config/config.hpp"
#include "general/error.hpp"
#include "models/ibm/IBMesh.hpp"
#include "models/ibm/IBModel.hpp"

#include <vector>

namespace ELFF {
namespace Models {

/**
 * @brief Immersed boundary coupling abstract base class for constitutive models
 * that take forces as inputs and give position and velocity as outputs
 */
class IBForceCoupled :
  public IBModel
{
protected:
  IBMesh mesh, mesh_next;

  void CopyCurrentToNext() { mesh_next = mesh; }

  virtual void ComputeNextPoints(std::vector<IBMesh::IBVertex> force,
                                 real_t dt) = 0;

public:
  IBForceCoupled(size_t NumberOfPoints)
    : mesh(NumberOfPoints)
    , mesh_next(NumberOfPoints) {};

  size_t GetNumberOfPoints()
  {
    // std::cout << "IBIBForceCoupledStructureModel::GetNumberOfPoint()\n";
    return mesh.GetNumberOfPoints();
  }

  virtual IBMesh& GetCurrent() { return mesh; }

  virtual IBMesh& GetNext(IBMesh::IBVertex* force,
                                   size_t n,
                                   real_t dt)
  {
    std::vector<IBMesh::IBVertex> force_wrapper(force, force + n);
    return GetNext(force_wrapper, dt);
  }

  virtual IBMesh& GetNext(std::vector<IBMesh::IBVertex> force,
                                   real_t dt)
  {
    ELFF_ASSERT(force.size() == mesh.GetPoints().size(),
                "Force length must match number of points.\n");

    ComputeNextPoints(force, dt);

    return mesh_next;
  }
};

} // namespace Models
} // namespace ELFF