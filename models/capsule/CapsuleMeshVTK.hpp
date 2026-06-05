#pragma once

#include "elff/models/capsule/CapsuleMesh.hpp"

#include "elff/io/CXX/vtkPolyData.hpp"

#include <array>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace ELFF {
namespace Models {

class CapsuleMeshVTK {
public:
  static IO::CXX::vtkPolyData to_vtk_polydata(const CapsuleMesh &mesh) {
    IO::CXX::vtkPolyData pd;

    const int nNodes = mesh.numNodes();
    const int nTriangles = mesh.numTriangles();

    if (nNodes > 0 && mesh.state.x.cols() != nNodes)
      throw std::invalid_argument(
          "CapsuleMeshVTK requires one position column per mesh node");

    pd.reserve_points(static_cast<size_t>(nNodes));
    pd.reserve_polygons(static_cast<size_t>(nTriangles));

    for (int i = 0; i < nNodes; ++i) {
      const Vec3 position = mesh.state.x.col(i);
      pd.add_point(static_cast<float>(position.x()),
                   static_cast<float>(position.y()),
                   static_cast<float>(position.z()));
    }

    for (const CapsuleTriangleTopology &triangle : mesh.triangles) {
      pd.add_polygon({ static_cast<int64_t>(triangle.nodes[0]),
                       static_cast<int64_t>(triangle.nodes[1]),
                       static_cast<int64_t>(triangle.nodes[2]) });
    }

    addVectorField(pd, mesh.state.v, "velocity");
    addVectorField(pd, mesh.state.f, "force");
    addVectorField(pd, mesh.state.normal, "normal");

    addScalarField(pd, mesh.state.meanCurv, "mean_curvature");
    addScalarField(pd, mesh.state.refCurv, "reference_curvature");
    addScalarField(pd, mesh.state.gaussCurv, "gaussian_curvature");

    return pd;
  }

private:
  static void addVectorField(IO::CXX::vtkPolyData &pd,
                             const Eigen::Matrix3Xd &field,
                             const char *name) {
    if (field.cols() != static_cast<Eigen::Index>(pd.number_of_points()))
      return;

    const int64_t fieldId = pd.add_pointdata_vector(name, 3);
    for (Eigen::Index i = 0; i < field.cols(); ++i) {
      pd.set_pointdata_vector3(fieldId,
                               static_cast<size_t>(i),
                               { field(0, i), field(1, i), field(2, i) });
    }
  }

  static void addScalarField(IO::CXX::vtkPolyData &pd,
                             const Eigen::VectorXd &field,
                             const char *name) {
    if (field.size() != static_cast<Eigen::Index>(pd.number_of_points()))
      return;

    const int64_t fieldId = pd.add_pointdata_scalar(name);
    std::vector<double> &data = pd.get_pointdata(fieldId);
    for (Eigen::Index i = 0; i < field.size(); ++i)
      data[static_cast<size_t>(i)] = field(i);
  }
};

} // namespace Models
} // namespace ELFF
