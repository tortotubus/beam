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

    addTriangleFields(pd, mesh);

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

  static void addTriangleFields(IO::CXX::vtkPolyData &pd,
                                const CapsuleMesh &mesh) {
    if (mesh.state.triGeom.size() != static_cast<size_t>(mesh.numTriangles()))
      return;

    const int64_t areaId = pd.add_celldata_scalar("triangle_area");
    const int64_t normalId = pd.add_celldata_vector("triangle_normal", 3);
    const int64_t stretchId = pd.add_celldata_vector("stretch", 2);
    const int64_t tensionId = pd.add_celldata_vector("tension", 2);

    std::vector<double> &area = pd.get_celldata(areaId);
    std::vector<double> &stretch = pd.get_celldata(stretchId);
    std::vector<double> &tension = pd.get_celldata(tensionId);

    for (int tid = 0; tid < mesh.numTriangles(); ++tid) {
      const CapsuleTriangleGeometry &geom =
          mesh.state.triGeom[static_cast<size_t>(tid)];

      area[static_cast<size_t>(tid)] = geom.area;
      pd.set_celldata_vector3(normalId,
                              static_cast<size_t>(tid),
                              { geom.normal.x(), geom.normal.y(),
                                geom.normal.z() });

      const size_t vectorOffset = 2 * static_cast<size_t>(tid);
      stretch[vectorOffset + 0] = geom.stretch.x();
      stretch[vectorOffset + 1] = geom.stretch.y();
      tension[vectorOffset + 0] = geom.tension.x();
      tension[vectorOffset + 1] = geom.tension.y();
    }
  }
};

} // namespace Models
} // namespace ELFF
