#pragma once

#include "elff/models/capsule/CapsuleGeometry.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace ELFF {
namespace Models {

struct CapsuleSphereSpec {
  double radius = 1.0;
  Vec3 center = Vec3::Zero();
  int refinements = 0;
};

struct CapsuleEllipsoidSpec {
  Vec3 radii = Vec3::Ones();
  Vec3 center = Vec3::Zero();
  int refinements = 0;
};

struct CapsuleBiconcaveSpec {
  double radius = 1.0;
  Vec3 center = Vec3::Zero();
  int refinements = 0;
};

class CapsuleMeshBuilder {
public:
  static CapsuleMesh icosahedron(double radius = 1.0,
                                 const Vec3& center = Vec3::Zero()) {
    CapsuleMesh mesh;
    mesh.nodeTopo.resize(12);
    mesh.triangles.resize(20);
    mesh.state.resize(12, 20);

    const double phi = 0.5 * (1.0 + std::sqrt(5.0));
    std::vector<Vec3> vertices = {
      {-1.0,  phi, 0.0}, { 1.0,  phi, 0.0},
      {-1.0, -phi, 0.0}, { 1.0, -phi, 0.0},
      {0.0, -1.0,  phi}, {0.0,  1.0,  phi},
      {0.0, -1.0, -phi}, {0.0,  1.0, -phi},
      { phi, 0.0, -1.0}, { phi, 0.0,  1.0},
      {-phi, 0.0, -1.0}, {-phi, 0.0,  1.0}
    };

    for (int i = 0; i < 12; ++i)
      mesh.state.x.col(i) = center + radius * vertices[i].normalized();

    const std::array<std::array<int,3>,20> faces{{
      {{0, 11, 5}}, {{0, 5, 1}}, {{0, 1, 7}}, {{0, 7, 10}},
      {{0, 10, 11}}, {{1, 5, 9}}, {{5, 11, 4}}, {{11, 10, 2}},
      {{10, 7, 6}}, {{7, 1, 8}}, {{3, 9, 4}}, {{3, 4, 2}},
      {{3, 2, 6}}, {{3, 6, 8}}, {{3, 8, 9}}, {{4, 9, 5}},
      {{2, 4, 11}}, {{6, 2, 10}}, {{8, 6, 7}}, {{9, 8, 1}}
    }};
    for (int i = 0; i < 20; ++i)
      mesh.triangles[i].nodes = faces[i];

    orientOutward(mesh);
    rebuildTopology(mesh);
    updateReferenceGeometry(mesh);
    return mesh;
  }

  static CapsuleMesh sphere(const CapsuleSphereSpec& spec = {}) {
    if (spec.radius <= 0.0)
      throw std::invalid_argument("Capsule sphere radius must be positive");
    if (spec.refinements < 0)
      throw std::invalid_argument("Capsule sphere refinements must be nonnegative");

    CapsuleMesh mesh = icosahedron(spec.radius, spec.center);
    for (int i = 0; i < spec.refinements; ++i) {
      refineTriangles(mesh);
      projectToSphere(mesh, spec.radius, spec.center);
    }

    orientOutward(mesh);
    rebuildTopology(mesh);
    updateReferenceGeometry(mesh);
    return mesh;
  }

  static CapsuleMesh ellipsoid(const CapsuleEllipsoidSpec& spec = {}) {
    if ((spec.radii.array() <= 0.0).any())
      throw std::invalid_argument("Capsule ellipsoid radii must be positive");

    CapsuleSphereSpec sphereSpec;
    sphereSpec.radius = 1.0;
    sphereSpec.center = Vec3::Zero();
    sphereSpec.refinements = spec.refinements;
    CapsuleMesh mesh = sphere(sphereSpec);

    for (int i = 0; i < mesh.numNodes(); ++i)
      mesh.state.x.col(i) =
        spec.center + spec.radii.asDiagonal() * mesh.state.x.col(i);

    orientOutward(mesh);
    rebuildTopology(mesh);
    updateReferenceGeometry(mesh);
    return mesh;
  }

  static CapsuleMesh biconcave(const CapsuleBiconcaveSpec& spec = {}) {
    if (spec.radius <= 0.0)
      throw std::invalid_argument("Capsule biconcave radius must be positive");

    CapsuleSphereSpec sphereSpec;
    sphereSpec.radius = spec.radius;
    sphereSpec.center = Vec3::Zero();
    sphereSpec.refinements = spec.refinements;
    CapsuleMesh mesh = sphere(sphereSpec);

    constexpr double c0 = 0.2072;
    constexpr double c1 = 2.0026;
    constexpr double c2 = -1.1228;
    for (int i = 0; i < mesh.numNodes(); ++i) {
      const double rho = std::min(
        1.0,
        std::sqrt(mesh.state.x(0,i) * mesh.state.x(0,i) +
                  mesh.state.x(2,i) * mesh.state.x(2,i)) / spec.radius);
      const double sign = mesh.state.x(1,i) >= 0.0 ? 1.0 : -1.0;
      mesh.state.x(1,i) =
        sign * 0.5 * spec.radius * std::sqrt(1.0 - rho * rho) *
        (c0 + c1 * rho * rho + c2 * rho * rho * rho * rho);
      mesh.state.x.col(i) += spec.center;
    }

    orientOutward(mesh);
    rebuildTopology(mesh);
    updateReferenceGeometry(mesh);
    return mesh;
  }

  static void rebuildTopology(CapsuleMesh& mesh) {
    mesh.nodeTopo.assign(mesh.numNodes(), CapsuleNodeTopology{});
    mesh.edges.clear();

    std::map<std::pair<int,int>, int> edgeIds;
    for (int tid = 0; tid < mesh.numTriangles(); ++tid) {
      auto& tri = mesh.triangles[tid];
      for (int a = 0; a < 3; ++a) {
        const int n0 = tri.nodes[a];
        const int n1 = tri.nodes[(a + 1) % 3];
        if (n0 < 0 || n1 < 0 || n0 >= mesh.numNodes() || n1 >= mesh.numNodes())
          throw std::runtime_error("Capsule triangle references an invalid node");

        const auto key = orderedEdge(n0, n1);
        auto it = edgeIds.find(key);
        if (it == edgeIds.end()) {
          CapsuleEdge edge;
          edge.nodes = {key.first, key.second};
          edge.triangles = {-1, -1};
          edge.restLength =
            (mesh.state.x.col(edge.nodes[0]) - mesh.state.x.col(edge.nodes[1])).norm();
          const int eid = static_cast<int>(mesh.edges.size());
          mesh.edges.push_back(edge);
          it = edgeIds.emplace(key, eid).first;
        }

        const int eid = it->second;
        tri.edges[a] = eid;
        auto& edge = mesh.edges[eid];
        if (edge.triangles[0] == -1)
          edge.triangles[0] = tid;
        else if (edge.triangles[1] == -1 && edge.triangles[0] != tid)
          edge.triangles[1] = tid;
        else if (edge.triangles[0] != tid && edge.triangles[1] != tid)
          throw std::runtime_error("Capsule edge belongs to more than two triangles");

        addUnique(mesh.nodeTopo[n0].neighbors, n1);
        addUnique(mesh.nodeTopo[n1].neighbors, n0);
        addUnique(mesh.nodeTopo[n0].edges, eid);
        addUnique(mesh.nodeTopo[n1].edges, eid);
      }

      for (const int nid : tri.nodes)
        addUnique(mesh.nodeTopo[nid].triangles, tid);
    }
  }

  static void orientOutward(CapsuleMesh& mesh) {
    if (mesh.numNodes() == 0 || mesh.numTriangles() == 0)
      return;

    Vec3 centroid = Vec3::Zero();
    for (int i = 0; i < mesh.numNodes(); ++i)
      centroid += mesh.state.x.col(i);
    centroid /= static_cast<double>(mesh.numNodes());

    for (auto& tri : mesh.triangles) {
      const Vec3 x0 = mesh.state.x.col(tri.nodes[0]);
      const Vec3 x1 = mesh.state.x.col(tri.nodes[1]);
      const Vec3 x2 = mesh.state.x.col(tri.nodes[2]);
      const Vec3 normal = (x1 - x0).cross(x2 - x0);
      const Vec3 triCentroid = (x0 + x1 + x2) / 3.0;
      if (normal.dot(triCentroid - centroid) < 0.0)
        std::swap(tri.nodes[1], tri.nodes[2]);
    }
  }

  static void validateClosedManifold(const CapsuleMesh& mesh) {
    for (const auto& edge : mesh.edges) {
      if (edge.triangles[0] < 0 || edge.triangles[1] < 0)
        throw std::runtime_error("Capsule mesh is not a closed manifold");
    }
  }

private:
  static std::pair<int,int> orderedEdge(int a, int b) {
    return a < b ? std::make_pair(a, b) : std::make_pair(b, a);
  }

  static void addUnique(std::vector<int>& values, int value) {
    if (std::find(values.begin(), values.end(), value) == values.end())
      values.push_back(value);
  }

  static void refineTriangles(CapsuleMesh& mesh) {
    std::map<std::pair<int,int>, int> midpointIds;
    const auto oldTriangles = mesh.triangles;
    std::vector<Vec3> positions(mesh.numNodes());
    for (int i = 0; i < mesh.numNodes(); ++i)
      positions[i] = mesh.state.x.col(i);

    std::vector<CapsuleTriangleTopology> newTriangles;
    newTriangles.reserve(oldTriangles.size() * 4);

    for (const auto& tri : oldTriangles) {
      const int a = tri.nodes[0];
      const int b = tri.nodes[1];
      const int c = tri.nodes[2];
      const int ab = midpointNode(positions, midpointIds, a, b);
      const int bc = midpointNode(positions, midpointIds, b, c);
      const int ca = midpointNode(positions, midpointIds, c, a);

      newTriangles.push_back({{a, ab, ca}, {}});
      newTriangles.push_back({{b, bc, ab}, {}});
      newTriangles.push_back({{c, ca, bc}, {}});
      newTriangles.push_back({{ab, bc, ca}, {}});
    }

    mesh.nodeTopo.resize(positions.size());
    mesh.triangles = std::move(newTriangles);
    mesh.state.resize(mesh.numNodes(), mesh.numTriangles());
    for (int i = 0; i < static_cast<int>(positions.size()); ++i)
      mesh.state.x.col(i) = positions[i];
    rebuildTopology(mesh);
  }

  static int midpointNode(std::vector<Vec3>& positions,
                          std::map<std::pair<int,int>, int>& midpointIds,
                          int a,
                          int b) {
    const auto key = orderedEdge(a, b);
    auto it = midpointIds.find(key);
    if (it != midpointIds.end())
      return it->second;

    const int id = static_cast<int>(positions.size());
    positions.push_back(0.5 * (positions[a] + positions[b]));
    midpointIds.emplace(key, id);
    return id;
  }

  static void projectToSphere(CapsuleMesh& mesh,
                              double radius,
                              const Vec3& center) {
    for (int i = 0; i < mesh.numNodes(); ++i) {
      const Vec3 radial = mesh.state.x.col(i) - center;
      const double norm = radial.norm();
      if (norm <= 0.0)
        throw std::runtime_error("Cannot project capsule node at sphere center");
      mesh.state.x.col(i) = center + radius * radial / norm;
    }
  }

  static void updateReferenceGeometry(CapsuleMesh& mesh) {
    CapsuleGeometryOps::updateTriangleGeometry(mesh);
    CapsuleGeometryOps::updateNodeNormals(mesh);
    CapsuleGeometryOps::updateVolume(mesh);
    mesh.state.initialVolume = mesh.state.volume;
  }
};

} // namespace Models
} // namespace ELFF
