#pragma once

#include "elff/models/capsule/Capsule.hpp"
#include "elff/models/capsule/CapsuleBending.hpp"
#include "elff/models/capsule/CapsuleMeshVTK.hpp"

#include "elff/models/ibm/IBVelocityCoupled.hpp"

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <vector>

namespace ELFF {
namespace Models {

class IBMCapsule
    : public IBVelocityCoupled
    , public Capsule {
public:
  using Capsule::mesh;

  explicit IBMCapsule(CapsuleMesh mesh,
                      std::unique_ptr<CapsuleMembraneLaw> law =
                          std::make_unique<CapsuleNeoHookeanLaw>())
      : IBVelocityCoupled(static_cast<size_t>(mesh.numNodes())),
        Capsule(std::move(mesh), std::move(law)) {
    initializeStructuralState();
    updateNodalMeasures();
    capsuleToIBMeshCurrent();
    capsuleToIBMeshNext();
  }

  explicit IBMCapsule(Capsule capsule)
      : IBVelocityCoupled(static_cast<size_t>(capsule.mesh().numNodes())),
        Capsule(std::move(capsule)) {
    updateNodalMeasures();
    capsuleToIBMeshCurrent();
    capsuleToIBMeshNext();
  }

  void append_vtk_polydata(
      std::vector<IO::CXX::vtkPolyData> &datasets) const override {
    datasets.push_back(CapsuleMeshVTK::to_vtk_polydata(Capsule::mesh()));
  }

  void pack_state(IBModelState &state) const override {
    static constexpr int64_t stateVersion = 1;

    state.ints.clear();
    state.reals.clear();
    state.bytes.clear();

    const CapsuleMesh &capsuleMesh = Capsule::mesh();
    const int nNodes = capsuleMesh.numNodes();
    const int nEdges = capsuleMesh.numEdges();
    const int nTriangles = capsuleMesh.numTriangles();
    const int membraneLawType = packMembraneLawType(Capsule::law());
    const int bendingLawType =
        Capsule::hasBendingLaw() ? packBendingLawType(Capsule::bendingLaw())
                                 : 0;

    state.ints = {stateVersion,
                  static_cast<int64_t>(nNodes),
                  static_cast<int64_t>(nEdges),
                  static_cast<int64_t>(nTriangles),
                  static_cast<int64_t>(membraneLawType),
                  static_cast<int64_t>(bendingLawType),
                  Capsule::volumeConservationEnabled() ? 1 : 0};

    state.reals.reserve(lawRealCount(membraneLawType) +
                        bendingLawRealCount(bendingLawType) + 5 +
                        static_cast<size_t>(nEdges) +
                        15 * static_cast<size_t>(nNodes) +
                        22 * static_cast<size_t>(nTriangles));

    packMembraneLaw(state, Capsule::law(), membraneLawType);
    packBendingLaw(state, bendingLawType);

    for (int d = 0; d < 3; ++d)
      state.reals.push_back(capsuleMesh.state.centroid[d]);
    state.reals.push_back(capsuleMesh.state.volume);
    state.reals.push_back(capsuleMesh.state.initialVolume);

    for (const CapsuleEdge &edge : capsuleMesh.edges)
      state.reals.push_back(edge.restLength);

    packMatrix(state, capsuleMesh.state.x);
    packMatrix(state, capsuleMesh.state.v);
    packMatrix(state, capsuleMesh.state.f);
    packMatrix(state, capsuleMesh.state.normal);
    packVector(state, capsuleMesh.state.meanCurv);
    packVector(state, capsuleMesh.state.refCurv);
    packVector(state, capsuleMesh.state.gaussCurv);

    for (const CapsuleTriangleReference &ref : capsuleMesh.triRef) {
      for (int c = 0; c < 2; ++c)
        for (int r = 0; r < 2; ++r)
          state.reals.push_back(ref.refShape(r, c));
      for (int c = 0; c < 2; ++c)
        for (int r = 0; r < 3; ++r)
          state.reals.push_back(ref.dNdXi(r, c));
      state.reals.push_back(ref.restArea);
    }

    for (const CapsuleTriangleGeometry &geom : capsuleMesh.state.triGeom) {
      state.reals.push_back(geom.area);
      for (int d = 0; d < 3; ++d)
        state.reals.push_back(geom.normal[d]);
      for (int d = 0; d < 3; ++d)
        state.reals.push_back(geom.centroid[d]);
      for (int d = 0; d < 2; ++d)
        state.reals.push_back(geom.stretch[d]);
      for (int d = 0; d < 2; ++d)
        state.reals.push_back(geom.tension[d]);
    }
  }

  void unpack_state(const IBModelState &state) override {
    static constexpr int64_t stateVersion = 1;

    ELFF_VERIFY(state.ints.size() == 7,
                "IBMCapsule::unpack_state(): invalid integer metadata.");
    ELFF_VERIFY(state.ints[0] == stateVersion,
                "IBMCapsule::unpack_state(): unsupported state version.");

    CapsuleMesh &capsuleMesh = Capsule::mesh();
    const int nNodes = capsuleMesh.numNodes();
    const int nEdges = capsuleMesh.numEdges();
    const int nTriangles = capsuleMesh.numTriangles();
    ELFF_VERIFY(state.ints[1] == nNodes,
                "IBMCapsule::unpack_state(): node count mismatch.");
    ELFF_VERIFY(state.ints[2] == nEdges,
                "IBMCapsule::unpack_state(): edge count mismatch.");
    ELFF_VERIFY(state.ints[3] == nTriangles,
                "IBMCapsule::unpack_state(): triangle count mismatch.");

    const int membraneLawType = static_cast<int>(state.ints[4]);
    const int bendingLawType = static_cast<int>(state.ints[5]);
    const bool volumeConservationEnabled = state.ints[6] != 0;

    const size_t expectedReals = lawRealCount(membraneLawType) +
                                 bendingLawRealCount(bendingLawType) + 5 +
                                 static_cast<size_t>(nEdges) +
                                 15 * static_cast<size_t>(nNodes) +
                                 22 * static_cast<size_t>(nTriangles);
    ELFF_VERIFY(state.reals.size() == expectedReals,
                "IBMCapsule::unpack_state(): invalid real buffer size.");

    size_t k = 0;
    unpackMembraneLaw(state, k, membraneLawType);
    unpackBendingLaw(state, k, bendingLawType);

    for (int d = 0; d < 3; ++d)
      capsuleMesh.state.centroid[d] = state.reals[k++];
    capsuleMesh.state.volume = state.reals[k++];
    capsuleMesh.state.initialVolume = state.reals[k++];

    for (CapsuleEdge &edge : capsuleMesh.edges)
      edge.restLength = state.reals[k++];

    unpackMatrix(state, k, capsuleMesh.state.x, 3, nNodes);
    unpackMatrix(state, k, capsuleMesh.state.v, 3, nNodes);
    unpackMatrix(state, k, capsuleMesh.state.f, 3, nNodes);
    unpackMatrix(state, k, capsuleMesh.state.normal, 3, nNodes);
    unpackVector(state, k, capsuleMesh.state.meanCurv, nNodes);
    unpackVector(state, k, capsuleMesh.state.refCurv, nNodes);
    unpackVector(state, k, capsuleMesh.state.gaussCurv, nNodes);

    ELFF_VERIFY(capsuleMesh.triRef.size() == static_cast<size_t>(nTriangles),
                "IBMCapsule::unpack_state(): reference triangle count mismatch.");
    for (CapsuleTriangleReference &ref : capsuleMesh.triRef) {
      for (int c = 0; c < 2; ++c)
        for (int r = 0; r < 2; ++r)
          ref.refShape(r, c) = state.reals[k++];
      for (int c = 0; c < 2; ++c)
        for (int r = 0; r < 3; ++r)
          ref.dNdXi(r, c) = state.reals[k++];
      ref.restArea = state.reals[k++];
    }

    ELFF_VERIFY(capsuleMesh.state.triGeom.size() ==
                    static_cast<size_t>(nTriangles),
                "IBMCapsule::unpack_state(): triangle geometry count mismatch.");
    for (CapsuleTriangleGeometry &geom : capsuleMesh.state.triGeom) {
      geom.area = state.reals[k++];
      for (int d = 0; d < 3; ++d)
        geom.normal[d] = state.reals[k++];
      for (int d = 0; d < 3; ++d)
        geom.centroid[d] = state.reals[k++];
      for (int d = 0; d < 2; ++d)
        geom.stretch[d] = state.reals[k++];
      for (int d = 0; d < 2; ++d)
        geom.tension[d] = state.reals[k++];
    }

    ELFF_VERIFY(k == state.reals.size(),
                "IBMCapsule::unpack_state(): trailing real data.");

    Capsule::setVolumeConservationEnabled(volumeConservationEnabled);
    updateNodalMeasures();
    capsuleToIBMeshCurrent();
    capsuleToIBMeshNext();
  }

protected:
  void ComputeMidpointForces() override {
    ibMeshToCapsule(IBVelocityCoupled::mesh_next);
    computeElasticForces();
    capsuleForcesToIBMesh(IBVelocityCoupled::mesh_next);
  }

  void ComputeNextPoints(std::vector<IBMesh::IBVertex> &velocity,
                         real_t dt) override {
    IBVelocityCoupled::ComputeNextPoints(velocity, dt);
    setIBMeshVelocity(IBVelocityCoupled::mesh, velocity);
    setCapsuleVelocity(velocity);
    ibMeshToCapsule(IBVelocityCoupled::mesh);
    setCapsuleVelocity(velocity);
    postAdvanceCorrection();
    updateNodalMeasures();
    capsuleToIBMeshCurrent();
  }

private:
  void capsuleToIBMeshCurrent() {
    capsuleToIBMesh(IBVelocityCoupled::mesh);
  }

  void capsuleToIBMeshNext() {
    capsuleToIBMesh(IBVelocityCoupled::mesh_next);
  }

  void capsuleToIBMesh(IBMesh &ibMesh) const {
    auto &points = ibMesh.GetPoints();
    auto &velocity = ibMesh.GetVelocity();
    const CapsuleMesh &capsuleMesh = Capsule::mesh();

    for (int nid = 0; nid < capsuleMesh.numNodes(); ++nid) {
      points[static_cast<size_t>(nid)] = toIBVertex(capsuleMesh.state.x.col(nid));
      velocity[static_cast<size_t>(nid)] =
          toIBVertex(capsuleMesh.state.v.col(nid));
    }
    capsuleForcesToIBMesh(ibMesh);
  }

  void ibMeshToCapsule(IBMesh &ibMesh) {
    CapsuleMesh &capsuleMesh = Capsule::mesh();
    const auto &points = ibMesh.GetPoints();
    const auto &velocity = ibMesh.GetVelocity();

    for (int nid = 0; nid < capsuleMesh.numNodes(); ++nid) {
      capsuleMesh.state.x.col(nid) = toVec3(points[static_cast<size_t>(nid)]);
      capsuleMesh.state.v.col(nid) =
          toVec3(velocity[static_cast<size_t>(nid)]);
    }
  }

  void setCapsuleVelocity(const std::vector<IBMesh::IBVertex> &velocity) {
    CapsuleMesh &capsuleMesh = Capsule::mesh();
    for (int nid = 0; nid < capsuleMesh.numNodes(); ++nid)
      capsuleMesh.state.v.col(nid) =
          toVec3(velocity[static_cast<size_t>(nid)]);
  }

  static void setIBMeshVelocity(
      IBMesh &ibMesh,
      const std::vector<IBMesh::IBVertex> &velocity) {
    auto &ibVelocity = ibMesh.GetVelocity();
    for (size_t nid = 0; nid < ibVelocity.size(); ++nid)
      ibVelocity[nid] = velocity[nid];
  }

  void capsuleForcesToIBMesh(IBMesh &ibMesh) const {
    auto &forces = ibMesh.GetForces();
    const auto &measures = ibMesh.GetMeasures();
    const CapsuleMesh &capsuleMesh = Capsule::mesh();
    for (int nid = 0; nid < capsuleMesh.numNodes(); ++nid) {
      const size_t index = static_cast<size_t>(nid);
      ELFF_ASSERT(measures[index] > 0.0,
                  "Capsule nodal measure must be positive.\n");
      forces[index] = toIBVertex(capsuleMesh.state.f.col(nid) / measures[index]);
    }
  }

  void updateNodalMeasures() {
    CapsuleGeometryOps::updateTriangleGeometry(Capsule::mesh());
    const std::vector<double> nodeArea =
        CapsuleBendingAssembler::computeNodeAreas(Capsule::mesh());
    std::vector<real_t> measures(nodeArea.begin(), nodeArea.end());
    SetNodalMeasures(measures);
  }

  static IBMesh::IBVertex toIBVertex(const Vec3 &value) {
    return { static_cast<real_t>(value.x()),
             static_cast<real_t>(value.y()),
             static_cast<real_t>(value.z()) };
  }

  static Vec3 toVec3(const IBMesh::IBVertex &value) {
    return Vec3(value.x, value.y, value.z);
  }

  static int packMembraneLawType(const CapsuleMembraneLaw &law) {
    if (dynamic_cast<const CapsuleNeoHookeanLaw *>(&law))
      return 1;
    if (dynamic_cast<const CapsuleSkalakLaw *>(&law))
      return 2;
    throw std::runtime_error("Unsupported capsule membrane law checkpoint type");
  }

  static int packBendingLawType(const CapsuleBendingLaw &law) {
    if (dynamic_cast<const CapsuleLinearBendingLaw *>(&law))
      return 1;
    if (dynamic_cast<const CapsuleHelfrichBendingLaw *>(&law))
      return 2;
    throw std::runtime_error("Unsupported capsule bending law checkpoint type");
  }

  static size_t lawRealCount(int lawType) {
    if (lawType == 1)
      return 1;
    if (lawType == 2)
      return 2;
    throw std::runtime_error("Unsupported capsule membrane law checkpoint type");
  }

  static size_t bendingLawRealCount(int lawType) {
    if (lawType == 0)
      return 0;
    if (lawType == 1 || lawType == 2)
      return 1;
    throw std::runtime_error("Unsupported capsule bending law checkpoint type");
  }

  static void packMembraneLaw(IBModelState &state,
                              const CapsuleMembraneLaw &law,
                              int lawType) {
    if (lawType == 1) {
      const auto &neo = dynamic_cast<const CapsuleNeoHookeanLaw &>(law);
      state.reals.push_back(neo.Es);
      return;
    }
    if (lawType == 2) {
      const auto &skalak = dynamic_cast<const CapsuleSkalakLaw &>(law);
      state.reals.push_back(skalak.Es);
      state.reals.push_back(skalak.C);
      return;
    }
    throw std::runtime_error("Unsupported capsule membrane law checkpoint type");
  }

  void packBendingLaw(IBModelState &state, int lawType) const {
    if (lawType == 0) {
      return;
    }
    if (lawType == 1) {
      const auto &linear =
          dynamic_cast<const CapsuleLinearBendingLaw &>(Capsule::bendingLaw());
      state.reals.push_back(linear.Eb);
      return;
    }
    if (lawType == 2) {
      const auto &helfrich =
          dynamic_cast<const CapsuleHelfrichBendingLaw &>(Capsule::bendingLaw());
      state.reals.push_back(helfrich.Eb);
      return;
    }
    throw std::runtime_error("Unsupported capsule bending law checkpoint type");
  }

  void unpackMembraneLaw(const IBModelState &state, size_t &k, int lawType) {
    if (lawType == 1) {
      Capsule::setLaw(std::make_unique<CapsuleNeoHookeanLaw>(state.reals[k++]));
      return;
    }
    if (lawType == 2) {
      const double elasticModulus = state.reals[k++];
      const double areaDilatationModulus = state.reals[k++];
      Capsule::setLaw(std::make_unique<CapsuleSkalakLaw>(
          elasticModulus, areaDilatationModulus));
      return;
    }
    throw std::runtime_error("Unsupported capsule membrane law checkpoint type");
  }

  void unpackBendingLaw(const IBModelState &state, size_t &k, int lawType) {
    if (lawType == 0) {
      Capsule::setBendingLaw(nullptr);
      return;
    }
    if (lawType == 1) {
      Capsule::setBendingLaw(
          std::make_unique<CapsuleLinearBendingLaw>(state.reals[k++]));
      return;
    }
    if (lawType == 2) {
      Capsule::setBendingLaw(
          std::make_unique<CapsuleHelfrichBendingLaw>(state.reals[k++]));
      return;
    }
    throw std::runtime_error("Unsupported capsule bending law checkpoint type");
  }

  template <typename Derived>
  static void packMatrix(IBModelState &state,
                         const Eigen::MatrixBase<Derived> &matrix) {
    for (Eigen::Index c = 0; c < matrix.cols(); ++c)
      for (Eigen::Index r = 0; r < matrix.rows(); ++r)
        state.reals.push_back(matrix(r, c));
  }

  template <typename Derived>
  static void packVector(IBModelState &state,
                         const Eigen::MatrixBase<Derived> &vector) {
    for (Eigen::Index i = 0; i < vector.size(); ++i)
      state.reals.push_back(vector(i));
  }

  template <typename MatrixType>
  static void unpackMatrix(const IBModelState &state,
                           size_t &k,
                           MatrixType &matrix,
                           int rows,
                           int cols) {
    matrix.resize(rows, cols);
    for (int c = 0; c < cols; ++c)
      for (int r = 0; r < rows; ++r)
        matrix(r, c) = state.reals[k++];
  }

  template <typename VectorType>
  static void unpackVector(const IBModelState &state,
                           size_t &k,
                           VectorType &vector,
                           int size) {
    vector.resize(size);
    for (int i = 0; i < size; ++i)
      vector(i) = state.reals[k++];
  }
};

} // namespace Models
} // namespace ELFF
