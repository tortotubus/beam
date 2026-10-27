#ifndef ELFF_FLUID_DLM_INTERNAL_H_
#define ELFF_FLUID_DLM_INTERNAL_H_

#include <elff/fluid/DLM.H>

#include <AMReX_BLassert.H>
#include <AMReX_BoxIterator.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <vector>

namespace ELFF::Fluid::DLM::Detail {

inline constexpr amrex::Real pi = amrex::Real(3.141592653589793238462643383279502884L);
inline constexpr amrex::Real tiny = amrex::Real(1.0e-30);

struct RigidMotion
{
    amrex::Real mass = 0.0;
    std::array<amrex::Real, 3> velocity{0.0, 0.0, 0.0};
    std::array<amrex::Real, 3> angular_velocity{0.0, 0.0, 0.0};
};

struct GeneralizedForce
{
    std::array<amrex::Real, 3> force{0.0, 0.0, 0.0};
    std::array<amrex::Real, 3> torque{0.0, 0.0, 0.0};
};

using FaceOwnerMasks = amrex::Array<std::unique_ptr<amrex::iMultiFab>, AMREX_SPACEDIM>;

enum class MarkerField
{
    Lambda,
    Residual,
    Search,
    Response
};

struct PackedFaceRef
{
    int local_fab = 0;
    int i = 0;
    int j = 0;
    int k = 0;
    amrex::Real weight = 0.0;
};

struct PackedConstraintContribution
{
    int constraint = 0;
    amrex::Real weight = 0.0;
};

struct PackedFaceRow
{
    int local_fab = 0;
    int i = 0;
    int j = 0;
    int k = 0;
    int begin = 0;
    int end = 0;
};

class ConstraintStencilOperator
{
public:
    ConstraintStencilOperator(FaceFields const& face_layout,
                              std::vector<ConstraintPoint> const& constraints);

    void interpolate(FaceFields const& face_values,
                     FaceOwnerMasks const& owner_masks,
                     std::vector<ConstraintPoint>& constraints,
                     bool into_response);

    void spread(LevelData const& state,
                std::vector<ConstraintPoint> const& constraints,
                MarkerField field,
                amrex::Real scale,
                FaceFields& face_rhs);

private:
    amrex::Array<std::vector<int>, AMREX_SPACEDIM> m_constraint_ids;
    amrex::Array<std::vector<int>, AMREX_SPACEDIM> m_constraint_offsets;
    amrex::Array<std::vector<PackedFaceRef>, AMREX_SPACEDIM> m_interpolation_entries;
    amrex::Array<std::vector<PackedFaceRow>, AMREX_SPACEDIM> m_face_rows;
    amrex::Array<std::vector<PackedConstraintContribution>, AMREX_SPACEDIM>
        m_spreading_entries;

#ifdef AMREX_USE_GPU
    amrex::Array<amrex::Gpu::DeviceVector<int>, AMREX_SPACEDIM> m_device_constraint_ids;
    amrex::Array<amrex::Gpu::DeviceVector<int>, AMREX_SPACEDIM> m_device_constraint_offsets;
    amrex::Array<amrex::Gpu::DeviceVector<PackedFaceRef>, AMREX_SPACEDIM>
        m_device_interpolation_entries;
    amrex::Array<amrex::Gpu::DeviceVector<PackedFaceRow>, AMREX_SPACEDIM> m_device_face_rows;
    amrex::Array<amrex::Gpu::DeviceVector<PackedConstraintContribution>, AMREX_SPACEDIM>
        m_device_spreading_entries;
    amrex::Gpu::DeviceVector<amrex::Real> m_device_constraint_values;
    amrex::Gpu::DeviceVector<amrex::Real> m_device_interpolation_sums;
#endif

    std::vector<amrex::Real> m_constraint_values;
    std::vector<amrex::Real> m_interpolation_sums;
};

amrex::Real
cellVolume(amrex::Geometry const& geom);

amrex::Real
minCellSize(amrex::Geometry const& geom);

amrex::IntVect
makeCell(int i, int j, int k);

int
cellK(amrex::IntVect const& iv);

std::array<amrex::Real, 3>
cellCenter(amrex::Geometry const& geom, amrex::IntVect const& iv);

std::array<amrex::Real, 3>
faceCenter(amrex::Geometry const& geom, amrex::IntVect const& iv, int component);

RigidSphere
optionSphere(DLMOptions const& options);

amrex::Real
distanceSquared(std::array<amrex::Real, 3> const& point, DLMOptions const& options);

amrex::Real
distanceSquared(std::array<amrex::Real, 3> const& point, BodyState const& body);

bool
sameCell(amrex::IntVect const& a, amrex::IntVect const& b);

int
domainIndex(amrex::Geometry const& geom, int dir, int idx);

void
addStencilEntry(std::vector<StencilEntry>& stencil, amrex::IntVect const& cell,
                amrex::Real weight);

std::vector<StencilEntry>
q1Stencil(amrex::Geometry const& geom, std::array<amrex::Real, 3> const& point);

int
domainStaggeredIndex(amrex::Geometry const& geom, int dir, int component, int idx);

std::vector<StencilEntry>
q1StaggeredStencil(amrex::Geometry const& geom,
                   std::array<amrex::Real, 3> const& point,
                   int component);

amrex::Real
q2Basis(int node, amrex::Real t);

bool
supportedDlmInterpolation(std::string const& interpolation);

bool
supportedMultiplierLayout(std::string const& layout);

std::vector<StencilEntry>
q2StaggeredStencil(amrex::Geometry const& geom,
                   std::array<amrex::Real, 3> const& point,
                   std::array<amrex::Real, 3> const& center_vector,
                   int component);

MarkerPoint
makeMarker(amrex::Geometry const& geom,
           DLMOptions const& options,
           RigidSphere const& sphere,
           std::size_t body_index,
           std::array<amrex::Real, 3> position,
           amrex::Real quadrature_weight,
           bool boundary);

void
appendBoundaryMarkers(amrex::Geometry const& geom,
                      DLMOptions const& options,
                      RigidSphere const& sphere,
                      std::size_t body_index,
                      std::vector<MarkerPoint>& markers);

void
appendInteriorMarkers(amrex::Geometry const& geom,
                      DLMOptions const& options,
                      RigidSphere const& sphere,
                      std::size_t body_index,
                      std::vector<MarkerPoint>& markers);

void
appendMarkers(amrex::Geometry const& geom,
              DLMOptions const& options,
              RigidSphere const& sphere,
              std::size_t body_index,
              std::vector<MarkerPoint>& markers);

std::vector<MarkerPoint>
buildMarkers(amrex::Geometry const& geom, DLMOptions const& options);

std::vector<MarkerPoint>
buildMarkers(amrex::Geometry const& geom,
             std::vector<RigidSphere> const& spheres,
             DLMOptions const& options);

ConstraintPoint
makeConstraint(amrex::Geometry const& geom,
               DLMOptions const& options,
               MarkerPoint const& marker,
               int component);

std::vector<ConstraintPoint>
buildConstraints(amrex::Geometry const& geom, DLMOptions const& options);

std::vector<ConstraintPoint>
buildConstraints(amrex::Geometry const& geom,
                 std::vector<MarkerPoint> const& markers,
                 std::vector<RigidSphere> const& spheres,
                 DLMOptions const& options);

std::vector<ConstraintPoint>
buildConstraints(amrex::Geometry const& geom,
                 std::vector<RigidSphere> const& spheres,
                 DLMOptions const& options);

void
writeConstraintCatalog(std::vector<ConstraintPoint> const& constraints,
                       std::string const& filename);

amrex::Real
bodyVolume(DLMOptions const& options);

amrex::Real
sphereVolume(amrex::Real radius);

void
initializeBodyState(BodyState& body, DLMOptions const& options);

void
initializeBodyState(BodyState& body,
                    RigidSphere const& sphere,
                    DLMOptions const& options);

std::vector<BodyState>
buildBodyStates(std::vector<RigidSphere> const& spheres, DLMOptions const& options);

void
ensureBodyState(BodyState& body, DLMOptions const& options);

void
syncPrimaryBody(State& dlm_state);

RigidMotion
bodyMotion(BodyState const& body);

amrex::Real
markerWeightAtCell(MarkerPoint const& marker, amrex::IntVect const& cell);

amrex::Real
constraintWeightAtFace(ConstraintPoint const& constraint, amrex::IntVect const& face);

void
defineFaceFieldsImpl(LevelData const& state, FaceFields& faces);

FaceOwnerMasks
makeFaceOwnerMasks(FaceFields const& faces, amrex::Geometry const& geom);

void
cellVectorToComponentFacesImpl(LevelData const& state,
                               amrex::MultiFab const& cell_vector,
                               FaceFields& faces);

void
componentFacesToCellVectorImpl(LevelData const& state,
                               FaceFields const& faces,
                               amrex::MultiFab& cell_vector);

void
fillRigidMask(LevelData& state, DLMOptions const& options);

void
fillRigidMask(LevelData& state,
              std::vector<BodyState> const& bodies,
              DLMOptions const& options);

RigidMotion
prescribedRigidMotion(DLMOptions const& options);

RigidMotion
prescribedRigidMotion(BodyState const& body);

RigidMotion
bestFitRigidMotionCells(LevelData const& state, DLMOptions const& options);

RigidMotion
bestFitRigidMotionMarkers(std::vector<MarkerPoint> const& markers);

RigidMotion
bestFitRigidMotionMarkers(std::vector<MarkerPoint> const& markers,
                          std::size_t body_index);

RigidMotion
rigidMotionCells(LevelData const& state, DLMOptions const& options);

RigidMotion
rigidMotionMarkers(std::vector<MarkerPoint> const& markers, DLMOptions const& options);

RigidMotion
rigidMotionMarkers(std::vector<MarkerPoint> const& markers,
                   BodyState const& body,
                   DLMOptions const& options);

std::array<amrex::Real, 3>
rigidVelocityAt(std::array<amrex::Real, 3> const& position,
                RigidMotion const& motion,
                DLMOptions const& options);

std::array<amrex::Real, 3>
rigidVelocityAt(std::array<amrex::Real, 3> const& position,
                RigidMotion const& motion,
                BodyState const& body);

void
interpolateMarkers(LevelData const& state, std::vector<MarkerPoint>& markers);

std::array<amrex::Real, 3>&
markerField(MarkerPoint& marker, MarkerField field);

std::array<amrex::Real, 3> const&
markerField(MarkerPoint const& marker, MarkerField field);

bool
hasNonzeroMarkerField(std::vector<MarkerPoint> const& markers, MarkerField field);

void
copyMarkerField(std::vector<MarkerPoint>& markers, MarkerField dst, MarkerField src);

void
resetMarkerUzawaFields(std::vector<MarkerPoint>& markers);

amrex::Real
markerDot(std::vector<MarkerPoint> const& markers, MarkerField a, MarkerField b);

GeneralizedForce
markerGeneralizedForce(std::vector<MarkerPoint> const& markers, MarkerField field);

amrex::Real&
constraintField(ConstraintPoint& constraint, MarkerField field);

amrex::Real const&
constraintField(ConstraintPoint const& constraint, MarkerField field);

bool
hasNonzeroConstraintField(std::vector<ConstraintPoint> const& constraints, MarkerField field);

void
copyConstraintField(std::vector<ConstraintPoint>& constraints, MarkerField dst, MarkerField src);

void
resetConstraintUzawaFields(std::vector<ConstraintPoint>& constraints);

amrex::Real
constraintDot(std::vector<ConstraintPoint> const& constraints, MarkerField a, MarkerField b);

GeneralizedForce
constraintGeneralizedForce(std::vector<ConstraintPoint> const& constraints, MarkerField field);

GeneralizedForce
constraintGeneralizedForce(std::vector<ConstraintPoint> const& constraints,
                           MarkerField field,
                           std::size_t body_index);

GeneralizedForce
initialBodyRhs(BodyState const& body,
               DLMOptions const& options,
               GeneralizedForce marker_rhs,
               amrex::Real dt);

RigidMotion
solveBodyVelocity(BodyState const& body,
                  GeneralizedForce const& rhs,
                  amrex::Real dt);

RigidMotion
solveInitialBodyMotion(BodyState const& body,
                       DLMOptions const& options,
                       std::vector<ConstraintPoint> const& constraints,
                       amrex::Real dt);

RigidMotion
solveSearchBodyMotion(BodyState const& body,
                      DLMOptions const& options,
                      std::vector<ConstraintPoint> const& constraints,
                      amrex::Real dt);

void
setBodyMotion(BodyState& body, RigidMotion const& motion);

void
addBodyMotion(BodyState& body, RigidMotion const& motion, amrex::Real scale);

void
updateLambdaAndResidual(std::vector<MarkerPoint>& markers, amrex::Real alpha);

void
updateSearch(std::vector<MarkerPoint>& markers, amrex::Real beta);

void
updateLambdaAndResidual(std::vector<ConstraintPoint>& constraints, amrex::Real alpha);

void
updateSearch(std::vector<ConstraintPoint>& constraints, amrex::Real beta);

void
spreadMarkerRhs(LevelData const& state,
                std::vector<MarkerPoint> const& markers,
                MarkerField field,
                amrex::Real scale,
                amrex::MultiFab& rhs_delta);

void
addDeltaToVelocity(LevelData& state,
                   amrex::MultiFab const& delta,
                   amrex::MultiFab& cumulative_delta,
                   amrex::Real scale);

void
addFaceDeltaToVelocity(LevelData& state,
                       FaceFields const& face_delta,
                       amrex::MultiFab& cell_delta,
                       amrex::MultiFab& cumulative_delta,
                       amrex::Real scale);

void
interpolateMarkerResponse(amrex::MultiFab const& delta,
                          std::vector<MarkerPoint>& markers);

void
subtractRigidResponseFromMarkers(std::vector<MarkerPoint>& markers,
                                 RigidMotion const& motion,
                                 DLMOptions const& options);

void
subtractRigidResponseFromConstraints(std::vector<ConstraintPoint>& constraints,
                                     RigidMotion const& motion,
                                     DLMOptions const& options);

void
subtractRigidResponseFromConstraints(std::vector<ConstraintPoint>& constraints,
                                     std::vector<RigidMotion> const& motions,
                                     std::vector<BodyState> const& bodies);

void
fillMarkerResiduals(std::vector<MarkerPoint>& markers,
                    RigidMotion const& motion,
                    DLMOptions const& options);

void
fillConstraintResiduals(std::vector<ConstraintPoint>& constraints,
                        RigidMotion const& motion,
                        DLMOptions const& options);

void
fillConstraintResiduals(std::vector<ConstraintPoint>& constraints,
                        std::vector<BodyState> const& bodies);

void
fillRigidVelocityDiagnostic(LevelData& state,
                            RigidMotion const& motion,
                            DLMOptions const& options);

void
fillRigidVelocityDiagnostic(LevelData& state,
                            std::vector<BodyState> const& bodies,
                            DLMOptions const& options);

void
enforceInteriorRigidVelocityCells(LevelData& state,
                                  std::vector<BodyState> const& bodies,
                                  amrex::MultiFab& cumulative_delta);

void
setDlmForceFromDelta(LevelData& state, amrex::MultiFab const& cumulative_delta, amrex::Real dt);

void
addDlmForceFromDelta(LevelData& state, amrex::MultiFab const& cumulative_delta, amrex::Real dt);

void
setHydrodynamicLoads(State& dlm_state, amrex::Real dt);

void
applyRigidVelocityCells(LevelData& state,
                        RigidMotion const& motion,
                        DLMOptions const& options,
                        amrex::Real dt);

void
applyCellMaskProjection(LevelData& state, DLMOptions const& options, amrex::Real dt);

void
ensureMarkerState(State& dlm_state, LevelData const& state, DLMOptions const& options);

void
ensureMarkerState(State& dlm_state,
                  LevelData const& state,
                  std::vector<RigidSphere> const& spheres,
                  DLMOptions const& options);

void
applyMarkerProjection(State& dlm_state,
                      LevelData& state,
                      DLMOptions const& options,
                      amrex::Real dt,
                      StaggeredFluidResponse const& fluid_response);

void
applyMarkerProjection(State& dlm_state,
                      LevelData& state,
                      std::vector<RigidSphere> const& spheres,
                      DLMOptions const& options,
                      amrex::Real dt,
                      StaggeredFluidResponse const& fluid_response);

} // namespace ELFF::Fluid::DLM::Detail

#endif
