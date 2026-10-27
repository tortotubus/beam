#include "DLMInternal.H"

#include <AMReX_BLProfiler.H>

#include <tuple>

namespace ELFF::Fluid::DLM::Detail {

amrex::Real
markerWeightAtCell(MarkerPoint const& marker, amrex::IntVect const& cell)
{
    for (StencilEntry const& entry : marker.stencil) {
        if (sameCell(entry.cell, cell)) {
            return entry.weight;
        }
    }
    return 0.0;
}

amrex::Real
constraintWeightAtFace(ConstraintPoint const& constraint, amrex::IntVect const& face)
{
    for (StencilEntry const& entry : constraint.stencil) {
        if (sameCell(entry.cell, face)) {
            return entry.weight;
        }
    }
    return 0.0;
}

void
defineFaceFieldsImpl(LevelData const& state, FaceFields& faces)
{
    for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
        faces[comp].define(state.faceVelocity()[comp].boxArray(),
                           state.faceVelocity()[comp].DistributionMap(),
                           1,
                           0,
                           amrex::MFInfo(),
                           state.faceVelocity()[comp].Factory());
    }
}

FaceOwnerMasks
makeFaceOwnerMasks(FaceFields const& faces, amrex::Geometry const& geom)
{
    FaceOwnerMasks masks;
    for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
        masks[comp] = faces[comp].OwnerMask(geom.periodicity());
    }
    return masks;
}

void
cellVectorToComponentFacesImpl(LevelData const& state,
                               amrex::MultiFab const& cell_vector,
                               FaceFields& faces)
{
    BL_PROFILE("ELFF::Fluid::DLM::cellVectorToComponentFaces");
    const amrex::Geometry& geom = state.geom();
    const amrex::Box& domain = geom.Domain();

    for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
        amrex::BoxArray cell_grids_on_face_layout =
            amrex::convert(faces[comp].boxArray(), amrex::IntVect::TheZeroVector());
        amrex::MultiFab cell_vector_on_face_layout(cell_grids_on_face_layout,
                                                   faces[comp].DistributionMap(),
                                                   AMREX_SPACEDIM,
                                                   1,
                                                   amrex::MFInfo(),
                                                   cell_vector.Factory());
        cell_vector_on_face_layout.setVal(0.0,
                                          0,
                                          AMREX_SPACEDIM,
                                          cell_vector_on_face_layout.nGrowVect());
        cell_vector_on_face_layout.ParallelCopy(cell_vector,
                                                0,
                                                0,
                                                AMREX_SPACEDIM,
                                                0,
                                                1,
                                                geom.periodicity());

        faces[comp].setVal(0.0);
        for (amrex::MFIter mfi(faces[comp], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.tilebox();
            auto const face = faces[comp].array(mfi);
            auto const vel = cell_vector_on_face_layout.const_array(mfi);
            const int lo = domain.smallEnd(comp);
            const int hi = domain.bigEnd(comp);
            const bool periodic = geom.isPeriodic(comp);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                int li = i;
                int lj = j;
                int lk = k;
                int ri = i;
                int rj = j;
                int rk = k;

                if (comp == 0) {
                    --li;
                    if (!periodic) {
                        if (li < lo) {
                            li = lo;
                        }
                        if (ri > hi) {
                            ri = hi;
                        }
                    }
                } else if (comp == 1) {
                    --lj;
                    if (!periodic) {
                        if (lj < lo) {
                            lj = lo;
                        }
                        if (rj > hi) {
                            rj = hi;
                        }
                    }
#if (AMREX_SPACEDIM == 3)
                } else {
                    --lk;
                    if (!periodic) {
                        if (lk < lo) {
                            lk = lo;
                        }
                        if (rk > hi) {
                            rk = hi;
                        }
                    }
#endif
                }

                face(i,j,k) = amrex::Real(0.5) * (vel(li,lj,lk,comp) + vel(ri,rj,rk,comp));
            });
        }
    }
}

void
componentFacesToCellVectorImpl(LevelData const& state,
                               FaceFields const& faces,
                               amrex::MultiFab& cell_vector)
{
    BL_PROFILE("ELFF::Fluid::DLM::componentFacesToCellVector");
    amrex::Array<amrex::MultiFab const*, AMREX_SPACEDIM> face_ptrs{
        AMREX_D_DECL(&faces[0], &faces[1], &faces[2])
    };
    cell_vector.setVal(0.0);
    amrex::average_face_to_cellcenter(cell_vector, face_ptrs, state.geom());
}

void
interpolateMarkers(LevelData const& state, std::vector<MarkerPoint>& markers)
{
    const std::size_t nvalues = markers.size() * 3;
    std::vector<amrex::Real> sums(nvalues, 0.0);

    for (amrex::MFIter mfi(state.velocity()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto const vel = state.velocity().const_array(mfi);
        for (amrex::BoxIterator bit(bx); bit.ok(); ++bit) {
            const amrex::IntVect iv = bit();
            const int i = iv[0];
            const int j = iv[1];
            const int k = cellK(iv);
            for (std::size_t marker_id = 0; marker_id < markers.size(); ++marker_id) {
                const amrex::Real weight = markerWeightAtCell(markers[marker_id], iv);
                if (weight <= 0.0) {
                    continue;
                }
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                    sums[3 * marker_id + dir] += weight * vel(i,j,k,dir);
                }
            }
        }
    }

    if (!sums.empty()) {
        amrex::ParallelDescriptor::ReduceRealSum(sums.data(), static_cast<int>(sums.size()));
    }

    for (std::size_t marker_id = 0; marker_id < markers.size(); ++marker_id) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            markers[marker_id].interpolated_velocity[dir] = sums[3 * marker_id + dir];
        }
    }
}

void
spreadMarkerRhs(LevelData const& state,
                std::vector<MarkerPoint> const& markers,
                MarkerField field,
                amrex::Real scale,
                amrex::MultiFab& rhs_delta)
{
    rhs_delta.setVal(0.0);
    const amrex::Real inv_cell_volume = amrex::Real(1.0) / cellVolume(state.geom());
    for (amrex::MFIter mfi(rhs_delta); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto const rhs = rhs_delta.array(mfi);
        for (amrex::BoxIterator bit(bx); bit.ok(); ++bit) {
            const amrex::IntVect iv = bit();
            std::array<amrex::Real, 3> force{0.0, 0.0, 0.0};

            for (MarkerPoint const& marker : markers) {
                const amrex::Real stencil_weight = markerWeightAtCell(marker, iv);
                if (stencil_weight <= 0.0) {
                    continue;
                }
                auto const& values = markerField(marker, field);
                const amrex::Real weight =
                    scale * marker.quadrature_weight * stencil_weight * inv_cell_volume;
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                    force[dir] += weight * values[dir];
                }
            }

            const int i = iv[0];
            const int j = iv[1];
            const int k = cellK(iv);
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                rhs(i,j,k,dir) = force[dir];
            }
        }
    }
}

void
addDeltaToVelocity(LevelData& state,
                   amrex::MultiFab const& delta,
                   amrex::MultiFab& cumulative_delta,
                   amrex::Real scale)
{
    for (amrex::MFIter mfi(state.velocity(), amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto const vel = state.velocity().array(mfi);
        auto const delta_arr = delta.const_array(mfi);
        auto const cumulative = cumulative_delta.array(mfi);
        amrex::ParallelFor(bx, AMREX_SPACEDIM,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
            const amrex::Real update = scale * delta_arr(i,j,k,n);
            vel(i,j,k,n) += update;
            cumulative(i,j,k,n) += update;
        });
    }
}

void
addFaceDeltaToVelocity(LevelData& state,
                       FaceFields const& face_delta,
                       amrex::MultiFab& cell_delta,
                       amrex::MultiFab& cumulative_delta,
                       amrex::Real scale)
{
    componentFacesToCellVectorImpl(state, face_delta, cell_delta);
    addDeltaToVelocity(state, cell_delta, cumulative_delta, scale);
}

void
interpolateMarkerResponse(amrex::MultiFab const& delta,
                          std::vector<MarkerPoint>& markers)
{
    const std::size_t nvalues = markers.size() * 3;
    std::vector<amrex::Real> sums(nvalues, 0.0);

    for (amrex::MFIter mfi(delta); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto const delta_arr = delta.const_array(mfi);
        for (amrex::BoxIterator bit(bx); bit.ok(); ++bit) {
            const amrex::IntVect iv = bit();
            const int i = iv[0];
            const int j = iv[1];
            const int k = cellK(iv);
            for (std::size_t marker_id = 0; marker_id < markers.size(); ++marker_id) {
                const amrex::Real weight = markerWeightAtCell(markers[marker_id], iv);
                if (weight <= 0.0) {
                    continue;
                }
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                    sums[3 * marker_id + dir] += weight * delta_arr(i,j,k,dir);
                }
            }
        }
    }

    if (!sums.empty()) {
        amrex::ParallelDescriptor::ReduceRealSum(sums.data(), static_cast<int>(sums.size()));
    }

    for (std::size_t marker_id = 0; marker_id < markers.size(); ++marker_id) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            markers[marker_id].response[dir] = sums[3 * marker_id + dir];
        }
    }
}

namespace {

struct PendingFaceContribution
{
    int local_fab = 0;
    int i = 0;
    int j = 0;
    int k = 0;
    int constraint = 0;
    amrex::Real weight = 0.0;
};

bool
sameFace(PendingFaceContribution const& a, PendingFaceContribution const& b)
{
    return a.local_fab == b.local_fab && a.i == b.i && a.j == b.j && a.k == b.k;
}

#ifdef AMREX_USE_GPU
template <typename T>
void
copyToDevice(std::vector<T> const& host, amrex::Gpu::DeviceVector<T>& device)
{
    device.resize(host.size());
    if (!host.empty()) {
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                              host.begin(), host.end(), device.begin());
    }
}
#endif

} // namespace

ConstraintStencilOperator::ConstraintStencilOperator(
    FaceFields const& face_layout,
    std::vector<ConstraintPoint> const& constraints)
    : m_constraint_values(constraints.size(), 0.0),
      m_interpolation_sums(constraints.size(), 0.0)
{
    BL_PROFILE("ELFF::Fluid::DLM::ConstraintStencilOperator::build");
    std::vector<std::vector<PackedFaceRef>> constraint_entries(constraints.size());
    amrex::Array<std::vector<PendingFaceContribution>, AMREX_SPACEDIM> pending;

    // Build the topology once.  A staggered face can occur in two local valid
    // FABs.  Keeping both copies here is intentional for J^T; OwnerMask later
    // selects exactly one copy for J.
    for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
        for (amrex::MFIter mfi(face_layout[comp]); mfi.isValid(); ++mfi) {
            amrex::Box const& bx = mfi.validbox();
            const int local_fab = mfi.LocalIndex();
            for (std::size_t constraint_id = 0; constraint_id < constraints.size();
                 ++constraint_id) {
                ConstraintPoint const& constraint = constraints[constraint_id];
                if (constraint.component != comp) {
                    continue;
                }
                for (StencilEntry const& entry : constraint.stencil) {
                    if (entry.weight == 0.0 || !bx.contains(entry.cell)) {
                        continue;
                    }
                    PackedFaceRef const ref{local_fab,
                                            entry.cell[0],
                                            entry.cell[1],
                                            cellK(entry.cell),
                                            entry.weight};
                    constraint_entries[constraint_id].push_back(ref);
                    pending[comp].push_back(PendingFaceContribution{
                        ref.local_fab, ref.i, ref.j, ref.k,
                        static_cast<int>(constraint_id), ref.weight});
                }
            }
        }
    }

    for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
        m_constraint_offsets[comp].push_back(0);
        for (std::size_t constraint_id = 0; constraint_id < constraints.size();
             ++constraint_id) {
            if (constraints[constraint_id].component != comp) {
                continue;
            }
            m_constraint_ids[comp].push_back(static_cast<int>(constraint_id));
            auto const& entries = constraint_entries[constraint_id];
            m_interpolation_entries[comp].insert(m_interpolation_entries[comp].end(),
                                                  entries.begin(), entries.end());
            m_constraint_offsets[comp].push_back(
                static_cast<int>(m_interpolation_entries[comp].size()));
        }

        auto& contributions = pending[comp];
        std::sort(contributions.begin(), contributions.end(),
                  [](PendingFaceContribution const& a,
                     PendingFaceContribution const& b) {
                      return std::tie(a.local_fab, a.k, a.j, a.i, a.constraint)
                           < std::tie(b.local_fab, b.k, b.j, b.i, b.constraint);
                  });

        std::size_t begin = 0;
        while (begin < contributions.size()) {
            std::size_t end = begin + 1;
            while (end < contributions.size()
                   && sameFace(contributions[begin], contributions[end])) {
                ++end;
            }
            const int row_begin = static_cast<int>(m_spreading_entries[comp].size());
            for (std::size_t index = begin; index < end; ++index) {
                m_spreading_entries[comp].push_back(PackedConstraintContribution{
                    contributions[index].constraint, contributions[index].weight});
            }
            m_face_rows[comp].push_back(PackedFaceRow{
                contributions[begin].local_fab,
                contributions[begin].i,
                contributions[begin].j,
                contributions[begin].k,
                row_begin,
                static_cast<int>(m_spreading_entries[comp].size())});
            begin = end;
        }

#ifdef AMREX_USE_GPU
        copyToDevice(m_constraint_ids[comp], m_device_constraint_ids[comp]);
        copyToDevice(m_constraint_offsets[comp], m_device_constraint_offsets[comp]);
        copyToDevice(m_interpolation_entries[comp], m_device_interpolation_entries[comp]);
        copyToDevice(m_face_rows[comp], m_device_face_rows[comp]);
        copyToDevice(m_spreading_entries[comp], m_device_spreading_entries[comp]);
#endif
    }

#ifdef AMREX_USE_GPU
    m_device_constraint_values.resize(constraints.size());
    m_device_interpolation_sums.resize(constraints.size());
#endif
}

void
ConstraintStencilOperator::interpolate(
    FaceFields const& face_values,
    FaceOwnerMasks const& owner_masks,
    std::vector<ConstraintPoint>& constraints,
    bool into_response)
{
    BL_PROFILE("ELFF::Fluid::DLM::ConstraintStencilOperator::interpolate");
    AMREX_ALWAYS_ASSERT(constraints.size() == m_interpolation_sums.size());
    std::fill(m_interpolation_sums.begin(), m_interpolation_sums.end(), amrex::Real(0.0));

#ifdef AMREX_USE_GPU
    if (amrex::Gpu::inLaunchRegion()) {
        if (!m_interpolation_sums.empty()) {
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                                  m_interpolation_sums.begin(), m_interpolation_sums.end(),
                                  m_device_interpolation_sums.begin());
        }
        for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
            const int nrows = static_cast<int>(m_constraint_ids[comp].size());
            if (nrows == 0) {
                continue;
            }
            auto const faces = face_values[comp].const_arrays();
            auto const owners = owner_masks[comp]->const_arrays();
            int const* ids = m_device_constraint_ids[comp].dataPtr();
            int const* offsets = m_device_constraint_offsets[comp].dataPtr();
            PackedFaceRef const* entries = m_device_interpolation_entries[comp].dataPtr();
            amrex::Real* sums = m_device_interpolation_sums.dataPtr();
            amrex::ParallelFor(nrows, [=] AMREX_GPU_DEVICE (int row) noexcept {
                amrex::Real sum = 0.0;
                for (int index = offsets[row]; index < offsets[row + 1]; ++index) {
                    PackedFaceRef const& entry = entries[index];
                    if (owners[entry.local_fab](entry.i, entry.j, entry.k) != 0) {
                        sum += entry.weight
                             * faces[entry.local_fab](entry.i, entry.j, entry.k);
                    }
                }
                sums[ids[row]] = sum;
            });
        }
        if (!m_interpolation_sums.empty()) {
            amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                             m_device_interpolation_sums.begin(),
                             m_device_interpolation_sums.end(),
                             m_interpolation_sums.begin());
        }
    } else
#endif
    {
        for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
            auto const faces = face_values[comp].const_arrays();
            auto const owners = owner_masks[comp]->const_arrays();
            auto const& ids = m_constraint_ids[comp];
            auto const& offsets = m_constraint_offsets[comp];
            auto const& entries = m_interpolation_entries[comp];
            const int nrows = static_cast<int>(ids.size());
#ifdef AMREX_USE_OMP
#pragma omp parallel for schedule(static)
#endif
            for (int row = 0; row < nrows; ++row) {
                amrex::Real sum = 0.0;
                for (int index = offsets[row]; index < offsets[row + 1]; ++index) {
                    PackedFaceRef const& entry = entries[index];
                    if (owners[entry.local_fab](entry.i, entry.j, entry.k) != 0) {
                        sum += entry.weight
                             * faces[entry.local_fab](entry.i, entry.j, entry.k);
                    }
                }
                m_interpolation_sums[ids[row]] = sum;
            }
        }
    }

    if (!m_interpolation_sums.empty()) {
        amrex::ParallelDescriptor::ReduceRealSum(
            m_interpolation_sums.data(), static_cast<int>(m_interpolation_sums.size()));
    }
    for (std::size_t constraint_id = 0; constraint_id < constraints.size(); ++constraint_id) {
        if (into_response) {
            constraints[constraint_id].response = m_interpolation_sums[constraint_id];
        } else {
            constraints[constraint_id].interpolated_value =
                m_interpolation_sums[constraint_id];
        }
    }
}

void
ConstraintStencilOperator::spread(
    LevelData const& state,
    std::vector<ConstraintPoint> const& constraints,
    MarkerField field,
    amrex::Real scale,
    FaceFields& face_rhs)
{
    BL_PROFILE("ELFF::Fluid::DLM::ConstraintStencilOperator::spread");
    AMREX_ALWAYS_ASSERT(constraints.size() == m_constraint_values.size());
    const amrex::Real inv_cell_volume = amrex::Real(1.0) / cellVolume(state.geom());
    for (std::size_t constraint_id = 0; constraint_id < constraints.size(); ++constraint_id) {
        ConstraintPoint const& constraint = constraints[constraint_id];
        m_constraint_values[constraint_id] = scale * constraint.quadrature_weight
                                           * inv_cell_volume
                                           * constraintField(constraint, field);
    }

    for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
        face_rhs[comp].setVal(0.0);
    }

#ifdef AMREX_USE_GPU
    if (amrex::Gpu::inLaunchRegion()) {
        if (!m_constraint_values.empty()) {
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                                  m_constraint_values.begin(), m_constraint_values.end(),
                                  m_device_constraint_values.begin());
        }
        amrex::Real const* values = m_device_constraint_values.dataPtr();
        for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
            const int nrows = static_cast<int>(m_face_rows[comp].size());
            if (nrows == 0) {
                continue;
            }
            auto const rhs = face_rhs[comp].arrays();
            PackedFaceRow const* rows = m_device_face_rows[comp].dataPtr();
            PackedConstraintContribution const* entries =
                m_device_spreading_entries[comp].dataPtr();
            amrex::ParallelFor(nrows, [=] AMREX_GPU_DEVICE (int row_index) noexcept {
                PackedFaceRow const& row = rows[row_index];
                amrex::Real sum = 0.0;
                for (int index = row.begin; index < row.end; ++index) {
                    PackedConstraintContribution const& entry = entries[index];
                    sum += entry.weight * values[entry.constraint];
                }
                rhs[row.local_fab](row.i, row.j, row.k) = sum;
            });
        }
    } else
#endif
    {
        for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
            auto const rhs = face_rhs[comp].arrays();
            auto const& rows = m_face_rows[comp];
            auto const& entries = m_spreading_entries[comp];
            const int nrows = static_cast<int>(rows.size());
#ifdef AMREX_USE_OMP
#pragma omp parallel for schedule(static)
#endif
            for (int row_index = 0; row_index < nrows; ++row_index) {
                PackedFaceRow const& row = rows[row_index];
                amrex::Real sum = 0.0;
                for (int index = row.begin; index < row.end; ++index) {
                    PackedConstraintContribution const& entry = entries[index];
                    sum += entry.weight * m_constraint_values[entry.constraint];
                }
                rhs[row.local_fab](row.i, row.j, row.k) = sum;
            }
        }
    }
}

} // namespace ELFF::Fluid::DLM::Detail
