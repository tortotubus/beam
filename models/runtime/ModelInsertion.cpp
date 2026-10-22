#include "elff/models/runtime/ModelInsertion.hpp"
#include "elff/contact/NarrowPhaseSegmentDistance.hpp"
#include "elff/models/ibm/IBCoupling.hpp"
#include <cmath>
#include <cstring>

namespace ELFF::Models {
using namespace Contact;

namespace {

InsertionResult failure(InsertionStatus status) { return {status}; }

/** */
AABB bounds(const LineSegment &segment) {
  AABB box;
  for (size_t d = 0; d < 3; ++d) {
    box.min[d] = std::min(segment.start[d], segment.end[d]);
    box.max[d] = std::max(segment.start[d], segment.end[d]);
  }
  return box;
}

} // namespace

bool ModelInsertion::available() const {
  return configured_ && !sealed_ && runtime_.mutable_registration();
}

InsertionStatus ModelInsertion::configure(const InsertionGridConfig &config,
    std::optional<std::array<real_t, 3>> exact_upper) {
  if (configured_ || sealed_ || !runtime_.mutable_registration())
    return InsertionStatus::InvalidLifecycle;
  if (!std::isfinite(config.cell_width) || config.cell_width <= 0)
    return InsertionStatus::InvalidInput;
  size_t count = 1;
  std::array<real_t, 3> upper;
  // Match the grid's storage bound without constructing invalid grid geometry.
  const auto limit = std::vector<std::vector<ShapeId>>().max_size();
  for (size_t d = 0; d < 3; ++d) {
    if (!config.dimensions[d] || config.dimensions[d] > limit / count)
      return InsertionStatus::InvalidInput;
    count *= config.dimensions[d];
    upper[d] = config.origin[d] + config.cell_width * config.dimensions[d];
    if (!std::isfinite(config.origin[d]) || !std::isfinite(upper[d]) ||
        upper[d] <= config.origin[d])
      return InsertionStatus::InvalidInput;
  }
  if (exact_upper)
    for (size_t d = 0; d < 3; ++d)
      if (!std::isfinite((*exact_upper)[d]) || (*exact_upper)[d] < config.origin[d] ||
          (*exact_upper)[d] > upper[d]) return InsertionStatus::InvalidInput;
  geometry_upper_ = exact_upper.value_or(upper);
  config_ = config;
  upper_ = upper;
  configured_ = true;
  return InsertionStatus::Accepted;
}

InsertionStatus ModelInsertion::validate(PolylineView geometry) const {
  if (!geometry.points || geometry.count < 2)
    return InsertionStatus::InvalidInput;
  for (size_t i = 0; i < geometry.count; ++i)
    for (size_t d = 0; d < 3; ++d) {
      const auto value = geometry.points[i][d];
      if (!std::isfinite(value))
        return InsertionStatus::InvalidInput;
      if (value < config_.origin[d] || value > geometry_upper_[d])
        return InsertionStatus::OutOfDomain;
    }
  return InsertionStatus::Accepted;
}

InsertionResult
ModelInsertion::check_clearance(const ContactGeometryProvider &provider,
                                real_t clearance) const {
  if (!available())
    return failure(InsertionStatus::InvalidLifecycle);
  return query(provider.contact_geometry(), clearance);
}

InsertionResult ModelInsertion::query(PolylineView geometry,
                                      real_t clearance) const {
  if (!std::isfinite(clearance) || clearance < 0)
    return failure(InsertionStatus::InvalidInput);
  const auto status = validate(geometry);
  if (status != InsertionStatus::Accepted)
    return failure(status);
  BroadPhaseLinkedCell grid(config_.origin, config_.cell_width,
                            config_.dimensions);
  struct Segment {
    LineSegment line;
    ModelId model;
    size_t index;
  };
  std::vector<Segment> segments;
  for (const auto &snapshot : snapshots_)
    for (size_t i = 0; i + 1 < snapshot.points.size(); ++i) {
      const LineSegment line{snapshot.points[i], snapshot.points[i + 1]};
      grid.insert(segments.size(), bounds(line));
      segments.push_back({line, snapshot.id, i});
    }
  const size_t first_candidate = segments.size();
  for (size_t i = 0; i + 1 < geometry.count; ++i) {
    const LineSegment line{geometry.points[i], geometry.points[i + 1]};
    auto box = bounds(line);
    for (size_t d = 0; d < 3; ++d) {
      // Wider arithmetic and outward rounding keep expanded bounds
      // conservative.
      const long double lo = static_cast<long double>(box.min[d]) - clearance;
      const long double hi = static_cast<long double>(box.max[d]) + clearance;
      box.min[d] =
          lo <= config_.origin[d]
              ? config_.origin[d]
              : std::max(
                    config_.origin[d],
                    std::nextafter(static_cast<real_t>(lo),
                                   -std::numeric_limits<real_t>::infinity()));
      box.max[d] = hi >= upper_[d]
                       ? upper_[d]
                       : std::min(upper_[d],
                                  std::nextafter(
                                      static_cast<real_t>(hi),
                                      std::numeric_limits<real_t>::infinity()));
    }
    grid.insert(segments.size(), box);
    segments.push_back({line, invalid_model_id, i});
  }
  NarrowPhaseSegmentDistance narrow;
  for (const auto &pair : grid.candidate_pairs(first_candidate)) {
    const auto &existing = segments[pair.first];
    const auto &candidate = segments[pair.second];
    if (narrow.within_clearance(candidate.line, existing.line, clearance))
      return {InsertionStatus::ClearanceViolation, invalid_model_id,
              existing.model, candidate.index, existing.index};
  }
  return {};
}

InsertionStatus ModelInsertion::register_model(ModelId id) {
  if (!available() || runtime_.ranks() != 1)
    return InsertionStatus::InvalidLifecycle;
  for (const auto &snapshot : snapshots_)
    if (snapshot.id == id)
      return InsertionStatus::InvalidInput;
  const auto *model = runtime_.model(id);
  if (!model)
    return InsertionStatus::InvalidInput;
  const auto *provider = dynamic_cast<const ContactGeometryProvider *>(model);
  if (!provider)
    return InsertionStatus::UnsupportedGeometry;
  const auto geometry = provider->contact_geometry();
  const auto status = validate(geometry);
  if (status != InsertionStatus::Accepted)
    return status;
  Snapshot prepared{id, {geometry.points, geometry.points + geometry.count}};
  snapshots_.push_back(std::move(prepared));
  return InsertionStatus::Accepted;
}

InsertionResult ModelInsertion::try_adopt(ModelType type,
                                          std::unique_ptr<Model> candidate,
                                          real_t clearance) {
  return adopt(type, std::move(candidate), clearance, nullptr, {});
}

InsertionResult ModelInsertion::try_adopt_ib(IBCoupling &coupling,
                                             ModelType type,
                                             std::unique_ptr<Model> candidate,
                                             real_t clearance,
                                             ib_model_info_t info) {
  return adopt(type, std::move(candidate), clearance, &coupling, info);
}

InsertionResult ModelInsertion::adopt(ModelType type,
                                      std::unique_ptr<Model> candidate,
                                      real_t clearance, IBCoupling *coupling,
                                      ib_model_info_t info) {
  if (!available() || runtime_.ranks() != 1)
    return failure(InsertionStatus::InvalidLifecycle);
  if (!candidate || type == ModelType::Invalid)
    return failure(InsertionStatus::InvalidInput);
  if (coupling) {
    if (&coupling->owning_runtime() != &runtime_)
      return failure(InsertionStatus::InvalidInput);
    auto *ib = dynamic_cast<IBModel *>(candidate.get());
    if (!ib ||
        ((dynamic_cast<IBForceCoupled *>(ib) != nullptr) ==
         (dynamic_cast<IBVelocityCoupled *>(ib) != nullptr)) ||
        (info.sparse != 0 && info.sparse != 1) || info.level < -1 ||
        info.node_count < 0)
      return failure(InsertionStatus::InvalidInput);
  }
  const auto *provider =
      dynamic_cast<const ContactGeometryProvider *>(candidate.get());
  if (!provider)
    return failure(InsertionStatus::UnsupportedGeometry);
  const auto geometry = provider->contact_geometry();
  auto result = query(geometry, clearance);
  if (!result.accepted())
    return result;
  const ModelId id = runtime_.size();
  Snapshot prepared{id, {geometry.points, geometry.points + geometry.count}};
  snapshots_.reserve(snapshots_.size() + 1);
  if (coupling) {
    auto *ib = dynamic_cast<IBModel *>(candidate.get());
    candidate.release();
    if (!coupling->adopt_owned_model(type, std::unique_ptr<IBModel>(ib), info))
      return failure(InsertionStatus::InvalidInput);
  } else if (!runtime_.adopt_owned_model(type, std::move(candidate))) {
    return failure(InsertionStatus::InvalidInput);
  }
  snapshots_.push_back(std::move(prepared));
  result.model_id = id;
  return result;
}

InsertionStatus ModelInsertion::agree_status(InsertionStatus local) const {
#ifdef ELFF_USE_MPI
  if (runtime_.ranks() > 1) {
    int low = static_cast<int>(local), high = low;
    MPI_Allreduce(MPI_IN_PLACE, &low, 1, MPI_INT, MPI_MIN, runtime_.communicator());
    MPI_Allreduce(MPI_IN_PLACE, &high, 1, MPI_INT, MPI_MAX, runtime_.communicator());
    if (high == static_cast<int>(InsertionStatus::PreparationFailed))
      return InsertionStatus::PreparationFailed;
    if (low != high) return InsertionStatus::CollectiveMismatch;
  }
#endif
  return local;
}

InsertionStatus ModelInsertion::register_model_collective(ModelId id) {
  auto status = available() ? InsertionStatus::Accepted : InsertionStatus::InvalidLifecycle;
  const auto* entry = runtime_.entry(id);
  if (status == InsertionStatus::Accepted) {
    if (!entry) status = InsertionStatus::InvalidInput;
    for (const auto& snapshot : snapshots_)
      if (snapshot.id == id) status = InsertionStatus::InvalidInput;
  }
  status = agree_status(status);
  if (status != InsertionStatus::Accepted) return status;
#ifdef ELFF_USE_MPI
  if (runtime_.ranks() > 1) {
    unsigned long long values[3] = {id, static_cast<unsigned long long>(entry->owner_pid), snapshots_.size()};
    auto same = true;
    for (auto value : values) {
      auto low = value, high = value;
      MPI_Allreduce(MPI_IN_PLACE, &low, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, runtime_.communicator());
      MPI_Allreduce(MPI_IN_PLACE, &high, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, runtime_.communicator());
      same = same && low == high;
    }
    if (!same) return InsertionStatus::CollectiveMismatch;
  }
#endif
  Snapshot prepared{id, {}};
  try {
    if (runtime_.rank() == entry->owner_pid) {
      auto* provider = dynamic_cast<const ContactGeometryProvider*>(entry->model.get());
      if (!provider) status = InsertionStatus::UnsupportedGeometry;
      else {
        auto geometry = provider->contact_geometry();
        status = validate(geometry);
        if (status == InsertionStatus::Accepted)
          prepared.points.assign(geometry.points, geometry.points + geometry.count);
      }
    }
  } catch (...) { status = InsertionStatus::PreparationFailed; }
  // Propagate owner validation before broadcasting its data.
#ifdef ELFF_USE_MPI
  if (runtime_.ranks() > 1) {
    int owner_status = static_cast<int>(status);
    MPI_Bcast(&owner_status, 1, MPI_INT, entry->owner_pid, runtime_.communicator());
    status = static_cast<InsertionStatus>(owner_status);
  }
#endif
  if (status != InsertionStatus::Accepted) return status;
  unsigned long long count = prepared.points.size();
#ifdef ELFF_USE_MPI
  if (runtime_.ranks() > 1)
    MPI_Bcast(&count, 1, MPI_UNSIGNED_LONG_LONG, entry->owner_pid, runtime_.communicator());
#endif
  try {
    if (count > prepared.points.max_size()) status = InsertionStatus::InvalidInput;
    else { prepared.points.resize(count); snapshots_.reserve(snapshots_.size() + 1); }
  } catch (...) { status = InsertionStatus::PreparationFailed; }
  status = agree_status(status);
  if (status != InsertionStatus::Accepted) return status;
#ifdef ELFF_USE_MPI
  if (runtime_.ranks() > 1) {
    // Bounded chunks avoid MPI's int count limit without allocating scratch buffers.
    auto* bytes = reinterpret_cast<char*>(prepared.points.data());
    size_t total = prepared.points.size() * sizeof(prepared.points[0]);
    for (size_t offset = 0; offset < total;) {
      int n = static_cast<int>(std::min<size_t>(total - offset, 4096));
      MPI_Bcast(bytes + offset, n, MPI_BYTE, entry->owner_pid, runtime_.communicator());
      offset += n;
    }
  }
#endif
  status = agree_status(validate({prepared.points.data(), prepared.points.size()}));
  if (status == InsertionStatus::Accepted) snapshots_.push_back(std::move(prepared));
  return status;
}

InsertionResult ModelInsertion::try_adopt_replicated(
    ModelType type, const CandidateFactory &factory, real_t clearance,
    int owner, IBCoupling *coupling, ib_model_info_t info) {
  const bool distributed = runtime_.ranks() > 1;
  InsertionStatus setup = InsertionStatus::Accepted;
  if (!available() ||
      next_attempt_ == std::numeric_limits<std::uint64_t>::max())
    setup = InsertionStatus::InvalidLifecycle;
  else if (!factory || type == ModelType::Invalid || owner < 0 ||
           owner >= runtime_.ranks() || !std::isfinite(clearance) ||
           clearance < 0 ||
           (coupling && (&coupling->owning_runtime() != &runtime_ ||
                         !IBCoupling::valid_info(info))))
    setup = InsertionStatus::InvalidInput;
  setup = agree_status(setup);
  if (setup != InsertionStatus::Accepted)
    return failure(setup);
#ifdef ELFF_USE_MPI
  if (distributed) {
    auto low = next_attempt_, high = next_attempt_;
    MPI_Allreduce(MPI_IN_PLACE, &low, 1, MPI_UINT64_T, MPI_MIN,
                  runtime_.communicator());
    MPI_Allreduce(MPI_IN_PLACE, &high, 1, MPI_UINT64_T, MPI_MAX,
                  runtime_.communicator());
    if (low != high)
      return failure(InsertionStatus::CollectiveMismatch);
  }
#endif
  const auto attempt = next_attempt_++;
  InsertionResult result;
  std::unique_ptr<Model> candidate;
  Snapshot prepared;
  std::vector<char> signature;
  ib_model_view_kind_t kind = IB_MODEL_VIEW_INVALID;
  try {
    candidate = factory(random_, attempt);
    if (!candidate)
      result.status = InsertionStatus::InvalidInput;
    else {
      const auto *provider =
          dynamic_cast<const ContactGeometryProvider *>(candidate.get());
      if (!provider)
        result.status = InsertionStatus::UnsupportedGeometry;
      else {
        const auto geometry = provider->contact_geometry();
        result = query(geometry, clearance);
        if (coupling) {
          kind =
              IBCoupling::describe(candidate.get(), runtime_.size(), info).kind;
          if (kind == IB_MODEL_VIEW_INVALID || !IBCoupling::valid_info(info))
            result.status = InsertionStatus::InvalidInput;
        }
        if (result.accepted() ||
            result.status == InsertionStatus::ClearanceViolation) {
          if (distributed) {
            // Exact comparison, not a hash: detect divergent setup even when both
            // ranks happen to report the same acceptance decision.
            auto append = [&](const auto &value) {
              const char *bytes = reinterpret_cast<const char *>(&value);
              signature.insert(signature.end(), bytes, bytes + sizeof(value));
            };
            append(type);
            append(owner);
            append(clearance);
            append(kind);
            append(info.sparse);
            append(info.level);
            append(info.node_count);
            for (size_t d = 0; d < 3; ++d) {
              append(config_.origin[d]);
              append(config_.dimensions[d]);
              append(geometry_upper_[d]);
            }
            append(config_.cell_width);
            append(runtime_.size());
            for (size_t i = 0; i < runtime_.size(); ++i) {
              append(runtime_.entry(i)->type);
              append(runtime_.entry(i)->owner_pid);
            }
            append(snapshots_.size());
            for (const auto &snapshot : snapshots_) {
              append(snapshot.id);
              append(snapshot.points.size());
              for (const auto &point : snapshot.points)
                for (auto x : point)
                  append(x);
            }
            if (coupling) {
              append(coupling->registrations_.size());
              for (const auto &registration : coupling->registrations_) {
                append(registration.id);
                append(registration.kind);
                append(registration.info.sparse);
                append(registration.info.level);
                append(registration.info.node_count);
              }
            }
            append(geometry.count);
            for (size_t i = 0; i < geometry.count; ++i)
              for (auto x : geometry.points[i])
                append(x);
            ModelState state;
            candidate->pack_state(state);
            append(state.ints.size());
            for (auto x : state.ints)
              append(x);
            append(state.reals.size());
            for (auto x : state.reals)
              append(x);
            append(state.bytes.size());
            for (auto x : state.bytes)
              append(x);
          }
          if (result.accepted()) {
            prepared = {runtime_.size(),
                        {geometry.points, geometry.points + geometry.count}};
            snapshots_.reserve(snapshots_.size() + 1);
            runtime_.reserve_model();
            if (coupling)
              coupling->registrations_.reserve(coupling->registrations_.size() +
                                               1);
          }
        }
      }
    }
  } catch (...) {
    result = failure(InsertionStatus::PreparationFailed);
  }
  auto status = agree_status(result.status);
  if (status != InsertionStatus::Accepted &&
      status != InsertionStatus::ClearanceViolation)
    return failure(status);
#ifdef ELFF_USE_MPI
  if (distributed) {
    unsigned long long size = signature.size(), root_size = size;
    MPI_Bcast(&root_size, 1, MPI_UNSIGNED_LONG_LONG, 0,
              runtime_.communicator());
    int mismatch = size != root_size;
    // Fixed stack storage avoids allocation failures between collectives.
    char buffer[65536];
    for (unsigned long long offset = 0; offset < root_size;
         offset += sizeof(buffer)) {
      const int count = static_cast<int>(
          std::min<unsigned long long>(sizeof(buffer), root_size - offset));
      if (runtime_.rank() == 0)
        std::memcpy(buffer, signature.data() + offset, count);
      MPI_Bcast(buffer, count, MPI_BYTE, 0, runtime_.communicator());
      if (size == root_size &&
          std::memcmp(buffer, signature.data() + offset, count) != 0)
        mismatch = 1;
    }
    MPI_Allreduce(MPI_IN_PLACE, &mismatch, 1, MPI_INT, MPI_MAX,
                  runtime_.communicator());
    if (mismatch)
      return failure(InsertionStatus::CollectiveMismatch);
  }
#endif
  if (!result.accepted())
    return result;
  const ModelId id = prepared.id;
  info.owner_rank = owner;
  if (runtime_.rank() == owner) {
    if (coupling) {
      auto *ib = dynamic_cast<IBModel *>(candidate.get());
      candidate.release();
      auto *adopted =
          coupling->adopt_owned_model(type, std::unique_ptr<IBModel>(ib), info);
      ELFF_ASSERT(adopted, "Prepared collective IB adoption failed.");
    } else {
      auto *adopted = runtime_.adopt_owned_model(type, std::move(candidate));
      ELFF_ASSERT(adopted, "Prepared collective model adoption failed.");
    }
  } else {
    const auto adopted = coupling ? coupling->add_remote_model(type, kind, info)
                                  : runtime_.add_remote_model(type, owner);
    ELFF_ASSERT(adopted == id, "Prepared collective remote adoption failed.");
  }
  snapshots_.push_back(std::move(prepared));
  result.model_id = id;
  return result;
}

int ModelInsertion::seal() {
  sealed_ = true;
  return 0;
}

} // namespace ELFF::Models
