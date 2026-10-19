#pragma once

#include "elff/c/models/ibm/IBRuntime.h"
#include "elff/io/CXX/vtkHDFPolyData.hpp"
#include "elff/io/CXX/vtkPolyData.hpp"
#include "elff/models/ibm/IBForceCoupled.hpp"
#include "elff/models/ibm/IBModel.hpp"
#include "elff/models/ibm/IBVelocityCoupled.hpp"
#include "elff/mpi/CommHandle.hpp"

#include <cstdint>
#include <climits>
#include <algorithm>
#include <cstdio>
#include <limits>
#include <memory>
#include <string>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace ELFF {
namespace Models {

/** Exact concrete model type stored in checkpoints. Values are persistent. */
enum class IBRuntimeModelType : int32_t
{
  Invalid = 0,
  EulerBeam = 1,
  EulerBeamPenalty = 2,
  EulerBeamGGL = 3,
  EulerBeamADDM = 4,
  EulerBeamHuang = 5,
  Capsule = 6,
  RigidBodyCircle = 7,
  RigidBodyFibre = 8,
  RigidBodySphere = 9,
  PinnedRigidBodyCircle = 10,
  PinnedRigidBodyFibre = 11,
  PinnedRigidBodySphere = 12,
  SpringCircle = 13,
  BendingRodImplicit = 14
};

struct IBRuntimeModelState
{
  IBRuntimeModelType type = IBRuntimeModelType::Invalid;
  IBModelState state;
};

/** Serializable, ordered snapshot of every model owned by a runtime. */
struct IBRuntimeState
{
  std::vector<IBRuntimeModelState> models;
};

/** Owns an ordered set of immersed-boundary models and their C views. */
class IBRuntime
{
public:
  struct Entry
  {
    IBRuntimeModelType type;
    int owner_pid;
    ib_model_view_kind_t coupling_kind;
    std::unique_ptr<IBModel> model;
    IBModel* base_model;
    void* coupling_handle;
    ib_model_info_t info = {0, 0, -1, 0};
  };

  IBRuntime() = default;

  template<typename T>
  T* adopt_owned_model(IBRuntimeModelType type, std::unique_ptr<T> model)
  {
    static_assert(std::is_base_of<IBModel, T>::value,
                  "runtime models must derive from IBModel");

    if (sealed_ || !model || type == IBRuntimeModelType::Invalid)
      return nullptr;

    T* result = model.get();
    IBModel* base = result;
    IBForceCoupled* force = dynamic_cast<IBForceCoupled*>(base);
    IBVelocityCoupled* velocity = dynamic_cast<IBVelocityCoupled*>(base);
    if ((force == nullptr) == (velocity == nullptr))
      return nullptr;

    const ib_model_view_kind_t kind =
      force ? IB_MODEL_VIEW_FORCE_COUPLED
            : IB_MODEL_VIEW_VELOCITY_COUPLED;
    void* coupling_handle =
      force ? static_cast<void*>(force) : static_cast<void*>(velocity);
    std::unique_ptr<IBModel> owner = std::move(model);
    entries_.push_back(
      {type, 0, kind, std::move(owner), base, coupling_handle});
    return result;
  }

  int rank() const { return local_rank(); }
  int ranks() const {
    int count = 1;
#ifdef ELFF_USE_MPI
    if (communicator_)
      MPI_Comm_size(communicator_.get(), &count);
#endif
    return count;
  }

  void add_remote_force_model(IBRuntimeModelType type, ib_model_info_t info) {
    if (sealed_) throw std::logic_error("Runtime is sealed");
    entries_.push_back({type, info.owner_rank,
      IB_MODEL_VIEW_FORCE_COUPLED, nullptr, nullptr, nullptr, info});
  }
  void add_remote_huang(ib_model_info_t info) {
    add_remote_force_model(IBRuntimeModelType::EulerBeamHuang, info);
  }
  void configure_last_model(ib_model_info_t info) {
    entries_.back().info = info;
    entries_.back().owner_pid = info.owner_rank;
  }

  int seal()
  {
    if (sealed_)
      return 0;

    views_.clear();
    views_.reserve(entries_.size());
    for (const auto& entry : entries_)
      views_.push_back({entry.coupling_kind, entry.coupling_handle, &entry.info});
    sealed_ = true;
    return 0;
  }

  int get_model_views(ib_model_view_array_t& views)
  {
    seal();
    views.count = views_.size();
    views.data = views_.empty() ? nullptr : views_.data();
    return 0;
  }

  bool sealed() const { return sealed_; }
  size_t size() const { return entries_.size(); }

#ifdef ELFF_USE_MPI
  void set_communicator(MPI_Comm comm)
  {
    if (!entries_.empty())
      throw std::logic_error("Set communicator before model construction");
    communicator_ =
      (comm != MPI_COMM_NULL) ? ELFF::MPI::CommHandle::duplicate(comm)
                              : ELFF::MPI::CommHandle();
  }
#endif

  int export_state(IBRuntimeState& state)
  {
    seal();
    const bool distributed = checkpoint_owner_pid() < 0;
    if (!distributed)
      for (const auto& entry : entries_)
        if (!entry.base_model) return -1;
    state.models.clear();
    state.models.resize(entries_.size());
    for (size_t i = 0; i < entries_.size(); ++i) {
      state.models[i].type = entries_[i].type;
      auto& packed = state.models[i].state;
      if (entries_[i].base_model)
        entries_[i].base_model->pack_state(packed);
#ifdef ELFF_USE_MPI
      if (distributed) {
        broadcast_vector(packed.ints, entries_[i].owner_pid);
        broadcast_vector(packed.reals, entries_[i].owner_pid);
        broadcast_vector(packed.bytes, entries_[i].owner_pid);
      }
#endif
    }
    return 0;
  }

  int import_state(const IBRuntimeState& state)
  {
    seal();
    if (state.models.size() != entries_.size())
      return -1;

    // Validate the complete ordered type sequence before mutating any model.
    for (size_t i = 0; i < entries_.size(); ++i) {
      if (state.models[i].type != entries_[i].type)
        return -1;
    }

    for (size_t i = 0; i < entries_.size(); ++i)
      if (entries_[i].base_model)
        entries_[i].base_model->unpack_state(state.models[i].state);
    return 0;
  }

  IO::CXX::vtkPolyData to_vtk_polydata()
  {
    seal();
    std::vector<IO::CXX::vtkPolyData> datasets;
    for (const auto& entry : entries_)
      if (entry.base_model) entry.base_model->append_vtk_polydata(datasets);
    return IO::CXX::vtkPolyData::append_many(datasets);
  }

  int write_polydata(const char* fname, bool overwrite = true)
  {
    seal();
    if (checkpoint_owner_pid() < 0) return -1;
    if (!should_checkpoint_locally())
      return 0;
    if (!fname)
      return -1;

    IO::CXX::vtkPolyData pd = to_vtk_polydata();
    IO::CXX::vtkHDFPolyData writer(std::string(fname), pd);
    writer.write_new_static(overwrite);
    return 0;
  }

  int write_checkpoint(const char* fname)
  {
    seal();
#ifdef ELFF_USE_MPI
    if (checkpoint_owner_pid() < 0) {
      // All owners contribute; rank zero publishes the existing ordered format.
      FILE* fp = local_rank() == 0 && fname ? std::fopen(fname, "wb") : nullptr;
      int rc = local_rank() == 0 && !fp ? -1 : 0;
      MPI_Bcast(&rc, 1, MPI_INT, 0, communicator_.get());
      if (rc != 0) return rc;
      rc = write_checkpoint(fp);
      if (fp && std::fclose(fp) != 0) rc = -1;
      MPI_Bcast(&rc, 1, MPI_INT, 0, communicator_.get());
      return rc;
    }
#endif
    if (!should_checkpoint_locally())
      return 0;
    if (!fname)
      return -1;

    FILE* fp = std::fopen(fname, "wb");
    if (!fp)
      return -1;
    const int rc = write_checkpoint(fp);
    if (std::fclose(fp) != 0)
      return -1;
    return rc;
  }

  int read_checkpoint(const char* fname)
  {
    seal();
    if (checkpoint_owner_pid() >= 0 && !should_checkpoint_locally())
      return 0;
    if (!fname)
      return -1;

    FILE* fp = std::fopen(fname, "rb");
    if (!fp)
      return -1;
    const int rc = read_checkpoint(fp);
    std::fclose(fp);
    return rc;
  }

  int write_checkpoint(FILE* fp)
  {
    seal();
    const bool distributed = checkpoint_owner_pid() < 0;
    if (!distributed && !should_checkpoint_locally()) return 0;
    IBRuntimeState state;
    if (export_state(state) != 0) return -1;
    if (distributed && local_rank() != 0) return 0;
    if (!fp) return -1;
    const FileHeader header = {
      file_magic, file_version, static_cast<int64_t>(state.models.size())};
    if (std::fwrite(&header, sizeof(header), 1, fp) != 1)
      return -1;

    for (const auto& model_state : state.models) {
      const int32_t type = static_cast<int32_t>(model_state.type);
      const int64_t ni = static_cast<int64_t>(model_state.state.ints.size());
      const int64_t nr = static_cast<int64_t>(model_state.state.reals.size());
      const int64_t nb = static_cast<int64_t>(model_state.state.bytes.size());

      if (std::fwrite(&type, sizeof(type), 1, fp) != 1 ||
          std::fwrite(&ni, sizeof(ni), 1, fp) != 1 ||
          std::fwrite(&nr, sizeof(nr), 1, fp) != 1 ||
          std::fwrite(&nb, sizeof(nb), 1, fp) != 1)
        return -1;
      if (ni > 0 &&
          std::fwrite(model_state.state.ints.data(), sizeof(int64_t), ni, fp) !=
            static_cast<size_t>(ni))
        return -1;
      if (nr > 0 &&
          std::fwrite(model_state.state.reals.data(), sizeof(real_t), nr, fp) !=
            static_cast<size_t>(nr))
        return -1;
      if (nb > 0 &&
          std::fwrite(model_state.state.bytes.data(), sizeof(char), nb, fp) !=
            static_cast<size_t>(nb))
        return -1;
    }
    return 0;
  }

  int read_checkpoint(FILE* fp)
  {
    seal();
    if (checkpoint_owner_pid() >= 0 && !should_checkpoint_locally())
      return 0;
    if (!fp)
      return -1;

    FileHeader header = {};
    if (std::fread(&header, sizeof(header), 1, fp) != 1)
      return -1;
    if (header.magic != file_magic)
      return -1;
    if (header.version != file_version) {
      std::fprintf(stderr,
                   "IBRuntime checkpoint version %d is unsupported; "
                   "expected version %d.\n",
                   header.version,
                   file_version);
      return -1;
    }
    if (header.nmodels < 0)
      return -1;
    if (static_cast<uint64_t>(header.nmodels) >
        static_cast<uint64_t>(std::numeric_limits<size_t>::max()))
      return -1;

    IBRuntimeState state;
    state.models.resize(static_cast<size_t>(header.nmodels));
    for (auto& model_state : state.models) {
      int32_t type = 0;
      int64_t ni = 0, nr = 0, nb = 0;
      if (std::fread(&type, sizeof(type), 1, fp) != 1 ||
          std::fread(&ni, sizeof(ni), 1, fp) != 1 ||
          std::fread(&nr, sizeof(nr), 1, fp) != 1 ||
          std::fread(&nb, sizeof(nb), 1, fp) != 1 || ni < 0 || nr < 0 ||
          nb < 0)
        return -1;
      if (static_cast<uint64_t>(ni) >
            static_cast<uint64_t>(std::numeric_limits<size_t>::max()) ||
          static_cast<uint64_t>(nr) >
            static_cast<uint64_t>(std::numeric_limits<size_t>::max()) ||
          static_cast<uint64_t>(nb) >
            static_cast<uint64_t>(std::numeric_limits<size_t>::max()))
        return -1;

      model_state.type = static_cast<IBRuntimeModelType>(type);
      model_state.state.ints.resize(static_cast<size_t>(ni));
      model_state.state.reals.resize(static_cast<size_t>(nr));
      model_state.state.bytes.resize(static_cast<size_t>(nb));
      if (ni > 0 &&
          std::fread(model_state.state.ints.data(), sizeof(int64_t), ni, fp) !=
            static_cast<size_t>(ni))
        return -1;
      if (nr > 0 &&
          std::fread(model_state.state.reals.data(), sizeof(real_t), nr, fp) !=
            static_cast<size_t>(nr))
        return -1;
      if (nb > 0 &&
          std::fread(model_state.state.bytes.data(), sizeof(char), nb, fp) !=
            static_cast<size_t>(nb))
        return -1;
    }
    return import_state(state);
  }

private:
#ifdef ELFF_USE_MPI
  template<class T>
  void broadcast_vector(std::vector<T>& values, int owner)
  {
    uint64_t size = values.size();
    MPI_Bcast(&size, 1, MPI_UINT64_T, owner, communicator_.get());
    values.resize(static_cast<size_t>(size));
    size_t bytes = values.size()*sizeof(T);
    char* data = reinterpret_cast<char*>(values.data());
    while (bytes) {
      int count = static_cast<int>(std::min(bytes, static_cast<size_t>(INT_MAX)));
      MPI_Bcast(data, count, MPI_BYTE, owner, communicator_.get());
      data += count;
      bytes -= count;
    }
  }
#endif

  struct FileHeader
  {
    uint64_t magic;
    int32_t version;
    int64_t nmodels;
  };

  static constexpr uint64_t file_magic = 0x494252554E54494DULL; // IBRUNTIM
  static constexpr int32_t file_version = 2;

  int local_rank() const
  {
#ifdef ELFF_USE_MPI
    if (communicator_) {
      int rank = 0;
      if (MPI_Comm_rank(communicator_.get(), &rank) != MPI_SUCCESS)
        return -1;
      return rank;
    }
#endif
    return 0;
  }

  int checkpoint_owner_pid() const
  {
    // Negative means model states must be collected from multiple owners.
    const int owner = entries_.empty() ? 0 : entries_.front().owner_pid;
    for (const auto& entry : entries_) {
      if (entry.owner_pid != owner)
        return -1;
    }
    return owner;
  }

  bool should_checkpoint_locally() const
  {
    return local_rank() == checkpoint_owner_pid();
  }

  std::vector<Entry> entries_;
  std::vector<ib_model_view_t> views_;
  bool sealed_ = false;
  ELFF::MPI::CommHandle communicator_;
};

} // namespace Models
} // namespace ELFF
