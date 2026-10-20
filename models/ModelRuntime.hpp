#pragma once

#include "elff/io/CXX/vtkHDFPolyData.hpp"
#include "elff/io/CXX/vtkPolyData.hpp"
#include "elff/models/Model.hpp"
#include "elff/models/RuntimeComponent.hpp"
#include "elff/general/error.hpp"
#include <typeindex>
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
enum class ModelType : int32_t
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

struct RuntimeModelState
{
  ModelType type = ModelType::Invalid;
  ModelState state;
};

/** Serializable, ordered snapshot of every model owned by a runtime. */
struct RuntimeState
{
  std::vector<RuntimeModelState> models;
};

using ModelId = size_t;
inline constexpr ModelId invalid_model_id = std::numeric_limits<ModelId>::max();

/** @brief Owns ordered models and optional components independently of coupling. */
class ModelRuntime
{
public:
  struct Entry {
    ModelType type;
    int owner_pid;
    std::unique_ptr<Model> model;
  };

  ModelRuntime() = default;
  virtual ~ModelRuntime() { components_.clear(); entries_.clear(); }
  ModelRuntime(const ModelRuntime&) = delete;
  ModelRuntime& operator=(const ModelRuntime&) = delete;

  /** @brief Attach one component of each exact type before sealing. */
  template<class T> T* attach_component() {
    static_assert(std::is_base_of_v<RuntimeComponent, T>);
    if (auto* existing = component<T>()) return existing;
    if (sealed_ || sealing_) return nullptr;
    auto owned = std::make_unique<T>(*this);
    T* result = owned.get();
    components_.emplace_back(std::type_index(typeid(T)), std::move(owned));
    return result;
  }
  template<class T> T* component() const {
    for (const auto& entry : components_)
      if (entry.first == std::type_index(typeid(T)))
        return static_cast<T*>(entry.second.get());
    return nullptr;
  }

  /** @brief Adopt a locally owned model. Its ID is its immutable insertion index. */
  template<class T>
  T* adopt_owned_model(ModelType type, std::unique_ptr<T> model) {
    static_assert(std::is_base_of_v<Model, T>);
    if (sealed_ || sealing_ || !model || type == ModelType::Invalid) return nullptr;
    T* result = model.get();
    entries_.push_back({type, rank(), std::move(model)});
    return result;
  }
  /** @brief Append a remote placeholder, or return invalid_model_id on invalid input. */
  ModelId add_remote_model(ModelType type, int owner) {
    if (sealed_ || sealing_ || type == ModelType::Invalid || owner < 0 ||
        owner >= ranks() || owner == rank()) return invalid_model_id;
    const ModelId id = entries_.size();
    entries_.push_back({type, owner, nullptr});
    return id;
  }
  const Entry* entry(ModelId id) const {
    return id < entries_.size() ? &entries_[id] : nullptr;
  }
  Model* model(ModelId id) const {
    const auto* e = entry(id);
    return e ? e->model.get() : nullptr;
  }
  /** @brief Reserve before a multi-part registration so its commit cannot allocate. */
  bool reserve_model() {
    if (sealed_ || sealing_) return false;
    entries_.reserve(entries_.size() + 1);
    return true;
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

  /** @brief Freeze models and components; repeated successful calls are harmless. */
  int seal() {
    if (sealed_) return 0;
    if (sealing_) return -1;
    sealing_ = true;
    try {
      for (auto& component : components_)
        if (component.second->seal() != 0) { sealing_ = false; return -1; }
    } catch (...) { sealing_ = false; throw; }
    sealing_ = false;
    sealed_ = true;
    return 0;
  }
  bool mutable_registration() const { return !sealed_ && !sealing_; }

  bool sealed() const { return sealed_; }
  size_t size() const { return entries_.size(); }

#ifdef ELFF_USE_MPI
  int set_communicator(MPI_Comm comm)
  {
    if (!mutable_registration() || !entries_.empty()) return -1;
    communicator_ =
      (comm != MPI_COMM_NULL) ? ELFF::MPI::CommHandle::duplicate(comm)
                              : ELFF::MPI::CommHandle();
    return 0;
  }
#endif

  int export_state(RuntimeState& state)
  {
    if (seal() != 0) return -1;
    const bool distributed = checkpoint_owner_pid() < 0;
    if (!distributed)
      for (const auto& entry : entries_)
        if (!entry.model) return -1;
    state.models.clear();
    state.models.resize(entries_.size());
    for (size_t i = 0; i < entries_.size(); ++i) {
      state.models[i].type = entries_[i].type;
      auto& packed = state.models[i].state;
      if (entries_[i].model)
        entries_[i].model->pack_state(packed);
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

  int import_state(const RuntimeState& state)
  {
    if (seal() != 0) return -1;
    if (state.models.size() != entries_.size())
      return -1;

    // Validate the complete ordered type sequence before mutating any model.
    for (size_t i = 0; i < entries_.size(); ++i) {
      if (state.models[i].type != entries_[i].type)
        return -1;
    }

    for (size_t i = 0; i < entries_.size(); ++i)
      if (entries_[i].model)
        entries_[i].model->unpack_state(state.models[i].state);
    return 0;
  }

  IO::CXX::vtkPolyData to_vtk_polydata()
  {
    if (seal() != 0) ELFF_ABORT("Unable to seal model runtime for output.");
    std::vector<IO::CXX::vtkPolyData> datasets;
    for (const auto& entry : entries_)
      if (entry.model) entry.model->append_vtk_polydata(datasets);
    return IO::CXX::vtkPolyData::append_many(datasets);
  }

  int write_polydata(const char* fname, bool overwrite = true)
  {
    if (seal() != 0) return -1;
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
    if (seal() != 0) return -1;
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
    if (seal() != 0) return -1;
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
    if (seal() != 0) return -1;
    const bool distributed = checkpoint_owner_pid() < 0;
    if (!distributed && !should_checkpoint_locally()) return 0;
    RuntimeState state;
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
    if (seal() != 0) return -1;
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
      ELFF_LOG("ModelRuntime checkpoint version " << header.version
               << " is unsupported; expected " << file_version << ".");
      return -1;
    }
    if (header.nmodels < 0)
      return -1;
    if (static_cast<uint64_t>(header.nmodels) >
        static_cast<uint64_t>(std::numeric_limits<size_t>::max()))
      return -1;

    RuntimeState state;
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

      model_state.type = static_cast<ModelType>(type);
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
  std::vector<std::pair<std::type_index, std::unique_ptr<RuntimeComponent>>> components_;
  bool sealing_ = false;
  bool sealed_ = false;
  ELFF::MPI::CommHandle communicator_;
};

} // namespace Models
} // namespace ELFF
