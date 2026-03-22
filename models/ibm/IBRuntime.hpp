#pragma once

#include "general/error.hpp"
#include "mpi/CommHandle.hpp"
#include "models/ibm/IBModel.hpp"

#include <cstdint>
#include <cstdio>
#include <vector>

namespace ELFF {
namespace Models {

struct IBRuntimeState
{
  std::vector<int> ids;
  std::vector<IBModelState> states;
};

class IBRuntime
{
public:
  struct Entry
  {
    int id;
    int pid;
    IBModel* model;
  };

  IBRuntime()
    : next_id(0)
  {}

  int register_model(IBModel& model)
  {
    for (const auto& entry : entries) {
      if (entry.model == &model)
        return entry.id;
    }

    const int id = next_id++;
    entries.push_back({id, -1, &model});
    return id;
  }

  int deregister_model(IBModel& model)
  {
    for (auto it = entries.begin(); it != entries.end(); ++it) {
      if (it->model == &model) {
        entries.erase(it);
        return 0;
      }
    }

    return -1;
  }

  int set_model_pid(IBModel& model, int pid)
  {
    for (auto& entry : entries) {
      if (entry.model == &model) {
        entry.pid = pid;
        return 0;
      }
    }

    return -1;
  }

#ifdef ELFF_USE_MPI
  void set_communicator(MPI_Comm comm)
  {
    communicator =
      (comm != MPI_COMM_NULL) ? ELFF::MPI::CommHandle::duplicate(comm)
                              : ELFF::MPI::CommHandle();
  }
#endif

  size_t size() const { return entries.size(); }

  void export_state(IBRuntimeState& state) const
  {
    state.ids.clear();
    state.states.clear();

    state.ids.reserve(entries.size());
    state.states.resize(entries.size());

    for (size_t i = 0; i < entries.size(); ++i) {
      state.ids.push_back(entries[i].id);
      entries[i].model->pack_state(state.states[i]);
    }
  }

  void import_state(const IBRuntimeState& state)
  {
    ELFF_ASSERT(state.ids.size() == entries.size(),
                "IBRuntime::import_state(): model count mismatch.\n");
    ELFF_ASSERT(state.states.size() == entries.size(),
                "IBRuntime::import_state(): state count mismatch.\n");

    for (size_t i = 0; i < state.ids.size(); ++i) {
      Entry* entry = find_entry(state.ids[i]);
      ELFF_ASSERT(entry,
                  "IBRuntime::import_state(): unknown model id in state.\n");
      entry->model->unpack_state(state.states[i]);
    }
  }

  int write_checkpoint(const char* fname) const
  {
    if (!should_checkpoint_locally())
      return 0;

    FILE* fp = std::fopen(fname, "wb");
    if (!fp)
      return -1;

    const int rc = write_checkpoint(fp);
    std::fclose(fp);
    return rc;
  }

  int read_checkpoint(const char* fname)
  {
    if (!should_checkpoint_locally())
      return 0;

    FILE* fp = std::fopen(fname, "rb");
    if (!fp)
      return -1;

    const int rc = read_checkpoint(fp);
    std::fclose(fp);
    return rc;
  }

  int write_checkpoint(FILE* fp) const
  {
    if (!should_checkpoint_locally())
      return 0;

    if (!fp)
      return -1;

    IBRuntimeState state;
    export_state(state);

    FileHeader header = {
      .magic = file_magic,
      .version = file_version,
      .nmodels = static_cast<int64_t>(state.ids.size())
    };

    if (std::fwrite(&header, sizeof(header), 1, fp) != 1)
      return -1;

    for (size_t i = 0; i < state.ids.size(); ++i) {
      const int32_t id = state.ids[i];
      const int64_t ni = static_cast<int64_t>(state.states[i].ints.size());
      const int64_t nr = static_cast<int64_t>(state.states[i].reals.size());
      const int64_t nb = static_cast<int64_t>(state.states[i].bytes.size());

      if (std::fwrite(&id, sizeof(id), 1, fp) != 1)
        return -1;
      if (std::fwrite(&ni, sizeof(ni), 1, fp) != 1)
        return -1;
      if (std::fwrite(&nr, sizeof(nr), 1, fp) != 1)
        return -1;
      if (std::fwrite(&nb, sizeof(nb), 1, fp) != 1)
        return -1;

      if (ni > 0 &&
          std::fwrite(state.states[i].ints.data(), sizeof(int64_t), ni, fp) !=
            static_cast<size_t>(ni))
        return -1;
      if (nr > 0 &&
          std::fwrite(state.states[i].reals.data(), sizeof(real_t), nr, fp) !=
            static_cast<size_t>(nr))
        return -1;
      if (nb > 0 &&
          std::fwrite(state.states[i].bytes.data(), sizeof(char), nb, fp) !=
            static_cast<size_t>(nb))
        return -1;
    }

    return 0;
  }

  int read_checkpoint(FILE* fp)
  {
    if (!should_checkpoint_locally())
      return 0;

    if (!fp)
      return -1;

    FileHeader header = {0};
    if (std::fread(&header, sizeof(header), 1, fp) != 1)
      return -1;

    if (header.magic != file_magic || header.version != file_version)
      return -1;

    IBRuntimeState state;
    state.ids.resize(static_cast<size_t>(header.nmodels));
    state.states.resize(static_cast<size_t>(header.nmodels));

    for (size_t i = 0; i < state.ids.size(); ++i) {
      int32_t id = -1;
      int64_t ni = 0, nr = 0, nb = 0;

      if (std::fread(&id, sizeof(id), 1, fp) != 1)
        return -1;
      if (std::fread(&ni, sizeof(ni), 1, fp) != 1)
        return -1;
      if (std::fread(&nr, sizeof(nr), 1, fp) != 1)
        return -1;
      if (std::fread(&nb, sizeof(nb), 1, fp) != 1)
        return -1;

      if (ni < 0 || nr < 0 || nb < 0)
        return -1;

      state.ids[i] = id;
      state.states[i].ints.resize(static_cast<size_t>(ni));
      state.states[i].reals.resize(static_cast<size_t>(nr));
      state.states[i].bytes.resize(static_cast<size_t>(nb));

      if (ni > 0 &&
          std::fread(state.states[i].ints.data(), sizeof(int64_t), ni, fp) !=
            static_cast<size_t>(ni))
        return -1;
      if (nr > 0 &&
          std::fread(state.states[i].reals.data(), sizeof(real_t), nr, fp) !=
            static_cast<size_t>(nr))
        return -1;
      if (nb > 0 &&
          std::fread(state.states[i].bytes.data(), sizeof(char), nb, fp) !=
            static_cast<size_t>(nb))
        return -1;
    }

    import_state(state);
    return 0;
  }

private:
  struct FileHeader
  {
    uint64_t magic;
    int32_t version;
    int64_t nmodels;
  };

  static constexpr uint64_t file_magic = 0x494252554E54494DULL; // "IBRUNTIM"
  static constexpr int32_t file_version = 1;

  Entry* find_entry(int id)
  {
    for (auto& entry : entries) {
      if (entry.id == id)
        return &entry;
    }
    return nullptr;
  }

  Entry* find_entry(IBModel& model)
  {
    for (auto& entry : entries) {
      if (entry.model == &model)
        return &entry;
    }
    return nullptr;
  }

  int checkpoint_owner_pid() const
  {
    int owner_pid = -1;

    for (const auto& entry : entries) {
      const int entry_pid = entry.pid >= 0 ? entry.pid : 0;
      if (owner_pid < 0) {
        owner_pid = entry_pid;
      } else {
        ELFF_ASSERT(owner_pid == entry_pid,
                    "IBRuntime::checkpoint_owner_pid(): mixed model owner "
                    "pids are not supported for checkpointing.\n");
      }
    }

    return owner_pid >= 0 ? owner_pid : 0;
  }

  int local_rank() const
  {
#ifdef ELFF_USE_MPI
    if (communicator) {
      int rank = 0;
      ELFF_ASSERT(MPI_Comm_rank(communicator.get(), &rank) == MPI_SUCCESS,
                  "IBRuntime::local_rank(): MPI_Comm_rank() failed.\n");
      return rank;
    }
#endif
    return 0;
  }

  bool should_checkpoint_locally() const
  {
    return local_rank() == checkpoint_owner_pid();
  }

private:
  std::vector<Entry> entries;
  int next_id;
  ELFF::MPI::CommHandle communicator;
};

} // namespace Models
} // namespace ELFF
