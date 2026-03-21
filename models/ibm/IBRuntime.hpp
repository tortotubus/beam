#pragma once

#include "general/error.hpp"
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
    entries.push_back({id, &model});
    return id;
  }

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
    FILE* fp = std::fopen(fname, "wb");
    if (!fp)
      return -1;

    const int rc = write_checkpoint(fp);
    std::fclose(fp);
    return rc;
  }

  int read_checkpoint(const char* fname)
  {
    FILE* fp = std::fopen(fname, "rb");
    if (!fp)
      return -1;

    const int rc = read_checkpoint(fp);
    std::fclose(fp);
    return rc;
  }

  int write_checkpoint(FILE* fp) const
  {
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

private:
  std::vector<Entry> entries;
  int next_id;
};

} // namespace Models
} // namespace ELFF
