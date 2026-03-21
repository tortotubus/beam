#pragma once

#include <vector>
#include <cstdint>
#include "config/config.hpp"

namespace ELFF {
namespace Models {

struct IBModelState {
    std::vector<int64_t> ints;
    std::vector<real_t> reals;
    std::vector<char> bytes;
};


class IBModel
{
public:
  virtual ~IBModel() = default;
  virtual void pack_state(IBModelState& s) const = 0;
  virtual void unpack_state(const IBModelState& s) = 0;
};

}
}