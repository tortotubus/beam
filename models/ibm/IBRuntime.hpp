#pragma once

// #include "IBForceCoupled.hpp"
// #include "IBVelocityCoupled.hpp"

#include <vector>

namespace ELFF {
namespace Models {

struct IBModelState
{
  std::vector<int> iarr;
  std::vector<real_t> rarr;
};

class IBModel
{
public:
  virtual void run(real_t dt) = 0;
  virtual void export_state(IBModelState* s) const = 0;
  virtual void import_state(const IBModelState* s) = 0;
  virtual ~IBModel() = default;
};

class IBRuntime
{
public:
  IBRuntime();

  int register_model(IBModel &model) {

  }

  void run(real_t dt) {

  }

private:
  std::vector<IBModel*> models;
  std::vector<int> global_ids;
};

}
}