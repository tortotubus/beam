#pragma once

#include <vector>
#include <cstdint>
#include "elff/config/config.hpp"
#include "elff/io/CXX/vtkPolyData.hpp"

namespace ELFF {
namespace Models {

/**
 * @brief Serialized storage container for immersed-boundary model state.
 *
 * Implementations of @ref IBModel use this structure to pack any persistent
 * integer, floating-point, and raw-byte data needed to restore a model
 * instance later.
 */
struct IBModelState {
    std::vector<int64_t> ints;  ///< Packed integer-valued state entries.
    std::vector<real_t> reals;  ///< Packed floating-point state entries.
    std::vector<char> bytes;    ///< Packed opaque byte data for custom state.
};

/**
 * @brief Abstract serialization interface for immersed-boundary models.
 *
 * Concrete IB model implementations expose their restartable state through
 * @ref pack_state and @ref unpack_state so callers can checkpoint and restore
 * model instances without depending on solver-specific data layouts.
 */
class IBModel
{
public:
  virtual ~IBModel() = default;

  /**
   * @brief Serialize the current model state into a caller-provided container.
   *
   * @param s Destination state container to populate.
   */
  virtual void pack_state(IBModelState& s) const = 0;

  /**
   * @brief Restore the model state from a previously packed container.
   *
   * @param s Source state container holding serialized model data.
   */
  virtual void unpack_state(const IBModelState& s) = 0;

  /**
   * @brief Optionally append PolyData datasets exported by this model.
   *
   * Models without a useful mesh/topology export keep the default no-op.
   */
  virtual void append_vtk_polydata(
    std::vector<IO::CXX::vtkPolyData>& datasets) const
  {
    (void)datasets;
  }
};

}
}
