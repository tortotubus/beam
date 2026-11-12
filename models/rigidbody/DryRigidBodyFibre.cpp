#include "elff/models/rigidbody/DryRigidBodyFibre.hpp"
#include <cmath>
#include <limits>

namespace ELFF::Models {
DryRigidBodyFibre::DryRigidBodyFibre(real_t length, real_t diameter, size_t nodes,
                                     real_t linear_density, const Vec3& center,
                                     const Vec3& velocity, const Quat& orientation,
                                     const Vec3& angular_velocity)
    : RigidBodyFibre(length, diameter, nodes, linear_density) {
  ELFF_VERIFY(std::isfinite(length) && std::isfinite(diameter) && std::isfinite(linear_density) &&
                  std::isfinite(mass()),
              "Fibre parameters must be finite.");
  set_initial_state(center, velocity, orientation, angular_velocity);
}

void DryRigidBodyFibre::pack_state(ModelState& s) const {
  ELFF_VERIFY(nodes() <= static_cast<size_t>(std::numeric_limits<int64_t>::max()) &&
                  time_iteration() <= static_cast<size_t>(std::numeric_limits<int64_t>::max()),
              "Fibre checkpoint integer overflow.");
  s = {};
  s.ints = {1, static_cast<int64_t>(nodes()), static_cast<int64_t>(time_iteration())};
  s.reals = {length(), diameter(), linear_density(), time()};
  for (const auto& v : {center_of_mass(), com_velocity()})
    s.reals.insert(s.reals.end(), v.begin(), v.end());
  s.reals.insert(s.reals.end(), pose().begin(), pose().end());
  s.reals.insert(s.reals.end(), angular_velocity_body().begin(), angular_velocity_body().end());
}

void DryRigidBodyFibre::unpack_state(const ModelState& s) {
  ELFF_VERIFY(s.ints.size() == 3 && s.ints[0] == 1 && s.ints[1] >= 2 &&
                  static_cast<uint64_t>(s.ints[1]) == nodes() && s.ints[2] >= 0 &&
                  static_cast<uint64_t>(s.ints[2]) <= std::numeric_limits<size_t>::max() &&
                  s.reals.size() == 17 && s.bytes.empty(),
              "Incompatible dry fibre state.");
  for (auto value : s.reals)
    ELFF_VERIFY(std::isfinite(value), "Non-finite dry fibre state.");
  ELFF_VERIFY(s.reals[0] == length() && s.reals[1] == diameter() &&
                  s.reals[2] == linear_density() && s.reals[3] >= 0,
              "Fibre geometry or time mismatch.");
  set_initial_state({s.reals[4], s.reals[5], s.reals[6]}, {s.reals[7], s.reals[8], s.reals[9]},
                    {s.reals[10], s.reals[11], s.reals[12], s.reals[13]},
                    {s.reals[14], s.reals[15], s.reals[16]});
  set_time_integration_state(s.reals[3], static_cast<size_t>(s.ints[2]));
}

void DryRigidBodyFibre::append_vtk_polydata(std::vector<IO::CXX::vtkPolyData>& datasets) const {
  IO::CXX::vtkPolyData data;
  for (const auto& p : mesh().world_points())
    data.add_point(static_cast<float>(p[0]), static_cast<float>(p[1]), static_cast<float>(p[2]));
  for (size_t i = 1; i < nodes(); ++i)
    data.add_line(static_cast<int64_t>(i - 1), static_cast<int64_t>(i));
  datasets.push_back(std::move(data));
}
} // namespace ELFF::Models
