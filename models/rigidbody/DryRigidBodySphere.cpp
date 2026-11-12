#include "elff/models/rigidbody/DryRigidBodySphere.hpp"
#include <cmath>
#include <limits>

namespace ELFF::Models {
namespace {
void finite_vector(const DryRigidBodySphere::Vec3& vector)
{
  for (auto x : vector)
    ELFF_VERIFY(std::isfinite(x), "Dry sphere vector must be finite.");
}
}

DryRigidBodySphere::DryRigidBodySphere(real_t radius, real_t density,
    size_t point_count, const Vec3& center, const Vec3& velocity,
    const Quat& orientation, const Vec3& angular_velocity)
  : RigidBodySphere(radius, point_count, density)
{
  ELFF_VERIFY(std::isfinite(radius) && std::isfinite(density)
              && std::isfinite(mass()), "Dry sphere parameters must be finite.");
  finite_vector(center);
  finite_vector(velocity);
  finite_vector(angular_velocity);
  real_t norm = 0;
  for (auto q : orientation) {
    ELFF_VERIFY(std::isfinite(q), "Invalid quaternion.");
    norm += q * q;
  }
  ELFF_VERIFY(std::isfinite(norm) && norm > 0, "Invalid quaternion norm.");
  const real_t inertia = (real_t{2} / 5) * mass() * radius * radius;
  for (auto omega : angular_velocity)
    ELFF_VERIFY(std::isfinite(inertia * omega), "Angular momentum overflow.");
  set_initial_state(center, velocity, orientation, angular_velocity);
}

void DryRigidBodySphere::pack_state(ModelState& s) const
{
  ELFF_VERIFY(point_count() <= static_cast<size_t>(std::numeric_limits<int64_t>::max())
              && time_iteration() <= static_cast<size_t>(std::numeric_limits<int64_t>::max()),
              "Sphere checkpoint integer overflow.");
  s = {};
  s.ints = {1, static_cast<int64_t>(point_count()), static_cast<int64_t>(time_iteration())};
  s.reals = {radius(), density(), time()};
  for (const auto& v : {center_of_mass(), com_velocity()})
    s.reals.insert(s.reals.end(), v.begin(), v.end());
  s.reals.insert(s.reals.end(), pose().begin(), pose().end());
  s.reals.insert(s.reals.end(), angular_velocity_body().begin(), angular_velocity_body().end());
}

void DryRigidBodySphere::unpack_state(const ModelState& s)
{
  ELFF_VERIFY(s.ints.size() == 3 && s.ints[0] == 1 && s.ints[1] > 0
              && static_cast<uint64_t>(s.ints[1]) == point_count() && s.ints[2] >= 0
              && static_cast<uint64_t>(s.ints[2]) <= std::numeric_limits<size_t>::max()
              && s.reals.size() == 16 && s.bytes.empty(), "Incompatible dry sphere state.");
  for (auto x : s.reals) ELFF_VERIFY(std::isfinite(x), "Non-finite dry sphere state.");
  ELFF_VERIFY(s.reals[0] == radius() && s.reals[1] == density() && s.reals[2] >= 0,
              "Dry sphere geometry or time mismatch.");
  real_t norm = 0;
  for (size_t i = 9; i < 13; ++i) norm += s.reals[i] * s.reals[i];
  ELFF_VERIFY(std::isfinite(norm) && norm > 0, "Invalid saved orientation.");
  const real_t inertia = (real_t{2} / 5) * mass() * radius() * radius();
  for (size_t i = 13; i < 16; ++i)
    ELFF_VERIFY(std::isfinite(inertia * s.reals[i]), "Saved angular momentum overflow.");
  set_initial_state({s.reals[3], s.reals[4], s.reals[5]},
    {s.reals[6], s.reals[7], s.reals[8]},
    {s.reals[9], s.reals[10], s.reals[11], s.reals[12]},
    {s.reals[13], s.reals[14], s.reals[15]});
  set_time_integration_state(s.reals[2], static_cast<size_t>(s.ints[2]));
}

void DryRigidBodySphere::append_vtk_polydata(std::vector<IO::CXX::vtkPolyData>& datasets) const
{
  IO::CXX::vtkPolyData data;
  for (const auto& p : mesh().world_points()) {
    data.add_point(static_cast<float>(p[0]), static_cast<float>(p[1]), static_cast<float>(p[2]));
  }
  for (size_t i = 0; i < point_count(); ++i) data.add_vertex(static_cast<int64_t>(i));
  datasets.push_back(std::move(data));
}

} // namespace ELFF::Models
