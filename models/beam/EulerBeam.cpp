#include "elff/models/beam/EulerBeam.hpp"

namespace ELFF {
namespace Models {

void
EulerBeam::solve()
{
  ELFF_ABORT("EulerBeam base class does not implement a solution method.\n");
}

void
EulerBeam::solve(std::array<real_t, 3> uniform_load)
{
  ELFF_ABORT("EulerBeam base class does not implement a solution method.\n");
}

void
EulerBeam::solve(std::vector<std::array<real_t, 3>> nonuniform_load)
{
  ELFF_ABORT("EulerBeam base class does not implement a solution method.\n");
}

void
EulerBeam::solve(real_t dt)
{
  ELFF_ABORT("EulerBeam base class does not implement a solution method.\n");
}

void
EulerBeam::solve(real_t dt, std::array<real_t, 3> uniform_load)
{
  ELFF_ABORT("EulerBeam base class does not implement a solution method.\n");
}

void
EulerBeam::solve(real_t dt, std::vector<std::array<real_t, 3>> nonuniform_load)
{
  ELFF_ABORT("EulerBeam base class does not implement a solution method.\n");
}

void
EulerBeam::set_time_dependent_boundary_conditions(EulerBeamTimeDependentBCs bcs)
{
  time_dependent_boundary_conditions = bcs;
  have_time_dependent_boundary_conditions = true;
}


EulerBeamMesh&
EulerBeam::get_mesh()
{
  return mesh;
}

const EulerBeamMesh&
EulerBeam::get_mesh() const
{
  return mesh;
}

EulerBeam::EulerBeam()
  : EI(1)
  , mu(0)
  , is_time_dependent(false)
  , time_iter(0)
  , t(0)
  , dimension(3)
  , mesh(20, 1.)
  , boundary_conditions(
      { .end = { left, right },
        .type = { clamped_bc, free_bc },
        .vals = { { .position = { 0., 0., 0. }, .slope = { 1., 0., 0. } },
                  { .position = { 0., 0., 0. }, .slope = { 0., 0., 0. } } } })
  , time_dependent_boundary_conditions()
  , have_time_dependent_boundary_conditions(false)
{
}

EulerBeam::EulerBeam(real_t length, real_t EI, size_t nodes, EulerBeamBCs bcs)
  : EI(EI)
  , mu(0)
  , is_time_dependent(false)
  , time_iter(0)
  , t(0)
  , dimension(3)
  , mesh(nodes, length)
  , boundary_conditions(bcs)
  , time_dependent_boundary_conditions()
  , have_time_dependent_boundary_conditions(false)
{
}

EulerBeam::EulerBeam(real_t length,
                     real_t EI,
                     real_t mu,
                     size_t nodes,
                     EulerBeamBCs bcs)
  : EI(EI)
  , mu(mu)
  , is_time_dependent(true)
  , time_iter(0)
  , t(0)
  , dimension(3)
  , mesh(nodes, length)
  , boundary_conditions(bcs)
  , time_dependent_boundary_conditions()
  , have_time_dependent_boundary_conditions(false)
{
}

EulerBeam::EulerBeam(real_t length,
                     real_t EI,
                     real_t mu,
                     size_t nodes,
                     EulerBeamTimeDependentBCs bcs)
  : EI(EI)
  , mu(mu)
  , is_time_dependent(true)
  , time_iter(0)
  , t(0)
  , dimension(3)
  , mesh(nodes, length)
  , boundary_conditions()
  , time_dependent_boundary_conditions(bcs)
  , have_time_dependent_boundary_conditions(true)
{
  for (size_t bi = 0; bi < 2; ++bi) {
    boundary_conditions.end[bi] = bcs.end[bi];
    boundary_conditions.type[bi] = bcs.type[bi];
    boundary_conditions.vals[bi] = {};
  }
}

} // namespace Models
} // namespace ELFF
