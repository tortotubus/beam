#include "beam.h"

#include "models/beam/IBEulerBeam.hpp"
#include "config/config.hpp"

using namespace beam;

extern "C"
{

  ib_beam_t ib_beam_new(double length,
                        double EI,
                        double mu,
                        int nodes,
                        double r_penalty)
  {

    EulerBeam::EulerBeamBCs bcs = { .end = { EulerBeam::left, EulerBeam::right },
                         .type = { EulerBeam::simple_bc, EulerBeam::free_bc } };

    return new IBEulerBeam(static_cast<real_t>(length),
                           static_cast<real_t>(EI),
                           static_cast<real_t>(mu),
                           static_cast<size_t>(nodes),
                           bcs,
                           static_cast<real_t>(r_penalty));
  }

  ib_beam_t ib_beam_new_theta(double length,
                              double EI,
                              double mu,
                              int nodes,
                              double r_penalty,
                              double theta)
  {

    real_t x0 = 0, y0 = 0;

    EulerBeamMesh ic_mesh(nodes, length);

    std::vector<std::array<real_t, 3>>& ic_centerline =
      ic_mesh.get_centerline();
    std::vector<std::array<real_t, 3>>& ic_slope = ic_mesh.get_slope();
    std::vector<std::array<real_t, 3>>& ic_velocity =
      ic_mesh.get_centerline_velocity();
    std::vector<real_t>& ic_s = ic_mesh.get_curvilinear_axis();

    for (size_t ni = 0; ni < nodes; ++ni) {
      real_t s = ic_s[ni];

      ic_centerline[ni][0] = x0 + (length - s) * std::cos(theta);
      ic_centerline[ni][1] = y0 + (length - s) * std::sin(theta);
      ic_centerline[ni][2] = 0.;

      // Tangent dX/ds:
      ic_slope[ni][0] = -std::cos(theta);
      ic_slope[ni][1] = -std::sin(theta);
      ic_slope[ni][2] = 0.;

      ic_velocity[ni][0] = 0.;
      ic_velocity[ni][1] = 0.;
      ic_velocity[ni][2] = 0.;
    }

    EulerBeam::EulerBeamBCs boundary_conditions = { .end = { EulerBeam::left, EulerBeam::right },
                                         .type = { EulerBeam::free_bc, EulerBeam::simple_bc },
                                         .vals = { {
                                                     .position = { 0, 0, 0 },
                                                   },
                                                   {} } };

    IBEulerBeam* beam = new IBEulerBeam(static_cast<real_t>(length),
                                        static_cast<real_t>(EI),
                                        static_cast<real_t>(mu),
                                        static_cast<size_t>(nodes),
                                        boundary_conditions,
                                        static_cast<real_t>(r_penalty));

    beam->apply_initial_condition(ic_mesh);
    return beam;
  }

  void ib_beam_destroy(ib_beam_t handle)
  {
    delete reinterpret_cast<IBEulerBeam*>(handle);
  }
}