#include <stdlib.h>

#include "elff/c/models/beam/IBEulerBeam.h"
#include "elff/c/models/ibm/IBForceCoupled.h"
#include "elff/c/models/ModelRuntime.h"

static ib_force_coupled_t beam;

int
main()
{
  model_runtime_t runtime = model_runtime_new();
  {
    const double L = 1.;
    const double EI = 1.;
    const double mu = 1.;
    const int Nb = 20;
    const double r_penalty = 1e4;
    const double theta = 0.;
    ib_euler_beam_bcs_t bcs = {0};
    bcs.end[0] = IB_EULER_BEAM_BC_LEFT;
    bcs.end[1] = IB_EULER_BEAM_BC_RIGHT;
    bcs.type[0] = IB_EULER_BEAM_BC_FREE;
    bcs.type[1] = IB_EULER_BEAM_BC_SIMPLE;
    const vertex_t origin = {0};
    beam = (ib_force_coupled_t)ib_euler_beam_new_theta(
      runtime, origin, bcs, L, EI, mu, Nb, r_penalty, theta);
  }

  double dt = 0.01;
  double t = 0.;

  for (int i = 0; i < 100; i++) {
    ib_mesh_t cur = ib_force_coupled_get_current(beam);
    vertex_t* force = calloc(cur.n, sizeof(vertex_t));
    for (int ni = 0; ni < cur.n; ++ni) {
      force[ni].x = 1.;
      force[ni].y = 0.;
      force[ni].z = 0.;
    }
    ib_mesh_t next = ib_force_coupled_get_next(beam, force, cur.n, dt);
    ib_mesh_free(&cur);
    ib_mesh_free(&next);
    free(force);
    t += dt;
  }

  model_runtime_delete(runtime);

  return 0;
}
