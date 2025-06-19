#define BGHOSTS 2

#include "navier-stokes/centered.h"

#include "ib_mesh.h"

#define POS_PBC_X(X)                                                           \
  ((u.x.boundary[left] != periodic_bc)                                         \
     ? (X)                                                                     \
     : (((X - (X0 + L0 / 2)) > L0 / 2.) ? (X) - L0 : (X)))
#define POS_PBC_Y(Y)                                                           \
  ((u.x.boundary[top] != periodic_bc)                                          \
     ? (Y)                                                                     \
     : (((Y - (Y0 + L0 / 2)) > L0 / 2.) ? (Y) - L0 : (Y)))
#define POS_PBC_Z(Z)                                                           \
  ((u.x.boundary[front] != periodic_bc)                                        \
     ? (Z)                                                                     \
     : (((Z - (Z0 + L0 / 2)) > L0 / 2.) ? (Z) - L0 : (Z)))

void generate_stencils (IBMesh *mesh) {
  for (size_t ni = 0; ni < mesh->n_nodes; ni++) {
    mesh->nodes[ni].stencil.n = 0;
    double delta = (L0/(1 << grid->maxdepth));
    for (int mi = -2; mi <= 2; mi++) {
      for (int mj = -2; mj <= 2; mj++) {
        #if dimension == 2
          Point point = locate(
            POS_PBC_X(mesh->nodes[ni].pos.x + mi*delta),
            POS_PBC_Y(mesh->nodes[ni].pos.y + mj*delta)
          ); 
        #else // dimension == 3
          for (int mk = -2; mk <= 2; mk++) {
            Point point = locate(
              POS_PBC_X(mesh->nodes[ni].pos.x + mi*delta),
              POS_PBC_Y(mesh->nodes[ni].pos.y + mj*delta),
              POS_PBC_Z(mesh->nodes[ni].pos.z + mk*delta)
            );
          }
        #endif
      }
    }
  }
}

trace void
lag2eul(vector forcing, IBMesh* mesh)
{
  for (int ni = 0; ni < mesh->n_nodes; ni++) {
    foreach_cache(mesh->nodes[ni].stencil)
    {
      if (point.level >= 0) {
        coord dist;
#if dimension == 1
        //
#elif dimension == 2
        dist.x = GENERAL_1DIST(x, mesh->nodes[ni].pos.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[ni].pos.y);
        if (fabs(dist.x) <= 2 * Delta && fabs(dist.y) <= 2 * Delta) {
          double weight = (1 + cos(.5 * pi * dist.x / Delta)) *
                          (1 + cos(.5 * pi * dist.y / Delta)) / (sq(4 * Delta));

          foreach_dimension() forcing.x[] += weight * mesh->nodes[ni].force.x;
        }
#else // dimension == 3
        dist.x = GENERAL_1DIST(x, mesh->nodes[ni].pos.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[ni].pos.y);
        dist.z = GENERAL_1DIST(z, mesh->nodes[ni].pos.z);
        if (fabs(dist.x) <= 2 * Delta && fabs(dist.y) <= 2 * Delta &&
            fabs(dist.z) <= 2 * Delta) {
          double weight = (1 + cos(.5 * pi * dist.x / Delta)) *
                          (1 + cos(.5 * pi * dist.y / Delta)) *
                          (1 + cos(.5 * pi * dist.z / Delta)) /
                          (cube(4 * Delta));

          foreach_dimension() forcing.x[] += weight * mesh->nodes[ni].force.x;
        }
#endif
      }
    }
  }
}

trace void
eul2lag(IBMesh* mesh)
{
  for (int ni = 0; ni < mesh->n_nodes; ni++) {
    foreach_dimension() mesh->nodes[ni].velocity.x = 0.;
    foreach_cache(mesh->nodes[ni].stencil)
    {
      if (point.level >= 0) {
        coord dist;
#if dimension == 1
        // Ignore
#elif dimension == 2
        dist.x = GENERAL_1DIST(x, mesh->nodes[ni].pos.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[ni].pos.y);
        if (fabs(dist.x) <= 2 * Delta && fabs(dist.y) <= 2 * Delta) {
          double weight = (1 + cos(.5 * pi * dist.x / Delta)) *
                          (1 + cos(.5 * pi * dist.y / Delta)) / 16.;
          foreach_dimension() mesh->nodes[ni].velocity.x += weight * u.x[];
        }
#else // dimension == 3
        dist.x = GENERAL_1DIST(x, mesh->nodes[ni].pos.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[ni].pos.y);
        dist.z = GENERAL_1DIST(z, mesh->nodes[ni].pos.z);
        if (fabs(dist.x) <= 2 * Delta && fabs(dist.y) <= 2 * Delta &&
            fabs(dist.z) <= 2 * Delta) {
          double weight = (1 + cos(.5 * pi * dist.x / Delta)) *
                          (1 + cos(.5 * pi * dist.y / Delta)) *
                          (1 + cos(.5 * pi * dist.z / Delta)) / 64.;
          foreach_dimension() mesh->nodes[ni].velocity.x += weight * u.x[];
        }
#endif
      }
    }
  }
}