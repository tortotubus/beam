vector u[];
scalar s[];
scalar p[];

#include "grid/quadtree.h"
#include "library/ibm/kernels.h"
#include "library/output_htg.h"

int
main()
{
  init_grid(1 << 7);
  L0 = 8;
  origin(-L0 / 2, -L0 / 2);

  periodic(right);
  periodic(top);

  if (pid() == 0) {
    printf("Period x: %d\n", Period.x);
    printf("Period y: %d\n", Period.y);
    printf("Period z: %d\n", Period.z);
  }

  foreach () {
    s[] = 0.;
    p[] = cell.pid;
    foreach_dimension()
    {
      u.x[] = 0.;
    }
  }

  coord c;

  c.x = L0;
  c.y = 0;
  c.z = 0;

#define D 2

  int isum = 0; 
  // foreach (reduction(+:isum))
  //   isum++;

  // foreach_cell() {
  //   if (fabs(c.x-x) <= Delta * D && fabs(c.y-y) <= Delta * D) {
  //     s[] = (1. + cos(pi * fabs(c.x-x) / (Delta * D))) *
  //           (1. + cos(pi * fabs(c.y-y) / (Delta * D))) / sq(2. * Delta * D);
  //   } else {
  //     isum++;
  //   }
  // } 

  // foreach_neighborhood(c,2) {
  //   s[] = 2.;
  //   isum++;
  // }

  if (isum != 0)
    printf("[pid %d] i: %d\n", pid(), isum);

  output_hdf_htg({ s, p }, { u }, "pbs_mpi", 0, 0);

  return 0;
}