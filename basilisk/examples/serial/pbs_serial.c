vector u[];
scalar s[];
scalar p[];

#include "grid/quadtree.h"
#include "library/ibm/kernels.h"
#include "library/output_htg.h"

int
main()
{
  init_grid(1 << 5);
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
    p[] = pid();
    foreach_dimension()
    {
      u.x[] = 0.;
    }
  }


  coord c;

  c.x = 2;
  c.y = 2;
  c.z = 0;

 
  foreach_neighbor_of_coord(c, 2)
  {
    s[] = -4; 
  }

  foreach_point(-2,-2) {
    s[] = 4;
    foreach_neighbor(1) {
      printf("%d\n", level);
    }
  }

  output_hdf_htg({ s, p }, { u }, "pbs_serial", 0, 0);

  return 0;
}