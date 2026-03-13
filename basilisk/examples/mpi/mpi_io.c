#define COMPRESSION 1

#include "grid/multigrid3D.h"
#include "library/io/output-vtk.h"

scalar pidf[];

int
main()
{
  dimensions(nx = 3, ny = 1, nz = 1);
  L0 = 3;
  init_grid(1 << 4);

  foreach () {
    pidf[] = pid();
  }

  foreach () {
    foreach_neighbor()
    {
      if (point.i == 0) {
        printf("mpiboundary\n");
      }
    }
  }

  {
    int i = 0;
    double t = 0.;
    output_hdf_imagedata();
  }

  return 0;
}