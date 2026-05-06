/* Mesh */
#include "grid/quadtree.h"

/* Immersed-boundary framework */
#include "library/ibm/IBMeshManager.h"
#include "library/ibm/IBAdapt.h"

/* Input/output headers */
#include "library/ibm/IBOutput.h"
#include "library/io/output-vtk.h"

coord circle(int n, int N, coord centre, double radius) {
  double rad = 2. * pi * ((double)n / (double)N);
  coord c = centre;
  c.x += radius * cos(rad);
  c.y += radius * sin(rad);
  return c;
}

int main() {
  init_grid(1 << 5);
  L0 = 4;
  origin(-L0 / 2., -L0 / 2., -L0 / 2.);

  int N_ib = 10;
  coord centre = {0};
  ibmeshmanager_init(1);
  ibmeshmanager_add_nodes(0, N_ib);
  foreach_ibnode() {
    coord c_pos = circle(node_id, N_ib, centre, 0.5);
    foreach_dimension() { ibval(npos.x) = c_pos.x; }
  }
  foreach_ibnode_per_ibmesh() {
    mesh->depth = 10;
    node->depth = 10;
  }

  adapt_wavelet_ibm(NULL, NULL, 5);

  {
    int i = 0; double t = 0.;
    output_hdf_htg(NULL,NULL,"adapt_wavelet_test");
    output_hdf_pd(NULL,NULL,"adapt_wavelet_test");
  }

  return 0;
}
