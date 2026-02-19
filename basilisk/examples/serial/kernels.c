#include "grid/quadtree.h"

#include "library/ibm/IBAdapt.h"
#include "library/ibm/IBMeshManager.h"

#include "library/ibm/IBOutput.h"
#include "library/io/vtk/vtkHDFHyperTreeGrid.h"

const int maxlevel = 5;
const int minlevel = 2;

scalar noisef[];
scalar levelf[];

coord
circle(int n, int N, coord centre, double radius)
{
  double rad = 2. * pi * ((double)n / (double)N);
  coord c = centre;
  c.x += radius * cos(rad);
  c.y += radius * sin(rad);
  return c;
}

int
main()
{
  init_grid(1 << minlevel);
  X0 = -1.85, Y0 = -4, L0 = 8;

  periodic(right);
  periodic(top);

  // Create mesh manager object with 0 meshes
  ibmeshmanager_init(0);

  // Add meshes
  int new_id = ibmeshmanager_add_mesh();

  // Add a single node to the first mesh
  const int N_circ = 100;
  const double r_circ = 0.15;
  const coord cen_circ = { 0 };

  ibmeshmanager_add_nodes(new_id, N_circ);

  foreach_ibmesh()
  {
    mesh->refinement_level = 10;
  }

  foreach_ibnode_per_ibmesh()
  {
    node->lagpos = circle(node_id, N_circ, cen_circ, r_circ);
  }

  for (int i = 0; i < maxlevel - minlevel; i++) {
    // Inject noise with noise
    foreach_cell()
    {
      noisef[] = noise();
    }

    // Refine and coarsen around nodes
    // adapt_wavelet_ibm({noisef}, (double[]){1e-1}, maxlevel, minlevel);
    adapt_wavelet_ibm(NULL, NULL, maxlevel, minlevel);
    printf("\n\n");
  }

  foreach () {
    levelf[] = point.level;
  }

  // Output
  vtkHDFHyperTreeGrid vtkhdf = vtk_HDF_hypertreegrid_init_static(
    { levelf, noisef }, NULL, "kernels_test.vtkhdf", true);
  vtk_HDF_hypertreegrid_close(&vtkhdf);


  output_ibnodes("kernels_test.vtk", 0, 0);

  ibmeshmanager_free();

  return 0;
}
