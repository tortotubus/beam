#include "grid/quadtree.h"

#include "library/ibm/IBAdapt.h"
#include "library/ibm/IBKernels.h"
#include "library/ibm/IBMeshManager.h"

#include "library/ibm/IBOutput.h"
#include "library/io/vtk/vtkHDFHyperTreeGrid.h"

const int maxlevel = 6;
const int minlevel = 5;
const int ibmlevel = 10;

scalar suppf[];
scalar weightf[];
scalar levelf[];
scalar pidf[];

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
  X0 = -4, Y0 = -4, L0 = 8;

  periodic(right);
  periodic(top);

  // Create mesh manager object with 0 meshes
  ibmeshmanager_init(0);

  // Add meshes
  int new_id = ibmeshmanager_add_mesh();

  // Add a single node to the first mesh
  const int N_circ = 101;
  const double r_circ = 0.15;
  const coord cen_circ = { 0.0 };

  ibmeshmanager_add_nodes(new_id, N_circ);

  foreach_ibmesh()
  {
    mesh->refinement_level = ibmlevel;
  }

  foreach_ibnode_per_ibmesh()
  {
    node->lagpos = circle(node_id, N_circ, cen_circ, r_circ);
  }

  // Refine and coarsen around nodes
  adapt_wavelet_ibm(NULL, NULL, maxlevel, minlevel);

  ibmeshmanager_init_stencil_caches();

  foreach () {
    levelf[] = point.level;
    suppf[] = 0.;
    weightf[] = 0.;
    pidf[] = pid();
  }

  foreach_ibnode()
  {
    peskin_cosine_kernel_dimensionless(node)
    {
      suppf[] = pid();
      weightf[] += weight;
    }
  }

  // Output
  vtkHDFHyperTreeGrid vtkhdf = vtk_HDF_hypertreegrid_init_static(
    { levelf, weightf, suppf, pidf }, NULL, "mpi_kernels_test.vtkhdf", true);
  vtk_HDF_hypertreegrid_close(&vtkhdf);

  output_ibnodes("mpi_kernels_test.vtk", 0, 0);

  ibmeshmanager_free();

  return 0;
}
