#include "grid/quadtree.h"

#include "library/ibm/IBAdapt.h"
#include "library/ibm/IBKernels.h"
#include "library/ibm/IBMeshManager.h"

#include "library/io/vtk/vtkHDFHyperTreeGrid.h"
#include "library/io/vtk/vtkHDFPolyData.h"
#include "library/io/vtk/vtkPolyData.h"

const int maxlevel = 6;
const int minlevel = 5;
const int ibmlevel = 10;

scalar suppf[];
scalar weightf[];
scalar levelf[];
scalar pidf[];
scalar gradientf[];

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
  const int N_circ = 100;
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

  //
  // Begin filling vtkPolyData
  //

  vtkPolyData pd = vtk_polydata_init(N_circ, N_circ, 0, 0, 0, 1);

  foreach_ibnode()
  {
#if dimension == 1
    vtk_polydata_add_point(&pd, node->lagpos.x, 0., 0.);
#elif dimension == 2
    vtk_polydata_add_point(&pd, node->lagpos.x, node->lagpos.y, 0.);
#else
    vtk_polydata_add_point(&pd, node->lagpos.x, node->lagpos.y, node->lagpos.z);
#endif
  }

  // foreach_ibnode()
  // {
  //   vtk_polydata_add_vertex(&pd, node_id);
  // }

  int64_t pntdata_id = vtk_polydata_add_pointdata_scalar(&pd, "pid");
  printf("pntdata_id: %d\n", pntdata_id);
  double *pntdata_arr = vtk_polydata_get_pointdata_data(&pd, pntdata_id);
  printf("pntdata_arr: %d\n", pntdata_arr);
  int pntdata_arr_idx = 0;

  foreach_ibnode() {
    pntdata_arr[pntdata_arr_idx] = pntdata_arr_idx; //double) node->pid;
    printf("pntdata_arr[%d] = %f\n", pntdata_arr_idx, pntdata_arr[pntdata_arr_idx]);
    pntdata_arr_idx++;
  }

  //
  // Write polydata
  //

  vtkHDFPolyData hdf_pd = vtk_HDF_polydata_init_static("vtkhdfpolydata.vtkhdf", true, &pd);
  vtk_HDF_polydata_close(&hdf_pd);

  vtk_polydata_free(&pd);

  ibmeshmanager_free();

  return 0;
}
