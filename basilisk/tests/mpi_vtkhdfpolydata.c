#include "grid/quadtree.h"

#include "library/ibm/IBAdapt.h"
#include "library/ibm/IBKernels.h"
#include "library/ibm/IBMeshManager.h"

// #include "library/io/vtk/vtkHDFPolyData.h"
// #include "library/io/vtk/vtkPolyData.h"

#include "library/io/output-vtk.h"

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
  const int N_circ = 20;
  const double r_circ = 0.15;
  coord cen_circ = { 0.0 };
  cen_circ.x = 3 * r_circ * pid();

  ibmeshmanager_add_nodes(new_id, N_circ);

  foreach_ibmesh()
  {
    node->depth = ibmlevel;
  }

  foreach_ibnode_per_ibmesh()
  {
    foreach_dimension()
    {
      ibval(npos.x) = circle(node_id, N_circ, cen_circ, r_circ).x;
    }
  }

  //
  // Begin filling vtkPolyData
  //

  //   vtkPolyData pd = vtk_polydata_init(N_circ, N_circ, 0, 0, 0, 1);

  //   foreach_ibnode()
  //   {
  // #if dimension == 1
  //     vtk_polydata_add_point(&pd, node->pos.x, 0., 0.);
  // #elif dimension == 2
  //     vtk_polydata_add_point(&pd, node->pos.x, node->pos.y, 0.);
  // #else
  //     vtk_polydata_add_point(&pd, node->pos.x, node->pos.y, node->pos.z);
  // #endif
  //   }

  //   foreach_ibnode()
  //   {
  //     vtk_polydata_add_vertex(&pd, node_id);
  //   }

  // int64_t pntdata_id = vtk_polydata_add_pointdata_scalar(&pd, "pid");
  // double* pntdata_arr = vtk_polydata_get_pointdata_data(&pd, pntdata_id);
  // int pntdata_arr_idx = 0;

  // foreach_ibnode()
  // {
  //   pntdata_arr[pntdata_arr_idx] = pntdata_arr_idx; // double) node->pid;
  //   pntdata_arr_idx++;
  // }

  //
  // Write polydata
  //

  // vtkHDFPolyData hdf_pd =
  // vtk_HDF_polydata_init_static("vtkhdfpolydata.vtkhdf", true, &pd);
  // vtk_HDF_polydata_close(&hdf_pd);

  // vtk_polydata_free(&pd);

  output_hdf_pd(NULL, 0, 0, false, true, false);

  ibmeshmanager_free();

  return 0;
}
