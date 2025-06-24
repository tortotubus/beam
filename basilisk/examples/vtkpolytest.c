#include "library/vtk/IO/HDF/vtkHDFPolyData.h"
#define dimension 3

int
main(void)
{
  vtkPolyData_t pd;
  vtkPolyData_init(&pd);

  // 1) Initialize 3 points in 3D
  if (vtkPoints_init(&pd.points, 3, 3) != 0) {
    fprintf(stderr, "Failed to alloc points\n");
    return 1;
  }

  // 2) Fill in the point coordinates (triangle in XY plane)
  float coords[9] = {
    0.0f, 0.0f, 0.0f, // point 0
    1.0f, 0.0f, 0.0f, // point 1
    0.0f, 1.0f, 0.0f  // point 2
  };
  memcpy(pd.points.data, coords, sizeof(coords));

  // 3) Create 3 VTK_VERTEX cells (each of size 1)
  int64_t vertex_sizes[3] = { 1, 1, 1 };
  if (vtkCellArray_init(&pd.vertices, 3, vertex_sizes) != 0) {
    fprintf(stderr, "Failed to alloc vertex cell array\n");
    vtkPolyData_free(&pd);
    return 1;
  }
  // connectivity just lists the point‐indices
  pd.vertices.connectivity[0] = 0;
  pd.vertices.connectivity[1] = 1;
  pd.vertices.connectivity[2] = 2;

  // --- At this point, 'pd' holds three points and three vertex‐cells. ---

  // Print the points
  for (int i = 0; i < pd.points.n_points; i++) {
    float* p = &pd.points.data[i * pd.points.n_components];
    printf("Point %d: (%f, %f, %f)\n", i, p[0], p[1], p[2]);
  }

  // Print the vertex cells using new offsets (length = n_cells+1)
  for (int64_t c = 0; c < pd.vertices.n_cells; c++) {
    int64_t start = pd.vertices.offsets[c];
    int64_t end = pd.vertices.offsets[c + 1];
    printf("Vertex cell %lld contains points:", (long long)c);
    for (int64_t idx = start; idx < end; idx++) {
      printf(" %lld", (long long)pd.vertices.connectivity[idx]);
    }
    printf("\n");
  }

  // Write out to HDF5
  vtkHDFPolyData hdf_pd = vtk_HDF_polydata_init(&pd, "polyvertstest.vtkhdf");
  vtk_HDF_polydata_close(&hdf_pd);

  // 4) Clean up
  vtkPolyData_free(&pd);
  return 0;
}
