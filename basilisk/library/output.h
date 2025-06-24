#ifndef MPI_SINGLE_FILE
#if _MPI
#define MPI_SINGLE_FILE 1
#else
#define MPI_SINGLE_FILE 0
#endif
#else
#if _MPI == 0
#define MPI_SINGLE_FILE 0
#endif
#endif

#include "./vtk/IO/HDF/vtkHDFHyperTreeGrid.h"

trace void
output_hdf_htg(scalar* scalar_list,
               vector* vector_list,
               const char* basename,
               int iter = i,
               double time = t)
{
#if _MPI
#if MPI_SINGLE_FILE
  char fname[64];
  sprintf(fname, "%s_%d.hdf", basename, iter);
  vtkHDFHyperTreeGrid vtk_hdf =
    vtk_HDF_hypertreegrid_init(scalar_list, vector_list, fname);
  vtk_HDF_hypertreegrid_close(&vtk_hdf);
#else
  char fname[64];
  sprintf(fname, "%s_pid_%d_%d.hdf", basename, pid(), iter);
  vtkHDFHyperTreeGrid vtk_hdf =
    vtk_HDF_hypertreegrid_init(scalar_list, vector_list, fname);
  vtk_HDF_hypertreegrid_close(&vtk_hdf);
#endif
#else
  char fname[64];
  sprintf(fname, "%s_%d.hdf", basename, iter);
  vtkHDFHyperTreeGrid vtk_hdf =
    vtk_HDF_hypertreegrid_init(scalar_list, vector_list, fname);
  vtk_HDF_hypertreegrid_close(&vtk_hdf);

  sprintf(fname, "%s_poly_%d.hdf", basename, iter);
#endif
}

#include "library/vtk/IO/HDF/vtkHDFPolyData.h"
#include "library/ibm.h"

trace void
output_hdf_polydata(char* basename, int iter = i, double time = t)
{
  int pd_dim = 3;
  vtkPolyData_t vtk_polydata;
  vtkPolyData_init(&vtk_polydata);

  // Count total number of points
  int nn_all = 0;
  for (int mi = 0; mi < ib_mesh_manager->nm; mi++) {
    nn_all += IBMESH(mi).nn;
  }
  vtkPoints_init(&vtk_polydata.points, nn_all, pd_dim);

  // Build points dataset
  int ni_all = 0;
  for (int mi = 0; mi < ib_mesh_manager->nm; mi++) {
    for (int ni = 0; ni < IBMESH(mi).nn; ni++) {
      vtk_polydata.points.data[ni_all + 0] = IBMESH(mi).nodes[ni].pos.x;
      vtk_polydata.points.data[ni_all + 1] = IBMESH(mi).nodes[ni].pos.y;
      vtk_polydata.points.data[ni_all + 2] = IBMESH(mi).nodes[ni].pos.z;
      ni_all += pd_dim;
    }
  }
  
  // Build vertices group
  vtk_polydata.vertices.n_cells = nn_all;
  vtk_polydata.vertices.connectivity = malloc(nn_all * sizeof(int64_t));
  vtk_polydata.vertices.offsets = malloc((nn_all + 1) * sizeof(int64_t));

  for (int vi = 0; vi <= nn_all; vi++) {
    if (vi < nn_all)
      vtk_polydata.vertices.connectivity[vi] = vi;

    vtk_polydata.vertices.offsets[vi] = vi;
  }

  char fname[64];
  sprintf(fname, "%s_polydata_%d.hdf", basename, iter);

  // Write out to HDF5
  vtkHDFPolyData hdf_pd = vtk_HDF_polydata_init(&vtk_polydata, fname);
  vtk_HDF_polydata_close(&hdf_pd);

  // Free built data
  vtkPolyData_free(&vtk_polydata);
}