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

trace void output_hdf_htg(
  scalar *scalar_list, 
  vector *vector_list,
  const char *basename, 
  int iter = i, 
  double time = t
) {
#if _MPI
#if MPI_SINGLE_FILE
    char fname[64];
    sprintf(fname, "%s_%d.hdf", basename, iter);
    vtkHDFHyperTreeGrid vtk_hdf = vtk_HDF_hypertreegrid_init(scalar_list, vector_list, fname);
    vtk_HDF_hypertreegrid_close(&vtk_hdf);
#else
    char fname[64];
    sprintf(fname, "%s_pid_%d_%d.hdf", basename, pid(), iter);
    vtkHDFHyperTreeGrid vtk_hdf = vtk_HDF_hypertreegrid_init(scalar_list, vector_list, fname);
    vtk_HDF_hypertreegrid_close(&vtk_hdf);
#endif 
#else
    char fname[64];
    sprintf(fname, "%s_%d.hdf", basename, iter);
    vtkHDFHyperTreeGrid vtk_hdf = vtk_HDF_hypertreegrid_init(scalar_list, vector_list, fname);
    vtk_HDF_hypertreegrid_close(&vtk_hdf);
#endif
}