#if _MPI
#define MPI_SINGLE_FILE 1
#else
#define MPI_SINGLE_FILE 0
#endif

#include "./vtk/vtkHDFHyperTreeGrid.h"

trace void output_hdf_htg (scalar* scalar_list,
                           vector* vector_list,
                           const char* basename,
                           int iter = i,
                           double time = t) {
#if _MPI
  for (scalar s in scalar_list)
    s.dirty = true;
  for (scalar s in vector_list)
    s.dirty = true;

  boundary (scalar_list);
  boundary (vector_list);

#if MPI_SINGLE_FILE
  char fname[64];
  sprintf (fname, "%s_%d.vtkhdf", basename, iter);
  vtkHDFHyperTreeGrid vtk_hdf =
    vtk_HDF_hypertreegrid_init (scalar_list, vector_list, fname);
  vtk_HDF_hypertreegrid_close (&vtk_hdf);
#else
  char fname[64];
  sprintf (fname, "%s_pid_%d_%d.vtkhdf", basename, pid (), iter);
  vtkHDFHyperTreeGrid vtk_hdf =
    vtk_HDF_hypertreegrid_init (scalar_list, vector_list, fname);
  vtk_HDF_hypertreegrid_close (&vtk_hdf);
#endif
#else
  char fname[64];
  sprintf (fname, "%s_%d.vtkhdf", basename, iter);
  vtkHDFHyperTreeGrid vtk_hdf =
    vtk_HDF_hypertreegrid_init (scalar_list, vector_list, fname);
  vtk_HDF_hypertreegrid_close (&vtk_hdf);

  sprintf (fname, "%s_poly_%d.vtkhdf", basename, iter);
#endif
}
