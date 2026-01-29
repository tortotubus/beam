#include <stdio.h>

int copy_file (const char* src, const char* dst) {
  FILE* in = fopen (src, "rb");
  if (!in) {
    perror ("fopen src");
    return -1;
  }

  FILE* out = fopen (dst, "wb");
  if (!out) {
    perror ("fopen dst");
    fclose (in);
    return -1;
  }

  char buf[1 << 16]; // 64 KB
  size_t n;

  while ((n = fread (buf, 1, sizeof (buf), in)) > 0) {
    if (fwrite (buf, 1, n, out) != n) {
      perror ("fwrite");
      fclose (in);
      fclose (out);
      return -1;
    }
  }

  if (ferror (in)) {
    perror ("fread");
    fclose (in);
    fclose (out);
    return -1;
  }

  fclose (in);
  fclose (out);
  return 0;
}

#include "library/vtk/vtkHDFHyperTreeGrid.h"

trace void output_hdf_htg (const char* basename,
                           scalar* scalar_list,
                           vector* vector_list,
                           int iter = i,
                           double time = t,
                           bool use_transient = true,
                           bool overwrite = true,
                           bool use_ab = false) {
#if _MPI
  for (scalar s in scalar_list)
    s.dirty = true;
  for (scalar s in vector_list)
    s.dirty = true;

  boundary (scalar_list);
  boundary (vector_list);
#endif

  if (use_ab && use_transient) {
    char fname_a[64];
    sprintf (fname_a, "%s_a.vtkhdf", basename);
    char fname_b[64];
    sprintf (fname_b, "%s_b.vtkhdf", basename);

    if (iter == 0) {
      vtkHDFHyperTreeGrid vtk_hdf = vtk_HDF_hypertreegrid_init_transient (
        scalar_list, vector_list, time, fname_a, overwrite);
      vtk_HDF_hypertreegrid_close (&vtk_hdf);
      if (pid () == 0)
        copy_file (fname_a, fname_b);
    } else {
      vtkHDFHyperTreeGrid vtk_hdf = vtk_HDF_hypertreegrid_append_transient (
        scalar_list, vector_list, time, fname_a);
      vtk_HDF_hypertreegrid_close (&vtk_hdf);
      if (pid () == 0)
        copy_file (fname_a, fname_b);
    }
  } else if (use_transient) {
    char fname[64];
    sprintf (fname, "%s.vtkhdf", basename);
    if (iter == 0) {
      vtkHDFHyperTreeGrid vtk_hdf = vtk_HDF_hypertreegrid_init_transient (
        scalar_list, vector_list, time, fname, overwrite);
      vtk_HDF_hypertreegrid_close (&vtk_hdf);
    } else {
      vtkHDFHyperTreeGrid vtk_hdf = vtk_HDF_hypertreegrid_append_transient (
        scalar_list, vector_list, time, fname);
      vtk_HDF_hypertreegrid_close (&vtk_hdf);
    }
  } else {

    char fname[128];
    sprintf (fname, "%s_%d.vtkhdf", basename, iter);
    vtkHDFHyperTreeGrid vtk_hdf = vtk_HDF_hypertreegrid_init_static (
      scalar_list, vector_list, fname, overwrite);
    vtk_HDF_hypertreegrid_close (&vtk_hdf);
 
    static int timeseries_len = 0;
    static char timeseries[1 << 26]; // 64 MB buffer

    // Tail that always lives *after* the last entry
    static const char timeseries_tail[] = ",\n\t]}";
    const int timeseries_delta = (int) strlen (timeseries_tail);

    char fname_timeseries[128];
    sprintf (
      fname_timeseries, "%s.vtkhdf.series", basename); // assuming fname exists

    if (!timeseries_len) {
      // First time: write full JSON with one entry and the tail
      timeseries_len += sprintf (timeseries,
                                 "{\n"
                                 "\t\"file-series-version\" : \"1.0\",\n"
                                 "\t\"files\" : [\n"
                                 "\t\t{ \"name\" : \"%s\", \"time\" : %f }%s",
                                 fname,
                                 time,
                                 timeseries_tail);
    } else {
      // We currently have ... { last entry }",\n\t]}"
      // The tail begins at timeseries_len - timeseries_delta
      int start = timeseries_len - timeseries_delta;

      // Overwrite the tail with: ",\n\t\t{ new entry }<tail again>"
      // This keeps all old content up to 'start', and appends a new entry
      timeseries_len =
        start + sprintf (timeseries + start,
                         ",\n\t\t{ \"name\" : \"%s\", \"time\" : %f }%s",
                         fname,
                         time,
                         timeseries_tail);
    }
 
    FILE* f = fopen (fname_timeseries, "w");
    if (f) {
      fwrite (timeseries, 1, (size_t) timeseries_len, f);
      fclose (f);
    }
  }
}
