
#include "library/io/output-common.h"

#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#include "library/io/vtk/vtkHDFHyperTreeGrid.h"

trace void output_hdf_htg_series (scalar* scalar_list,
                                  vector* vector_list,
                                  const char* basename,
                                  int iter,
                                  double time,
                                  bool overwrite) {
  char pname[128];
  sprintf (pname, "%s-vtkhdf-series", basename);
  create_path (pname);

  char fname[128];
  sprintf (fname, "%s-vtkhdf-series/%s_%d.vtkhdf", basename, basename, iter);

  static int timeseries_len = 0;
  static char timeseries[1 << 24]; // 64 MB buffer

  // Tail that always lives *after* the last entry
  static const char timeseries_tail[] = ",\n\t]}";
  const int timeseries_delta = (int) strlen (timeseries_tail);

  char fname_timeseries[128];
  sprintf (
    fname_timeseries, "%s.vtkhdf.series", basename); // assuming fname exists

  if (!timeseries_len) {
    timeseries_len +=
      sprintf (timeseries,
               "{\n"
               "\t\"file-series-version\" : \"1.0\",\n"
               "\t\"files\" : [\n"
               "\t\t{ \"iter\" : %d, \"time\" : %f, \"name\" : \"%s\" }%s",
               iter,
               time,
               fname,
               timeseries_tail);
  } else {
    // We currently have ... { last entry }",\n\t]}"
    // The tail begins at timeseries_len - timeseries_delta
    int start = timeseries_len - timeseries_delta;

    // Overwrite the tail with: ",\n\t\t{ new entry }<tail again>"
    // This keeps all old content up to 'start', and appends a new entry
    timeseries_len =
      start +
      sprintf (timeseries + start,
               ",\n\t\t{ \"iter\" : %d, \"time\" : %f, \"name\" : \"%s\" }%s",
               iter,
               time,
               fname,
               timeseries_tail);
  }

  FILE* f = fopen (fname_timeseries, "w");
  if (f) {
    fwrite (timeseries, 1, (size_t) timeseries_len, f);
    fclose (f);
  }

  vtkHDFHyperTreeGrid vtk_hdf = vtk_HDF_hypertreegrid_init_static (
    scalar_list, vector_list, fname, overwrite);
  vtk_HDF_hypertreegrid_close (&vtk_hdf);
}

trace void output_hdf_htg (scalar* scalar_list = NULL,
                           vector* vector_list = NULL,
                           const char* basename = NULL,
                           int iter = i,
                           double time = t,
                           bool use_transient = false,
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

  scalar* slist = scalar_list;
  vector* vlist = vector_list;

  if (!slist) {
    for (scalar s in all) {
      if (s.v.x.i == s.i) {
        vector v = s.v;
        vlist = vectors_add (vlist, v);
      } else if (s.v.x.i < 0) {
        slist = list_add (slist, s);
      }
    }
  }

  char basename_buff[64];

  if (!basename) {
    snprintf (
      basename_buff, sizeof (basename_buff), "%s-htg", get_unique_basename ());
  } else {
    snprintf (
      basename_buff, sizeof (basename_buff), "%s-htg", basename);
  }

  output_hdf_htg_series (
    slist, vlist, basename_buff, iter, time, overwrite);
}
