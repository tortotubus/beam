#include "library/io/vtk/vtkHDFHyperTreeGrid.h"

#include "library/io/output-common.h"

#include <ctype.h>
#include <string.h>
#include <stdlib.h>

trace void output_hdf_htg_series (scalar* scalar_list,
                                  vector* vector_list,
                                  const char* basename,
                                  int iter,
                                  double time,
                                  bool overwrite) {
// #if _MPI
//   for (scalar s in scalar_list)
//     s.dirty = true;
//   for (scalar s in vector_list)
//     s.dirty = true;
//   boundary (scalar_list);
//   boundary (vector_list);
// #endif

  char pname[128];
  sprintf (pname, "%s-vtkhdf-series", basename);
  create_path (pname);

  char fname[128];
  sprintf (fname, "%s-vtkhdf-series/%s_%d.vtkhdf", basename, basename, iter);

  static int timeseries_len = 0;
  static char timeseries[1 << 16]; // 64 MB buffer

  // Tail that always lives *after* the last entry
  static const char timeseries_tail[] = ",\n\t]}";
  const int timeseries_delta = (int) strlen (timeseries_tail);

  char fname_timeseries[128];
  sprintf (
    fname_timeseries, "%s.vtkhdf.series", basename); // assuming fname exists

  if (!timeseries_len) {
#ifdef OUTPUT_DUMP_H
    if (checkpoint_state.restored) {
      // Read the existing timeseries file into our buffer
      FILE* fp = fopen (fname_timeseries, "r");

      char line[8192];
      size_t line_no = 1;
      bool found = false;

      while (fgets (line, sizeof line, fp) && !found) {
        if (line_no <= 3) {
          timeseries_len += sprintf (timeseries + timeseries_len, "%s", line);
        } else {
          int iter_val;
          // matches: \t\t{ "iter" : 268, ...
          if (sscanf (line, "\t\t{ \"iter\" : %d,", &iter_val) == 1) {
            if (iter_val >= checkpoint_state.restored_i) {
              found = true;
            }
          }
          // keep the line
          timeseries_len += sprintf (timeseries + timeseries_len, "%s", line);
        }
        line_no++;
      }
      timeseries_len +=
        sprintf (timeseries + timeseries_len - strlen(",\n"), "%s", timeseries_tail);
      fclose (fp);

      // We currently have ... { last entry }",\n\t]}"
      // The tail begins at timeseries_len - timeseries_delta
      int start = timeseries_len - timeseries_delta;

      // Overwrite the tail with: ",\n\t\t{ new entry }<tail again>"
      // This keeps all old content up to 'start', and appends a new entry
      timeseries_len =
        start +
        sprintf (timeseries + start,
                 "\t\t{ \"iter\" : %d, \"time\" : %f, \"name\" : \"%s\" }%s",
                 iter,
                 time,
                 fname,
                 timeseries_tail);

    } else {
      // First time: write full JSON with one entry and the tail
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
    }
#else
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
#endif
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

trace void output_hdf_htg (scalar* scalar_list = all,
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

  char basename_buff[32];

  if (!basename) {
    snprintf (
      basename_buff, sizeof (basename_buff), "%s", get_executable_name ());
  } else {
    snprintf (basename_buff, sizeof (basename_buff), "%s", basename);
  }

  output_hdf_htg_series (
    scalar_list, vector_list, basename_buff, iter, time, overwrite);

  // if (use_ab && use_transient) {

  //   static bool first_iter = true;

  //   char fname_a[64];
  //   sprintf (fname_a, "%s_a.vtkhdf", basename_buff);
  //   char fname_b[64];
  //   sprintf (fname_b, "%s_b.vtkhdf", basename_buff);

  //   if (first_iter) {
  //     first_iter = false;
  //     vtkHDFHyperTreeGrid vtk_hdf = vtk_HDF_hypertreegrid_init_transient (
  //       scalar_list, vector_list, time, fname_a, overwrite);
  //     vtk_HDF_hypertreegrid_close (&vtk_hdf);
  //     if (pid () == 0)
  //       copy_file (fname_a, fname_b);
  //   } else {
  //     vtkHDFHyperTreeGrid vtk_hdf = vtk_HDF_hypertreegrid_append_transient (
  //       scalar_list, vector_list, time, fname_a);
  //     vtk_HDF_hypertreegrid_close (&vtk_hdf);
  //     if (pid () == 0)
  //       copy_file (fname_a, fname_b);
  //   }
  // } else if (use_transient) {

  //   static bool first_iter = true;

  //   char fname[64];
  //   sprintf (fname, "%s.vtkhdf", basename_buff);
  //   if (first_iter && !checkpoint_state.restored) {
  //     first_iter = false;
  //     vtkHDFHyperTreeGrid vtk_hdf = vtk_HDF_hypertreegrid_init_transient (
  //       scalar_list, vector_list, time, fname, overwrite);
  //     vtk_HDF_hypertreegrid_close (&vtk_hdf);
  //   } else {
  //     vtkHDFHyperTreeGrid vtk_hdf = vtk_HDF_hypertreegrid_append_transient (
  //       scalar_list, vector_list, time, fname);
  //     vtk_HDF_hypertreegrid_close (&vtk_hdf);
  //   }
  // } else {
  //   output_hdf_htg_series (
  //     scalar_list, vector_list, basename_buff, iter, time, overwrite);
  // }
}
