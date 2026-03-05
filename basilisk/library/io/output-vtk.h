
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
  static char timeseries[1 << 16]; // 64 MB buffer

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
}

#include "library/io/vtk/vtkHDFPolyData.h"
#include "library/io/vtk/vtkPolyData.h"
#include "library/ibm/IBMeshManager.h"

trace void output_hdf_pd_series (const char* basename,
                                 int iter,
                                 double time,
                                 bool overwrite) {

  char fname[128];
  sprintf (fname, "%s-vtkhdf-series/%s_%d.vtkhdf", basename, basename, iter);

  if (pid () == 0) {
    char pname[128];
    sprintf (pname, "%s-vtkhdf-series", basename);
    create_path (pname);

    static int timeseries_len = 0;
    static char timeseries[1 << 16]; // 64 MB buffer

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
  }

  // Init polydata

  size_t n_points_local = 0;
  foreach_ibnode () {
    ibnode_update_pid (node);
    if (node->pid == pid ())
      n_points_local++;
  }

  size_t n_points = n_points_local;
  size_t n_vertices = n_points;
  size_t n_lines = 0;
  size_t n_strips = 0;
  size_t n_polygons = 0;
  size_t n_pointdata = 4;

  vtkPolyData vtk_pd = vtk_polydata_init (
    n_points, n_vertices, n_lines, n_strips, n_polygons, n_pointdata);

  // Populate polydata
  foreach_ibnode () {
    if (node->pid == pid ()) {
      coord pos = node->lagpos;
      vtk_polydata_add_point (&vtk_pd, pos.x, pos.y, pos.z);
    }
  }

  for (size_t i = 0; i < n_points; i++) {
    vtk_polydata_add_vertex (&vtk_pd, i);
  }

  size_t ncomp = dimension;

  int64_t force_id =
    vtk_polydata_add_pointdata_vector (&vtk_pd, "force", ncomp);
  double* force_data = vtk_polydata_get_pointdata_data (&vtk_pd, force_id);

  int64_t lagvel_id =
    vtk_polydata_add_pointdata_vector (&vtk_pd, "lagvel", ncomp);
  double* lagvel_data = vtk_polydata_get_pointdata_data (&vtk_pd, lagvel_id);

  int64_t eulvel_id =
    vtk_polydata_add_pointdata_vector (&vtk_pd, "eulvel", ncomp);
  double* eulvel_data = vtk_polydata_get_pointdata_data (&vtk_pd, eulvel_id);

  int64_t pid_id = vtk_polydata_add_pointdata_scalar (&vtk_pd, "node_pid");
  double* pid_data = vtk_polydata_get_pointdata_data (&vtk_pd, pid_id);

  {
    size_t i = 0;
    foreach_ibnode () {
      if (node->pid == pid ()) {
#if dimension >= 1
        force_data[i * dimension + 0] = node->force.x;
        lagvel_data[i * dimension + 0] = node->lagvel.x;
        eulvel_data[i * dimension + 0] = node->eulvel.x;
#endif
#if dimension >= 2
        force_data[i * dimension + 1] = node->force.y;
        lagvel_data[i * dimension + 1] = node->lagvel.y;
        eulvel_data[i * dimension + 1] = node->eulvel.y;
#endif
#if dimension >= 3
        force_data[i * dimension + 2] = node->force.z;
        lagvel_data[i * dimension + 2] = node->lagvel.z;
        eulvel_data[i * dimension + 2] = node->eulvel.z;
#endif
        pid_data[i] = node->pid;

        i++;
      }
    }
  }

  // Write vtkhdf polydata
  vtkHDFPolyData vtk_hdf_pd =
    vtk_HDF_polydata_init_static (fname, true, &vtk_pd);
  vtk_HDF_polydata_close (&vtk_hdf_pd);

  // Free vtkPolyData
  vtk_polydata_free (&vtk_pd);
}

trace void output_hdf_pd (const char* basename = NULL,
                          int iter = i,
                          double time = t,
                          bool use_transient = false,
                          bool overwrite = true,
                          bool use_ab = false) {
  char basename_buff[32];

  if (!basename) {
    snprintf (
      basename_buff, sizeof (basename_buff), "%s-pd", get_executable_name ());
  } else {
    snprintf (basename_buff, sizeof (basename_buff), "%s-pd", basename);
  }

  output_hdf_pd_series (basename_buff, iter, time, overwrite);
}