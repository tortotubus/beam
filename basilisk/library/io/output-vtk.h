
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
      basename_buff, sizeof (basename_buff), "%s-hdf", get_executable_name ());
  } else {
    snprintf (basename_buff, sizeof (basename_buff), "%s-hdf", basename);
  }

  output_hdf_htg_series (
    scalar_list, vector_list, basename_buff, iter, time, overwrite);
}

#include "library/io/vtk/vtkHDFPolyData.h"
#include "library/io/vtk/vtkPolyData.h"
#include "library/ibm/IBMeshManager.h"

trace void output_hdf_pd_series (IBscalar* slist = all,
                                 IBvector* vlist = NULL,
                                 const char* basename,
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
  }

  // Init polydata

#if _MPI
  size_t n_points_local = 0;
  foreach_ibnode () {
    ibnode_update_pid (node);
    if (node->pid == pid ()) {
      n_points_local++;
    }
  }
#else
  size_t n_points_local = 0;
  foreach_ibnode () {
    n_points_local++;
    node->pid = 0.;
  }
#endif

  size_t n_points = n_points_local;
  size_t n_vertices = n_points;
  size_t n_lines = 0;
  size_t n_strips = 0;
  size_t n_polygons = 0;

  int n_scalars = 0, n_vectors = 0;
  foreach_ibscalar (slist) {
    n_scalars++;
  }
  foreach_ibvector (vlist) {
    n_vectors++;
  }
  size_t n_pointdata = n_scalars + n_vectors + 3;

  vtkPolyData vtk_pd = vtk_polydata_init (
    n_points, n_vertices, n_lines, n_strips, n_polygons, n_pointdata);

  // Populate polydata
  foreach_ibnode () {
    if (node->pid == pid ()) {
      coord pos = node->pos;
      vtk_polydata_add_point (&vtk_pd, pos.x, pos.y, pos.z);
    }
  }

  for (size_t i = 0; i < n_points; i++) {
    vtk_polydata_add_vertex (&vtk_pd, i);
  }

  size_t ncomp = dimension;

  double** scalar_data = NULL;
  double** vector_data = NULL;
  int64_t* scalar_ids = NULL;
  int64_t* vector_ids = NULL;
  int64_t pid_id = vtk_polydata_add_pointdata_scalar (&vtk_pd, "pid");
  double* pid_data = vtk_polydata_get_pointdata_data (&vtk_pd, pid_id);
  int64_t vel_id = vtk_polydata_add_pointdata_vector (&vtk_pd, "vel", dimension);
  double* vel_data = vtk_polydata_get_pointdata_data (&vtk_pd, vel_id);
  int64_t f_id = vtk_polydata_add_pointdata_vector (&vtk_pd, "f", dimension);
  double* f_data = vtk_polydata_get_pointdata_data (&vtk_pd, f_id);

  scalar_ids = malloc (sizeof (int64_t) * n_scalars);
  scalar_data = malloc (sizeof (double*) * n_scalars);
  vector_ids = malloc (sizeof (int64_t) * n_vectors);
  vector_data = malloc (sizeof (double*) * n_vectors);

  {
    int i;

    i = 0;
    foreach_ibscalar (slist) {
      scalar_ids[i] = vtk_polydata_add_pointdata_scalar (&vtk_pd, ibname (s));
      scalar_data[i] = vtk_polydata_get_pointdata_data (&vtk_pd, scalar_ids[i]);
      i++;
    }

    i = 0;
    foreach_ibvector (vlist) {
      char* vector_name;
      size_t trunc_len = (size_t) (strlen (ibname (v.x)) - 2);
      vector_name = malloc ((trunc_len + 1) * sizeof (char));
      strncpy (vector_name, ibname (v.x), trunc_len);
      vector_name[trunc_len] = '\0';
      vector_ids[i] = vtk_polydata_add_pointdata_scalar (&vtk_pd, vector_name);
      vector_data[i] = vtk_polydata_get_pointdata_data (&vtk_pd, vector_ids[i]);
      free (vector_name);
      i++;
    }

    i = 0;
    foreach_ibnode () {
      if (node->pid == pid ()) {

        // scalar ibfields
        size_t si = 0;
        foreach_ibscalar (slist) {
          scalar_data[si][i] = ibval (s);
          si++;
        }

        // vector ibfields
        size_t vi = 0;
        foreach_ibvector (vlist) {
          vector_data[vi][i * dimension + 0] = ibval (v.x);
#if dimension >= 2
          vector_data[vi][i * dimension + 1] = ibval (v.y);
#endif
#if dimension >= 3
          vector_data[vi][i * dimension + 2] = ibval (v.z);
#endif
          vi++;
        }

        // ib nodal values

        pid_data[i] = node->pid;
        vel_data[i * dimension + 0] = node->vel.x;
        f_data[i * dimension + 0] = node->f.x;
#if dimension >= 2
        vel_data[i * dimension + 1] = node->vel.y;
        f_data[i * dimension + 1] = node->f.y;
#endif
#if dimension >= 3
        vel_data[i * dimension + 2] = node->vel.z;
        f_data[i * dimension + 2] = node->f.z;
#endif

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

trace void output_hdf_pd (IBscalar* slist = iball,
                          IBvector* vlist = NULL,
                          const char* basename = NULL,
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

  output_hdf_pd_series (slist, vlist, basename_buff, iter, time, overwrite);
}

#include "library/io/vtk/vtkHDFImageData.h"

trace void output_hdf_imagedata_series (scalar* scalar_list,
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

  vtkHDFImageData vtk_hdf_imagedata =
    vtk_HDF_imagedata_init_static (scalar_list, vector_list, fname, overwrite);

  vtk_HDF_imagedata_close (&vtk_hdf_imagedata);
}

trace void output_hdf_imagedata (scalar* slist = all,
                                 vector* vlist = NULL,
                                 const char* basename = NULL,
                                 int iter = i,
                                 double time = t,
                                 bool use_transient = false,
                                 bool overwrite = true,
                                 bool use_ab = false) {
  char basename_buff[32];

  if (!basename) {
    snprintf (
      basename_buff, sizeof (basename_buff), "%s-id", get_executable_name ());
  } else {
    snprintf (basename_buff, sizeof (basename_buff), "%s-id", basename);
  }

  output_hdf_imagedata_series (
    slist, vlist, basename_buff, iter, time, overwrite);

  // char fname[128];
  // sprintf (fname, "%s-id_%i.vtkhdf", basename_buff, iter);

  // vtkHDFImageData vtk_hdf_imagedata =
  //   vtk_HDF_imagedata_init_static (slist, vlist, fname, overwrite);

  // vtk_HDF_imagedata_close (&vtk_hdf_imagedata);
}
