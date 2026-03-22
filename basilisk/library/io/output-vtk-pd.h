#include "library/io/output-common.h"

#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#include "library/io/vtk/vtkHDFPolyData.h"
#include "library/io/vtk/vtkPolyData.h"
#include "library/io/vtk/vtkTimeSeries.h"
#include "library/ibm/IBMeshManager.h"

trace void output_hdf_pd_series (IBscalar* scalar_list = NULL,
                                 IBvector* vector_list = NULL,
                                 const char* basename,
                                 int iter,
                                 double time,
                                 bool overwrite) {

  char pname[128];
  sprintf (pname, "%s/vtkhdf-pd-series", basename);

  char fname[128];
  sprintf (fname, "%s/vtkhdf-pd-series/polydata_%d.vtkhdf", basename, iter);

  char series_entry[128];
  sprintf (series_entry, "vtkhdf-pd-series/polydata_%d.vtkhdf", iter);

  char series_filename[128];
  sprintf (series_filename, "%s/pd.vtkhdf.series", basename);

  if (pid () == 0) {
    assert (!create_path (basename));
    assert (!create_path (pname));
    output_vtk_series (series_entry, series_filename, iter, time);
  }

IBscalar* slist = scalar_list;
IBvector* vlist = vector_list;

if (!slist && !vlist) {
  foreach_ibscalar (iball) {
    if (_ibattribute[s.i].v.x.i == s.i)
      vlist = ibvectors_add (vlist, _ibattribute[s.i].v);
    else if (_ibattribute[s.i].v.x.i < 0)
      slist = iblist_add (slist, s);
  }
}


  // Init polydata

#if _MPI
  ibmeshmanager_update_pid ();
  size_t n_points_local = 0;

  foreach_ibnode () {
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

  int n_scalars = iblist_len(slist);
  int n_vectors = ibvectors_len(vlist);
  size_t n_pointdata = n_scalars + n_vectors;

  vtkPolyData vtk_pd = vtk_polydata_init (
    n_points, n_vertices, n_lines, n_strips, n_polygons, n_pointdata);

  // Populate polydata
  foreach_ibnode () {
    if (node->pid == pid ()) {
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
      vector_ids[i] = vtk_polydata_add_pointdata_vector (&vtk_pd, vector_name, dimension);
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

trace void output_hdf_pd (IBscalar* scalar_list = NULL,
                          IBvector* vector_list = NULL,
                          const char* basename = NULL,
                          int iter = i,
                          double time = t,
                          bool use_transient = false,
                          bool overwrite = true,
                          bool use_ab = false) {

  char basename_buff[64];

  if (!basename) {
    snprintf (
      basename_buff, sizeof (basename_buff), "%s", get_unique_basename ());
  } else {
    snprintf (basename_buff, sizeof (basename_buff), "%s", basename);
  }

  output_hdf_pd_series (scalar_list, vector_list, basename_buff, iter, time, overwrite);
}
