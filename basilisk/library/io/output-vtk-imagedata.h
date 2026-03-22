#include "library/io/output-common.h"

#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#include "library/io/vtk/vtkHDFImageData.h"
#include "library/io/vtk/vtkTimeSeries.h"

trace void output_hdf_imagedata_series (scalar* scalar_list,
                                        vector* vector_list,
                                        const char* basename,
                                        int iter,
                                        double time,
                                        bool overwrite) {
  char pname[128];
  sprintf (pname, "%s/vtkhdf-imagedata-series", basename);

  char fname[128];
  sprintf (fname, "%s/vtkhdf-imagedata-series/imagedata_%d.vtkhdf", basename, iter);

  char series_entry[128];
  sprintf (series_entry, "vtkhdf-imagedata-series/imagedata_%d.vtkhdf", iter);

  char series_filename[128];
  sprintf (series_filename, "%s/imagedata.vtkhdf.series", basename);

  if (pid () == 0) {
    assert (!create_path (basename));
    assert (!create_path (pname));
    output_vtk_series (series_entry, series_filename, iter, time);
  }

  vtkHDFImageData vtk_hdf_imagedata =
    vtk_HDF_imagedata_init_static (scalar_list, vector_list, fname, overwrite);

  vtk_HDF_imagedata_close (&vtk_hdf_imagedata);
}

trace void output_hdf_imagedata (scalar* scalar_list = NULL,
                                 vector* vector_list = NULL,
                                 const char* basename = NULL,
                                 int iter = i,
                                 double time = t,
                                 bool use_transient = false,
                                 bool overwrite = true,
                                 bool use_ab = false) {

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
      basename_buff, sizeof (basename_buff), "%s", get_unique_basename ());
  } else {
    snprintf (basename_buff, sizeof (basename_buff), "%s", basename);
  }

  output_hdf_imagedata_series (
    slist, vlist, basename_buff, iter, time, overwrite);
}
