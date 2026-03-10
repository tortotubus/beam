#include "vtkHDF.h"

typedef struct {
  vtkHDF vtk_hdf;
  hid_t grp_fielddata_id;
  hid_t grp_pointdata_id;
  hid_t grp_celldata_id;
} vtkHDFImageData;

void vtk_HDF_imagedata_close (vtkHDFImageData* vtk_hdf_imagedata) {
  if (vtk_hdf_imagedata->grp_celldata_id >= 0 &&
      H5Iis_valid (vtk_hdf_imagedata->grp_celldata_id) > 0) {
    H5Gclose (vtk_hdf_imagedata->grp_celldata_id);
    vtk_hdf_imagedata->grp_celldata_id = H5I_INVALID_HID;
  }

  if (vtk_hdf_imagedata->grp_fielddata_id >= 0 &&
      H5Iis_valid (vtk_hdf_imagedata->grp_fielddata_id) > 0) {
    H5Gclose (vtk_hdf_imagedata->grp_fielddata_id);
    vtk_hdf_imagedata->grp_fielddata_id = H5I_INVALID_HID;
  }

  if (vtk_hdf_imagedata->grp_pointdata_id >= 0 &&
      H5Iis_valid (vtk_hdf_imagedata->grp_pointdata_id) > 0) {
    H5Gclose (vtk_hdf_imagedata->grp_pointdata_id);
    vtk_hdf_imagedata->grp_pointdata_id = H5I_INVALID_HID;
  }

  vtk_HDF_close (&vtk_hdf_imagedata->vtk_hdf);
}

vtkHDFImageData vtk_HDF_imagedata_init_static (scalar* scalar_list,
                                               vector* vector_list,
                                               const char* fname,
                                               bool overwrite) {
#if _MPI
  vtkHDF vtk_hdf = vtk_HDF_init_MPIIO (fname, overwrite);
#else
  vtkHDF vtk_hdf = vtk_HDF_init (fname, overwrite);
#endif

  vtkHDFImageData vtk_hdf_imagedata = {.vtk_hdf = vtk_hdf,
                                       .grp_fielddata_id = H5I_INVALID_HID,
                                       .grp_pointdata_id = H5I_INVALID_HID,
                                       .grp_celldata_id = H5I_INVALID_HID};

  /*
   * Attribute: /VTKHDF/Type
   * Datatype: H5T_C_S1 / char
   * Dimension: {1}
   */
  {
    vtk_HDF_write_type_attribute ("ImageData",
                                  vtk_hdf_imagedata.vtk_hdf.grp_vtkhdf_id,
                                  &vtk_hdf_imagedata.vtk_hdf);
  }

  /*
   * Attribute: /VTKHDF/Direction
   * Datatype: H5T_IEEE_F64LE / double
   * Dimension: {9}
   */
  {
    // Value for the Direction attribute
    const double dir_value[] = {1, 0, 0, 0, 1, 0, 0, 0, 1};

    // Dimensions for the attribute
    const hsize_t dims_attr[1] = {9};

    vtk_HDF_write_attribute (
      "Direction",                             /* attribute_name */
      &dir_value,                              /* data */
      H5T_IEEE_F64LE,                          /* dtype_id */
      vtk_hdf_imagedata.vtk_hdf.grp_vtkhdf_id, /* group_id */
      dims_attr,                               /* dims */
      &vtk_hdf_imagedata.vtk_hdf /* vtkHDFHyperTreeGrid object/struct */
    );
  }

  /*
   * Attribute: /VTKHDF/Origin
   * Datatype: H5T_IEEE_F64LE / double
   * Dimension: {3}
   */
  {
    // Value for the Origin attribute
#if dimension == 1
    const double origin_value[] = {X0, 0., 0.};
#elif dimension == 2
    const double origin_value[] = {X0, Y0, 0.};
#else // dimension == 3
    const double origin_value[] = {X0, Y0, Z0};
#endif

    // Dimensions for the attribute
    const hsize_t dims_attr[1] = {3};

    vtk_HDF_write_attribute (
      "Origin",                                /* attribute_name */
      &origin_value,                           /* data */
      H5T_IEEE_F64LE,                          /* dtype_id */
      vtk_hdf_imagedata.vtk_hdf.grp_vtkhdf_id, /* group_id */
      dims_attr,                               /* dims */
      &vtk_hdf_imagedata.vtk_hdf /* vtkHDFHyperTreeGrid object/struct */
    );
  }

  /*
   * Attribute: /VTKHDF/Spacing
   * Datatype: H5T_IEEE_F64LE / double
   * Dimension: {3}
   */
  {
    const double Delta = L0 / N;

#if dimension == 1
    const double spacing_value[] = {Delta, 0., 0.};
#elif dimension == 2
    const double spacing_value[] = {Delta, Delta, 0.};
#else
    const double spacing_value[] = {Delta, Delta, Delta};
#endif

    // Dimensions for the attribute
    const hsize_t dims_attr[1] = {3};

    vtk_HDF_write_attribute (
      "Spacing",                               /* attribute_name */
      &spacing_value,                          /* data */
      H5T_IEEE_F64LE,                          /* dtype_id */
      vtk_hdf_imagedata.vtk_hdf.grp_vtkhdf_id, /* group_id */
      dims_attr,                               /* dims */
      &vtk_hdf_imagedata.vtk_hdf /* vtkHDFHyperTreeGrid object/struct */
    );
  }

  /*
   * Attribute: /VTKHDF/Version
   * Datatype: H5T_NATIVE_INT64 / int64_t
   * Dimension: {2}
   */
  {
    // Value of the VTKHDF version attribute
    const int64_t vers_value[2] = {2, 4};

    // Dimensions of the attribute
    const hsize_t dims_attr[1] = {2};

    vtk_HDF_write_attribute (
      "Version",                               /* attribute_name */
      vers_value,                              /* data */
      H5T_NATIVE_INT64,                        /* dtype_id */
      vtk_hdf_imagedata.vtk_hdf.grp_vtkhdf_id, /* group_id */
      dims_attr,                               /* dims */
      &vtk_hdf_imagedata.vtk_hdf               /* vtkHDF object/struct */
    );
  }

  /*
   * Attribute: /VTKHDF/WholeExtent
   * Datatype: H5T_NATIVE_INT64 / int64_t
   * Dimension: {6}
   */
  {
    // Value of the VTKHDF version attribute

#if dimension == 1
    const int64_t whole_extent_val[] = {
      0, (int64_t) ((1 << depth ()) * Dimensions.x), 0, 0, 0, 0};
#elif dimension == 2
    const int64_t whole_extent_val[] = {
      0,
      (int64_t) ((1 << depth ()) * Dimensions.x),
      0,
      (int64_t) ((1 << depth ()) * Dimensions.y),
      0,
      0};
#else // dimension == 3
    const int64_t whole_extent_val[] = {
      0,
      (int64_t) ((1 << depth ()) * Dimensions.x),
      0,
      (int64_t) ((1 << depth ()) * Dimensions.y),
      0,
      (int64_t) ((1 << depth ()) * Dimensions.z)};
#endif

    // Dimensions of the attribute
    const hsize_t dims_attr[1] = {6};

    vtk_HDF_write_attribute (
      "WholeExtent",                           /* attribute_name */
      whole_extent_val,                        /* data */
      H5T_NATIVE_INT64,                        /* dtype_id */
      vtk_hdf_imagedata.vtk_hdf.grp_vtkhdf_id, /* group_id */
      dims_attr,                               /* dims */
      &vtk_hdf_imagedata.vtk_hdf               /* vtkHDF object/struct */
    );
  }

  // return vtk_hdf_imagedata;

  {
    // Create the new CellData group inside group VTKHDF
    vtk_hdf_imagedata.grp_celldata_id =
      H5Gcreate2 (vtk_hdf_imagedata.vtk_hdf.grp_vtkhdf_id,
                  "CellData",
                  H5P_DEFAULT,
                  H5P_DEFAULT,
                  H5P_DEFAULT);
    vtk_HDF_check_object (&vtk_hdf_imagedata.vtk_hdf,
                          vtk_hdf_imagedata.grp_celldata_id);
  }

  {
    for (scalar s in scalar_list) {

#if dimension == 1
      const size_t nx = (size_t) ((1 << depth ()) * Dimensions.x);
      const size_t n = nx;

#elif dimension == 2
      const size_t nx = (size_t) ((1 << depth ()) * Dimensions.x);
      const size_t ny = (size_t) ((1 << depth ()) * Dimensions.y);
      const size_t n = nx * ny;

#else // dimension == 3
      const size_t nx = (size_t) ((1 << depth ()) * Dimensions.x);
      const size_t ny = (size_t) ((1 << depth ()) * Dimensions.y);
      const size_t nz = (size_t) ((1 << depth ()) * Dimensions.z);
      const size_t n = nx * ny * nz;
#endif

      float* s_data = malloc (n * sizeof (float));
      if (!s_data) {
        perror ("malloc(s_data)");
        exit (1);
      }

      foreach (serial) {
#if dimension == 1
        int i = point.i - GHOSTS;

        size_t vtk_idx = (size_t) i;

#elif dimension == 2
        int i = point.i - GHOSTS;
        int j = point.j - GHOSTS;
        size_t vtk_idx = (size_t) i + (size_t) nx * (size_t) j;

#else
        int i = point.i - GHOSTS;
        int j = point.j - GHOSTS;
        int k = point.k - GHOSTS;
        size_t vtk_idx =
          (size_t) i + (size_t) nx * ((size_t) j + (size_t) ny * (size_t) k);
#endif

        s_data[vtk_idx] = (float) val (s);
      }

#if dimension == 1
      hsize_t dims[] = {nx};
      int rank = 1;
#elif dimension == 2
      hsize_t dims[] = {ny, nx};
      int rank = 2;
#else
      hsize_t dims[] = {nz, ny, nx};
      int rank = 3;
#endif

      vtk_HDF_write_dataset (s.name,
                             s_data,
                             H5T_IEEE_F32LE,
                             vtk_hdf_imagedata.grp_celldata_id,
                             rank,
                             dims,
                             &vtk_hdf_imagedata.vtk_hdf);

      free (s_data);
    }
  }

  // {
  //   // Create the new PointData group inside group VTKHDF
  //   vtk_hdf_imagedata.grp_pointdata_id =
  //     H5Gcreate2 (vtk_hdf_imagedata.vtk_hdf.grp_vtkhdf_id,
  //                 "PointData",
  //                 H5P_DEFAULT,
  //                 H5P_DEFAULT,
  //                 H5P_DEFAULT);
  //   vtk_HDF_check_object (&vtk_hdf_imagedata.vtk_hdf,
  //                         vtk_hdf_imagedata.grp_pointdata_id);
  // }

  // {
  //   // Create the new FieldData group inside group VTKHDF
  //   vtk_hdf_imagedata.grp_fielddata_id =
  //     H5Gcreate2 (vtk_hdf_imagedata.vtk_hdf.grp_vtkhdf_id,
  //                 "FieldData",
  //                 H5P_DEFAULT,
  //                 H5P_DEFAULT,
  //                 H5P_DEFAULT);
  //   vtk_HDF_check_object (&vtk_hdf_imagedata.vtk_hdf,
  //                         vtk_hdf_imagedata.grp_fielddata_id);
  // }

  return vtk_hdf_imagedata;
}
