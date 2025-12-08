#include "vtkHDFPolyData.hpp"

namespace beam {
namespace io {
namespace C {

void
vtk_HDF_polydata_close(vtkHDFPolyData* vtk_hdf_pd)
{
  if (vtk_hdf_pd->grp_celldata_id >= 0)
    H5Gclose(vtk_hdf_pd->grp_celldata_id);

  if (vtk_hdf_pd->grp_vertices_id >= 0)
    H5Gclose(vtk_hdf_pd->grp_vertices_id);

  if (vtk_hdf_pd->grp_lines_id >= 0)
    H5Gclose(vtk_hdf_pd->grp_lines_id);

  if (vtk_hdf_pd->grp_strips_id >= 0)
    H5Gclose(vtk_hdf_pd->grp_strips_id);

  if (vtk_hdf_pd->grp_polygons_id >= 0)
    H5Gclose(vtk_hdf_pd->grp_polygons_id);

  if (vtk_hdf_pd->grp_steps_id >= 0)
    H5Gclose(vtk_hdf_pd->grp_steps_id);

  if (vtk_hdf_pd->grp_celldata_offsets_id >= 0)
    H5Gclose(vtk_hdf_pd->grp_celldata_offsets_id);

  if (vtk_hdf_pd->grp_pointdata_offsets_id >= 0)
    H5Gclose(vtk_hdf_pd->grp_pointdata_offsets_id);

  vtk_HDF_close(&vtk_hdf_pd->vtk_hdf);
}

vtkHDFPolyData
vtk_HDF_polydata_init_static(const char* fname,
                             bool overwrite,
                             vtkPolyData* vtk_pd)
{
  vtkHDFPolyData vtk_hdf_pd = {
    .grp_vertices_id = H5I_INVALID_HID,         /**< Vertices group */
    .grp_lines_id = H5I_INVALID_HID,            /**< Lines group */
    .grp_strips_id = H5I_INVALID_HID,           /**< Strips group */
    .grp_polygons_id = H5I_INVALID_HID,         /**< Polygons group */
    .grp_celldata_id = H5I_INVALID_HID,         /**< CellData group */
    .grp_pointdata_id = H5I_INVALID_HID,        /**< PointData group */
    .grp_steps_id = H5I_INVALID_HID,            /**< Steps group */
    .grp_celldata_offsets_id = H5I_INVALID_HID, /**< CellData Offsets group */
    .grp_pointdata_offsets_id = H5I_INVALID_HID /**< PointData Offsets group */
  };

  typedef struct
  {
    vtkHDF vtk_hdf;                 /**< Parent */
    hid_t grp_vertices_id;          /**< /VTKHDF/Vertices group */
    hid_t grp_lines_id;             /**< /VTKHDF/Lines group */
    hid_t grp_strips_id;            /**< /VTKHDF/Strips group */
    hid_t grp_polygons_id;          /**< /VTKHDF/Polygons group */
    hid_t grp_celldata_id;          /**< /VTKHDF/CellData group */
    hid_t grp_pointdata_id;         /**< /VTKHDF/PointData group */
    hid_t grp_steps_id;             /**< /VTKHDF/Steps group */
    hid_t grp_celldata_offsets_id;  /**< /VTKHDF/Steps/CellDataOffsets group */
    hid_t grp_pointdata_offsets_id; /**< /VTKHDF/Steps/PointDataOffsets group */
  } vtkHDFPolyData;

  vtk_hdf_pd.vtk_hdf = vtk_HDF_init(fname, overwrite);

  vtk_hdf_pd.grp_celldata_id = H5Gcreate2(vtk_hdf_pd.vtk_hdf.file_id,
                                          "/VTKHDF/CellData",
                                          H5P_DEFAULT,
                                          H5P_DEFAULT,
                                          H5P_DEFAULT);

  vtk_HDF_check_object(&vtk_hdf_pd.vtk_hdf, vtk_hdf_pd.grp_celldata_id);

  vtk_hdf_pd.grp_pointdata_id = H5Gcreate2(vtk_hdf_pd.vtk_hdf.file_id,
                                           "/VTKHDF/PointData",
                                           H5P_DEFAULT,
                                           H5P_DEFAULT,
                                           H5P_DEFAULT);

  vtk_HDF_check_object(&vtk_hdf_pd.vtk_hdf, vtk_hdf_pd.grp_pointdata_id);

  vtk_hdf_pd.grp_vertices_id = H5Gcreate2(vtk_hdf_pd.vtk_hdf.file_id,
                                          "/VTKHDF/Vertices",
                                          H5P_DEFAULT,
                                          H5P_DEFAULT,
                                          H5P_DEFAULT);
  vtk_HDF_check_object(&vtk_hdf_pd.vtk_hdf, vtk_hdf_pd.grp_vertices_id);

  vtk_hdf_pd.grp_lines_id = H5Gcreate2(vtk_hdf_pd.vtk_hdf.file_id,
                                       "/VTKHDF/Lines",
                                       H5P_DEFAULT,
                                       H5P_DEFAULT,
                                       H5P_DEFAULT);

  vtk_HDF_check_object(&vtk_hdf_pd.vtk_hdf, vtk_hdf_pd.grp_lines_id);

  vtk_hdf_pd.grp_strips_id = H5Gcreate2(vtk_hdf_pd.vtk_hdf.file_id,
                                        "/VTKHDF/Strips",
                                        H5P_DEFAULT,
                                        H5P_DEFAULT,
                                        H5P_DEFAULT);

  vtk_HDF_check_object(&vtk_hdf_pd.vtk_hdf, vtk_hdf_pd.grp_lines_id);

  vtk_hdf_pd.grp_polygons_id = H5Gcreate2(vtk_hdf_pd.vtk_hdf.file_id,
                                          "/VTKHDF/Polygons",
                                          H5P_DEFAULT,
                                          H5P_DEFAULT,
                                          H5P_DEFAULT);

  vtk_HDF_check_object(&vtk_hdf_pd.vtk_hdf, vtk_hdf_pd.grp_polygons_id);

  /*
   * Attribute: /VTKHDF/Version
   */
  {
    int64_t attribute_data[] = { 2, 4 };
    hsize_t attribute_dims[] = { 2 };
    const char* attribute_name = "Version";
    hid_t attribute_datatype = H5T_STD_I64LE;
    hid_t attribute_group = vtk_hdf_pd.vtk_hdf.grp_vtkhdf_id;

    vtk_HDF_write_attribute(attribute_name,
                            attribute_data,
                            attribute_datatype,
                            attribute_group,
                            attribute_dims,
                            &vtk_hdf_pd.vtk_hdf);
  }

  /*
   * Attribute: /VTKHDF/Type
   */
  {
    const char type_name[] = "PolyData";
    hid_t attribute_group = vtk_hdf_pd.vtk_hdf.grp_vtkhdf_id;
    vtk_HDF_write_type_attribute(
      type_name, attribute_group, &vtk_hdf_pd.vtk_hdf);
  }

  /*
   * Dataset: /VTKHDF/Points
   */
  {
    float* dataset_data = vtk_pd->points;
    hsize_t dataset_dims[] = { 0, 3 };
    dataset_dims[0] = vtk_polydata_number_of_points(vtk_pd);

    const char* dataset_name = "Points";
    hid_t dataset_datatype = H5T_IEEE_F32LE;
    hid_t dataset_group = vtk_hdf_pd.vtk_hdf.grp_vtkhdf_id;

    vtk_HDF_write_dataset(dataset_name,
                          dataset_data,
                          dataset_datatype,
                          dataset_group,
                          2,
                          dataset_dims,
                          &vtk_hdf_pd.vtk_hdf);
  }

  /*
   * Dataset: /VTKHDF/NumberOfPoints
   */
  {
    int64_t dataset_data[] = { 0 };
    dataset_data[0] = vtk_polydata_number_of_points(vtk_pd);
    hsize_t dataset_dims[] = { 1 };
    const char* dataset_name = "NumberOfPoints";
    hid_t dataset_datatype = H5T_STD_I64LE;
    hid_t dataset_group = vtk_hdf_pd.vtk_hdf.grp_vtkhdf_id;

    vtk_HDF_write_dataset(dataset_name,
                          dataset_data,
                          dataset_datatype,
                          dataset_group,
                          1,
                          dataset_dims,
                          &vtk_hdf_pd.vtk_hdf);
  }

  {
    /*
     * Dataset: /VTKHDF/Lines/Connectivity
     */
    {
      int64_t* dataset_data = vtk_pd->lines_connectivity;
      hsize_t dataset_dims[] = { vtk_pd->n_lines_connectivity };

      const char* dataset_name = "Connectivity";
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_lines_id;

      vtk_HDF_write_dataset(dataset_name,
                            dataset_data,
                            dataset_datatype,
                            dataset_group,
                            1,
                            dataset_dims,
                            &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Lines/NumberOfCells
     */
    {
      int64_t dataset_data[] = { 0 };
      dataset_data[0] = vtk_polydata_number_of_lines(vtk_pd);
      hsize_t dataset_dims[] = { 1 };
      const char* dataset_name = "NumberOfCells";
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_lines_id;

      vtk_HDF_write_dataset(dataset_name,
                            dataset_data,
                            dataset_datatype,
                            dataset_group,
                            1,
                            dataset_dims,
                            &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Lines/NumberOfConnectivityIds
     */
    {
      int64_t dataset_data[] = { (int64_t)vtk_pd->n_lines_connectivity };
      hsize_t dataset_dims[] = { 1 };
      const char* dataset_name = "NumberOfConnectivityIds";
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_lines_id;

      vtk_HDF_write_dataset(dataset_name,
                            dataset_data,
                            dataset_datatype,
                            dataset_group,
                            1,
                            dataset_dims,
                            &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Lines/Offsets
     */
    {
      int64_t* dataset_data = vtk_pd->lines_offsets;
      hsize_t dataset_dims[] = { vtk_pd->n_lines_offsets };
      const char* dataset_name = "Offsets";
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_lines_id;

      vtk_HDF_write_dataset(dataset_name,
                            dataset_data,
                            dataset_datatype,
                            dataset_group,
                            1,
                            dataset_dims,
                            &vtk_hdf_pd.vtk_hdf);
    }
  }
  {
    /*
     * Dataset: /VTKHDF/Polygons/Connectivity
     */
    {
      int64_t* dataset_data = vtk_pd->polygons_connectivity;
      hsize_t dataset_dims[] = { vtk_pd->n_polygons_connectivity };
      const char* dataset_name = "Connectivity";
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_polygons_id;

      vtk_HDF_write_dataset(dataset_name,
                            dataset_data,
                            dataset_datatype,
                            dataset_group,
                            1,
                            dataset_dims,
                            &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Polygons/NumberOfCells
     */
    {
      int64_t dataset_data[] = { 0 };
      dataset_data[0] = vtk_polydata_number_of_polygons(vtk_pd);
      hsize_t dataset_dims[] = { 1 };
      const char* dataset_name = "NumberOfCells";
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_polygons_id;

      vtk_HDF_write_dataset(dataset_name,
                            dataset_data,
                            dataset_datatype,
                            dataset_group,
                            1,
                            dataset_dims,
                            &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Polygons/NumberOfConnectivityIds
     */
    {
      int64_t dataset_data[] = { (int64_t)vtk_pd->n_polygons_connectivity };
      hsize_t dataset_dims[] = { 1 };
      const char* dataset_name = "NumberOfConnectivityIds";
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_polygons_id;

      vtk_HDF_write_dataset(dataset_name,
                            dataset_data,
                            dataset_datatype,
                            dataset_group,
                            1,
                            dataset_dims,
                            &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Polygons/Offsets
     */
    {
      int64_t* dataset_data = vtk_pd->polygons_offsets;
      hsize_t dataset_dims[] = { vtk_pd->n_polygons_offsets };
      const char* dataset_name = "Offsets";
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_polygons_id;

      vtk_HDF_write_dataset(dataset_name,
                            dataset_data,
                            dataset_datatype,
                            dataset_group,
                            1,
                            dataset_dims,
                            &vtk_hdf_pd.vtk_hdf);
    }
  }

  {
    /*
     * Dataset: /VTKHDF/Strips/Connectivity
     */
    {
      int64_t* dataset_data = vtk_pd->strips_connectivity;
      hsize_t dataset_dims[] = { vtk_pd->n_strips_connectivity };
      const char* dataset_name = "Connectivity";
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_strips_id;

      vtk_HDF_write_dataset(dataset_name,
                            dataset_data,
                            dataset_datatype,
                            dataset_group,
                            1,
                            dataset_dims,
                            &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Strips/NumberOfCells
     */
    {
      int64_t dataset_data[] = { 0 };
      dataset_data[0] = vtk_polydata_number_of_strips(vtk_pd);
      hsize_t dataset_dims[] = { 1 };
      const char* dataset_name = "NumberOfCells";
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_strips_id;

      vtk_HDF_write_dataset(dataset_name,
                            dataset_data,
                            dataset_datatype,
                            dataset_group,
                            1,
                            dataset_dims,
                            &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Strips/NumberOfConnectivityIds
     */
    {
      int64_t dataset_data[] = { (int64_t)vtk_pd->n_strips_connectivity };
      hsize_t dataset_dims[] = { 1 };
      const char* dataset_name = "NumberOfConnectivityIds";
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_strips_id;

      vtk_HDF_write_dataset(dataset_name,
                            dataset_data,
                            dataset_datatype,
                            dataset_group,
                            1,
                            dataset_dims,
                            &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Strips/Offsets
     */
    {
      int64_t* dataset_data = vtk_pd->strips_offsets;
      hsize_t dataset_dims[] = { 1 };
      const char* dataset_name = "Offsets";
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_strips_id;

      vtk_HDF_write_dataset(dataset_name,
                            dataset_data,
                            dataset_datatype,
                            dataset_group,
                            1,
                            dataset_dims,
                            &vtk_hdf_pd.vtk_hdf);
    }
  }

  {
    /*
     * Dataset: /VTKHDF/Vertices/Connectivity
     */
    {
      int64_t* dataset_data = vtk_pd->vertices_connectivity;
      hsize_t dataset_dims[] = { vtk_pd->n_vertices_connectivity };
      const char* dataset_name = "Connectivity";
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_vertices_id;

      vtk_HDF_write_dataset(dataset_name,
                            dataset_data,
                            dataset_datatype,
                            dataset_group,
                            1,
                            dataset_dims,
                            &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Vertices/NumberOfCells
     */
    {
      int64_t dataset_data[] = { 0 };
      dataset_data[0] = vtk_polydata_number_of_vertices(vtk_pd);
      hsize_t dataset_dims[] = { 1 };
      const char* dataset_name = "NumberOfCells";
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_vertices_id;

      vtk_HDF_write_dataset(dataset_name,
                            dataset_data,
                            dataset_datatype,
                            dataset_group,
                            1,
                            dataset_dims,
                            &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Vertices/NumberOfConnectivityIds
     */
    {
      int64_t dataset_data[] = { (int64_t)vtk_pd->n_vertices_connectivity };
      hsize_t dataset_dims[] = { 1 };
      const char* dataset_name = "NumberOfConnectivityIds";
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_vertices_id;

      vtk_HDF_write_dataset(dataset_name,
                            dataset_data,
                            dataset_datatype,
                            dataset_group,
                            1,
                            dataset_dims,
                            &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Vertices/Offsets
     */
    {
      int64_t* dataset_data = vtk_pd->vertices_offsets;
      hsize_t dataset_dims[] = { vtk_pd->n_vertices_offsets };
      const char* dataset_name = "Offsets";
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_vertices_id;

      vtk_HDF_write_dataset(dataset_name,
                            dataset_data,
                            dataset_datatype,
                            dataset_group,
                            1,
                            dataset_dims,
                            &vtk_hdf_pd.vtk_hdf);
    }
  }
  return vtk_hdf_pd;
}

vtkHDFPolyData
vtk_HDF_polydata_init_transient(const char* fname, bool overwrite, vtkPolyData* vtk_pd, float time)
{
  vtkHDFPolyData vtk_hdf_pd = {
    .grp_vertices_id = H5I_INVALID_HID,         /**< Vertices group */
    .grp_lines_id = H5I_INVALID_HID,            /**< Lines group */
    .grp_strips_id = H5I_INVALID_HID,           /**< Strips group */
    .grp_polygons_id = H5I_INVALID_HID,         /**< Polygons group */
    .grp_celldata_id = H5I_INVALID_HID,         /**< CellData group */
    .grp_pointdata_id = H5I_INVALID_HID,        /**< PointData group */
    .grp_steps_id = H5I_INVALID_HID,            /**< Steps group */
    .grp_celldata_offsets_id = H5I_INVALID_HID, /**< CellData Offsets group */
    .grp_pointdata_offsets_id = H5I_INVALID_HID /**< PointData Offsets group */
  };

  vtk_hdf_pd.vtk_hdf = vtk_HDF_init(fname, overwrite);

  vtk_hdf_pd.grp_celldata_id = H5Gcreate2(vtk_hdf_pd.vtk_hdf.file_id,
                                          "/VTKHDF/CellData",
                                          H5P_DEFAULT,
                                          H5P_DEFAULT,
                                          H5P_DEFAULT);

  vtk_HDF_check_object(&vtk_hdf_pd.vtk_hdf, vtk_hdf_pd.grp_celldata_id);

  vtk_hdf_pd.grp_pointdata_id = H5Gcreate2(vtk_hdf_pd.vtk_hdf.file_id,
                                           "/VTKHDF/PointData",
                                           H5P_DEFAULT,
                                           H5P_DEFAULT,
                                           H5P_DEFAULT);

  vtk_HDF_check_object(&vtk_hdf_pd.vtk_hdf, vtk_hdf_pd.grp_pointdata_id);

  vtk_hdf_pd.grp_vertices_id = H5Gcreate2(vtk_hdf_pd.vtk_hdf.file_id,
                                          "/VTKHDF/Vertices",
                                          H5P_DEFAULT,
                                          H5P_DEFAULT,
                                          H5P_DEFAULT);
  vtk_HDF_check_object(&vtk_hdf_pd.vtk_hdf, vtk_hdf_pd.grp_vertices_id);

  vtk_hdf_pd.grp_lines_id = H5Gcreate2(vtk_hdf_pd.vtk_hdf.file_id,
                                       "/VTKHDF/Lines",
                                       H5P_DEFAULT,
                                       H5P_DEFAULT,
                                       H5P_DEFAULT);

  vtk_HDF_check_object(&vtk_hdf_pd.vtk_hdf, vtk_hdf_pd.grp_lines_id);

  vtk_hdf_pd.grp_strips_id = H5Gcreate2(vtk_hdf_pd.vtk_hdf.file_id,
                                        "/VTKHDF/Strips",
                                        H5P_DEFAULT,
                                        H5P_DEFAULT,
                                        H5P_DEFAULT);

  vtk_HDF_check_object(&vtk_hdf_pd.vtk_hdf, vtk_hdf_pd.grp_lines_id);

  vtk_hdf_pd.grp_polygons_id = H5Gcreate2(vtk_hdf_pd.vtk_hdf.file_id,
                                          "/VTKHDF/Polygons",
                                          H5P_DEFAULT,
                                          H5P_DEFAULT,
                                          H5P_DEFAULT);

  vtk_HDF_check_object(&vtk_hdf_pd.vtk_hdf, vtk_hdf_pd.grp_polygons_id);

  vtk_hdf_pd.grp_steps_id = H5Gcreate2(vtk_hdf_pd.vtk_hdf.file_id,
                                       "/VTKHDF/Steps",
                                       H5P_DEFAULT,
                                       H5P_DEFAULT,
                                       H5P_DEFAULT);

  vtk_HDF_check_object(&vtk_hdf_pd.vtk_hdf, vtk_hdf_pd.grp_steps_id);

  vtk_hdf_pd.grp_celldata_offsets_id =
    H5Gcreate2(vtk_hdf_pd.vtk_hdf.file_id,
               "/VTKHDF/Steps/CellDataOffsets",
               H5P_DEFAULT,
               H5P_DEFAULT,
               H5P_DEFAULT);

  vtk_HDF_check_object(&vtk_hdf_pd.vtk_hdf, vtk_hdf_pd.grp_celldata_offsets_id);

  vtk_hdf_pd.grp_pointdata_offsets_id =
    H5Gcreate2(vtk_hdf_pd.vtk_hdf.file_id,
               "/VTKHDF/Steps/PointDataOffsets",
               H5P_DEFAULT,
               H5P_DEFAULT,
               H5P_DEFAULT);

  vtk_HDF_check_object(&vtk_hdf_pd.vtk_hdf,
                       vtk_hdf_pd.grp_pointdata_offsets_id);

  /*
   * Attribute: /VTKHDF/Version
   */
  {
    int64_t attribute_data[] = { 2, 4 };
    hsize_t attribute_dims[] = { 2 };
    const char* attribute_name = "Version";
    hid_t attribute_datatype = H5T_STD_I64LE;
    hid_t attribute_group = vtk_hdf_pd.vtk_hdf.grp_vtkhdf_id;

    vtk_HDF_write_attribute(attribute_name,
                            attribute_data,
                            attribute_datatype,
                            attribute_group,
                            attribute_dims,
                            &vtk_hdf_pd.vtk_hdf);
  }

  /*
   * Attribute: /VTKHDF/Type
   */
  {
    const char type_name[] = "PolyData";
    hid_t attribute_group = vtk_hdf_pd.vtk_hdf.grp_vtkhdf_id;
    vtk_HDF_write_type_attribute(
      type_name, attribute_group, &vtk_hdf_pd.vtk_hdf);
  }

  /*
   * Dataset: /VTKHDF/Points
   */
  {
    const char* dataset_name = "Points";
    float* dataset_data = vtk_pd->points;
    hid_t dataset_datatype = H5T_IEEE_F32LE;
    hid_t dataset_group = vtk_hdf_pd.vtk_hdf.grp_vtkhdf_id;
    hsize_t dataset_dims[] = { 0, 3 };
    dataset_dims[0] = vtk_polydata_number_of_points(vtk_pd);
    hsize_t dataset_max_dims[] = { H5S_UNLIMITED, 3 };
    hsize_t dataset_chunk_dims[] = { 4000, 3 };

    vtk_HDF_write_chunked_dataset(dataset_name,
                                  dataset_data,
                                  dataset_datatype,
                                  dataset_group,
                                  2,
                                  dataset_dims,
                                  dataset_max_dims,
                                  dataset_chunk_dims,
                                  &vtk_hdf_pd.vtk_hdf);
  }

  /*
   * Dataset: /VTKHDF/NumberOfPoints
   */
  {
    const char* dataset_name = "NumberOfPoints";
    int64_t dataset_data[] = { 0 };
    dataset_data[0] = vtk_polydata_number_of_points(vtk_pd);
    hid_t dataset_datatype = H5T_STD_I64LE;
    hid_t dataset_group = vtk_hdf_pd.vtk_hdf.grp_vtkhdf_id;
    hsize_t dataset_dims[] = { 1 };
    hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
    hsize_t dataset_chunk_dims[] = { 4000 };

    vtk_HDF_write_chunked_dataset(dataset_name,
                                  dataset_data,
                                  dataset_datatype,
                                  dataset_group,
                                  1,
                                  dataset_dims,
                                  dataset_max_dims,
                                  dataset_chunk_dims,
                                  &vtk_hdf_pd.vtk_hdf);
  }

  /*
   * Group: /VTKHDF/Lines
   */
  {
    /*
     * Dataset: /VTKHDF/Lines/Connectivity
     */
    {
      const char* dataset_name = "Connectivity";
      int64_t* dataset_data = vtk_pd->lines_connectivity;
      hsize_t dataset_dims[] = { vtk_pd->n_lines_connectivity };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_lines_id;

      vtk_HDF_write_chunked_dataset(dataset_name,
                                    dataset_data,
                                    dataset_datatype,
                                    dataset_group,
                                    1,
                                    dataset_dims,
                                    dataset_max_dims,
                                    dataset_chunk_dims,
                                    &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Lines/NumberOfCells
     */
    {
      const char* dataset_name = "NumberOfCells";
      int64_t dataset_data[] = { 0 };
      dataset_data[0] = vtk_polydata_number_of_lines(vtk_pd);
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_lines_id;
      hsize_t dataset_dims[] = { 1 };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_write_chunked_dataset(dataset_name,
                                    dataset_data,
                                    dataset_datatype,
                                    dataset_group,
                                    1,
                                    dataset_dims,
                                    dataset_max_dims,
                                    dataset_chunk_dims,
                                    &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Lines/NumberOfConnectivityIds
     */
    {
      const char* dataset_name = "NumberOfConnectivityIds";
      int64_t dataset_data[] = { (int64_t)vtk_pd->n_lines_connectivity };
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_lines_id;
      hsize_t dataset_dims[] = { 1 };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_write_chunked_dataset(dataset_name,
                                    dataset_data,
                                    dataset_datatype,
                                    dataset_group,
                                    1,
                                    dataset_dims,
                                    dataset_max_dims,
                                    dataset_chunk_dims,
                                    &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Lines/Offsets
     */
    {
      const char* dataset_name = "Offsets";
      int64_t* dataset_data = vtk_pd->lines_offsets;
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_lines_id;
      hsize_t dataset_dims[] = { vtk_pd->n_lines_offsets };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_write_chunked_dataset(dataset_name,
                                    dataset_data,
                                    dataset_datatype,
                                    dataset_group,
                                    1,
                                    dataset_dims,
                                    dataset_max_dims,
                                    dataset_chunk_dims,
                                    &vtk_hdf_pd.vtk_hdf);
    }
  }

  /*
   * Group: /VTKHDF/Polygons
   */
  {
    /*
     * Dataset: /VTKHDF/Polygons/Connectivity
     */
    {
      const char* dataset_name = "Connectivity";
      int64_t* dataset_data = vtk_pd->polygons_connectivity;
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_polygons_id;
      hsize_t dataset_dims[] = { vtk_pd->n_polygons_connectivity };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_write_chunked_dataset(dataset_name,
                                    dataset_data,
                                    dataset_datatype,
                                    dataset_group,
                                    1,
                                    dataset_dims,
                                    dataset_max_dims,
                                    dataset_chunk_dims,
                                    &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Polygons/NumberOfCells
     */
    {
      const char* dataset_name = "NumberOfCells";
      int64_t dataset_data[] = { 0 };
      dataset_data[0] = vtk_polydata_number_of_polygons(vtk_pd);
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_polygons_id;
      hsize_t dataset_dims[] = { 1 };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_write_chunked_dataset(dataset_name,
                                    dataset_data,
                                    dataset_datatype,
                                    dataset_group,
                                    1,
                                    dataset_dims,
                                    dataset_max_dims,
                                    dataset_chunk_dims,
                                    &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Polygons/NumberOfConnectivityIds
     */
    {
      const char* dataset_name = "NumberOfConnectivityIds";
      int64_t dataset_data[] = { (int64_t)vtk_pd->n_polygons_connectivity };
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_polygons_id;
      hsize_t dataset_dims[] = { 1 };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_write_chunked_dataset(dataset_name,
                                    dataset_data,
                                    dataset_datatype,
                                    dataset_group,
                                    1,
                                    dataset_dims,
                                    dataset_max_dims,
                                    dataset_chunk_dims,
                                    &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Polygons/Offsets
     */
    {
      const char* dataset_name = "Offsets";
      int64_t* dataset_data = vtk_pd->polygons_offsets;
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_polygons_id;
      hsize_t dataset_dims[] = { vtk_pd->n_polygons_offsets };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_write_chunked_dataset(dataset_name,
                                    dataset_data,
                                    dataset_datatype,
                                    dataset_group,
                                    1,
                                    dataset_dims,
                                    dataset_max_dims,
                                    dataset_chunk_dims,
                                    &vtk_hdf_pd.vtk_hdf);
    }
  }

  /*
   * Group: /VTKHDF/Strips
   */
  {
    /*
     * Dataset: /VTKHDF/Strips/Connectivity
     */
    {
      const char* dataset_name = "Connectivity";
      int64_t* dataset_data = vtk_pd->strips_connectivity;
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_strips_id;
      hsize_t dataset_dims[] = { vtk_pd->n_strips_connectivity };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_write_chunked_dataset(dataset_name,
                                    dataset_data,
                                    dataset_datatype,
                                    dataset_group,
                                    1,
                                    dataset_dims,
                                    dataset_max_dims,
                                    dataset_chunk_dims,
                                    &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Strips/NumberOfCells
     */
    {
      const char* dataset_name = "NumberOfCells";
      int64_t dataset_data[] = { 0 };
      dataset_data[0] = vtk_polydata_number_of_strips(vtk_pd);
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_strips_id;
      hsize_t dataset_dims[] = { 1 };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_write_chunked_dataset(dataset_name,
                                    dataset_data,
                                    dataset_datatype,
                                    dataset_group,
                                    1,
                                    dataset_dims,
                                    dataset_max_dims,
                                    dataset_chunk_dims,
                                    &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Strips/NumberOfConnectivityIds
     */
    {
      const char* dataset_name = "NumberOfConnectivityIds";
      int64_t dataset_data[] = { (int64_t)vtk_pd->n_strips_connectivity };
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_strips_id;
      hsize_t dataset_dims[] = { 1 };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_write_chunked_dataset(dataset_name,
                                    dataset_data,
                                    dataset_datatype,
                                    dataset_group,
                                    1,
                                    dataset_dims,
                                    dataset_max_dims,
                                    dataset_chunk_dims,
                                    &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Strips/Offsets
     */
    {
      const char* dataset_name = "Offsets";
      int64_t* dataset_data = vtk_pd->strips_offsets;
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_strips_id;
      hsize_t dataset_dims[] = { 1 };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_write_chunked_dataset(dataset_name,
                                    dataset_data,
                                    dataset_datatype,
                                    dataset_group,
                                    1,
                                    dataset_dims,
                                    dataset_max_dims,
                                    dataset_chunk_dims,
                                    &vtk_hdf_pd.vtk_hdf);
    }
  }

  /*
   * Group: /VTKHDF/Vertices
   */
  {
    /*
     * Dataset: /VTKHDF/Vertices/Connectivity
     */
    {
      const char* dataset_name = "Connectivity";
      int64_t* dataset_data = vtk_pd->vertices_connectivity;
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_vertices_id;
      hsize_t dataset_dims[] = { vtk_pd->n_vertices_connectivity };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_write_chunked_dataset(dataset_name,
                                    dataset_data,
                                    dataset_datatype,
                                    dataset_group,
                                    1,
                                    dataset_dims,
                                    dataset_max_dims,
                                    dataset_chunk_dims,
                                    &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Vertices/NumberOfCells
     */
    {
      const char* dataset_name = "NumberOfCells";
      int64_t dataset_data[] = { 0 };
      dataset_data[0] = vtk_polydata_number_of_vertices(vtk_pd);
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_vertices_id;
      hsize_t dataset_dims[] = { 1 };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_write_chunked_dataset(dataset_name,
                                    dataset_data,
                                    dataset_datatype,
                                    dataset_group,
                                    1,
                                    dataset_dims,
                                    dataset_max_dims,
                                    dataset_chunk_dims,
                                    &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Vertices/NumberOfConnectivityIds
     */
    {
      const char* dataset_name = "NumberOfConnectivityIds";
      int64_t dataset_data[] = { (int64_t)vtk_pd->n_vertices_connectivity };
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_vertices_id;
      hsize_t dataset_dims[] = { 1 };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_write_chunked_dataset(dataset_name,
                                    dataset_data,
                                    dataset_datatype,
                                    dataset_group,
                                    1,
                                    dataset_dims,
                                    dataset_max_dims,
                                    dataset_chunk_dims,
                                    &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Vertices/Offsets
     */
    {
      const char* dataset_name = "Offsets";
      int64_t* dataset_data = vtk_pd->vertices_offsets;
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_vertices_id;
      hsize_t dataset_dims[] = { vtk_pd->n_vertices_offsets };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_write_chunked_dataset(dataset_name,
                                    dataset_data,
                                    dataset_datatype,
                                    dataset_group,
                                    1,
                                    dataset_dims,
                                    dataset_max_dims,
                                    dataset_chunk_dims,
                                    &vtk_hdf_pd.vtk_hdf);
    }
  }

  /*
   * Group: /VTKHDF/Steps
   */
  {

    /*
     * Attribute: /VTKHDF/Steps/NSteps
     */
    {
      int64_t attribute_data = 1;
      const char* attribute_name = "NSteps";
      hid_t attribute_datatype = H5T_STD_I64LE;
      hid_t attribute_group = vtk_hdf_pd.grp_steps_id;

      vtk_HDF_write_scalar_attribute(attribute_name,
                                     &attribute_data,
                                     attribute_datatype,
                                     attribute_group,
                                     &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Steps/CellOffsets
     */
    {
      const char* dataset_name = "CellOffsets";
      int64_t dataset_data[] = { 0, 0, 0, 0 };
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_steps_id;
      hsize_t dataset_dims[] = { 1, 4 };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED, 4 };
      hsize_t dataset_chunk_dims[] = { 4000, 4 };

      vtk_HDF_write_chunked_dataset(dataset_name,
                                    dataset_data,
                                    dataset_datatype,
                                    dataset_group,
                                    2,
                                    dataset_dims,
                                    dataset_max_dims,
                                    dataset_chunk_dims,
                                    &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Steps/ConnectivityIdOffsets
     */
    {
      const char* dataset_name = "ConnectivityIdOffsets";
      int64_t dataset_data[] = { 0, 0, 0, 0 };
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_steps_id;
      hsize_t dataset_dims[] = { 1, 4 };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED, 4 };
      hsize_t dataset_chunk_dims[] = { 4000, 4 };

      vtk_HDF_write_chunked_dataset(dataset_name,
                                    dataset_data,
                                    dataset_datatype,
                                    dataset_group,
                                    2,
                                    dataset_dims,
                                    dataset_max_dims,
                                    dataset_chunk_dims,
                                    &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Steps/NumberOfParts
     */
    {
      const char* dataset_name = "NumberOfParts";
      int64_t dataset_data[] = { 1 };
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_steps_id;
      hsize_t dataset_dims[] = { 1 };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_write_chunked_dataset(dataset_name,
                                    dataset_data,
                                    dataset_datatype,
                                    dataset_group,
                                    1,
                                    dataset_dims,
                                    dataset_max_dims,
                                    dataset_chunk_dims,
                                    &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Steps/PartOffsets
     */
    {
      const char* dataset_name = "PartOffsets";
      int64_t dataset_data[] = { 0 };
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_steps_id;
      hsize_t dataset_dims[] = { 1 };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_write_chunked_dataset(dataset_name,
                                    dataset_data,
                                    dataset_datatype,
                                    dataset_group,
                                    1,
                                    dataset_dims,
                                    dataset_max_dims,
                                    dataset_chunk_dims,
                                    &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Steps/PointOffsets
     */
    {
      const char* dataset_name = "PointOffsets";
      int64_t dataset_data[] = { 0 };
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_steps_id;
      hsize_t dataset_dims[] = { 1 };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_write_chunked_dataset(dataset_name,
                                    dataset_data,
                                    dataset_datatype,
                                    dataset_group,
                                    1,
                                    dataset_dims,
                                    dataset_max_dims,
                                    dataset_chunk_dims,
                                    &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Steps/Values
     */
    {
      const char* dataset_name = "Values";
      float dataset_data[] = { time };
      hid_t dataset_datatype = H5T_IEEE_F32LE;
      hid_t dataset_group = vtk_hdf_pd.grp_steps_id;
      hsize_t dataset_dims[] = { 1 };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_write_chunked_dataset(dataset_name,
                                    dataset_data,
                                    dataset_datatype,
                                    dataset_group,
                                    1,
                                    dataset_dims,
                                    dataset_max_dims,
                                    dataset_chunk_dims,
                                    &vtk_hdf_pd.vtk_hdf);
    }
  }

  return vtk_hdf_pd;
}

vtkHDFPolyData
vtk_HDF_polydata_append_transient(const char* fname,
                                  vtkPolyData* vtk_pd,
                                  float time)
{
  vtkHDFPolyData vtk_hdf_pd = {
    .grp_vertices_id = H5I_INVALID_HID,         /**< Vertices group */
    .grp_lines_id = H5I_INVALID_HID,            /**< Lines group */
    .grp_strips_id = H5I_INVALID_HID,           /**< Strips group */
    .grp_polygons_id = H5I_INVALID_HID,         /**< Polygons group */
    .grp_celldata_id = H5I_INVALID_HID,         /**< CellData group */
    .grp_pointdata_id = H5I_INVALID_HID,        /**< PointData group */
    .grp_steps_id = H5I_INVALID_HID,            /**< Steps group */
    .grp_celldata_offsets_id = H5I_INVALID_HID, /**< CellData Offsets group */
    .grp_pointdata_offsets_id = H5I_INVALID_HID /**< PointData Offsets group */
  };

  vtk_hdf_pd.vtk_hdf = vtk_HDF_open(fname);

  vtk_hdf_pd.grp_celldata_id =
    H5Gopen2(vtk_hdf_pd.vtk_hdf.grp_vtkhdf_id, "/VTKHDF/CellData", H5P_DEFAULT);

  vtk_HDF_check_object(&vtk_hdf_pd.vtk_hdf, vtk_hdf_pd.grp_celldata_id);

  vtk_hdf_pd.grp_pointdata_id = H5Gopen2(
    vtk_hdf_pd.vtk_hdf.grp_vtkhdf_id, "/VTKHDF/PointData", H5P_DEFAULT);

  vtk_HDF_check_object(&vtk_hdf_pd.vtk_hdf, vtk_hdf_pd.grp_pointdata_id);

  vtk_hdf_pd.grp_vertices_id =
    H5Gopen2(vtk_hdf_pd.vtk_hdf.grp_vtkhdf_id, "/VTKHDF/Vertices", H5P_DEFAULT);
  vtk_HDF_check_object(&vtk_hdf_pd.vtk_hdf, vtk_hdf_pd.grp_vertices_id);

  vtk_hdf_pd.grp_lines_id =
    H5Gopen2(vtk_hdf_pd.vtk_hdf.grp_vtkhdf_id, "/VTKHDF/Lines", H5P_DEFAULT);

  vtk_HDF_check_object(&vtk_hdf_pd.vtk_hdf, vtk_hdf_pd.grp_lines_id);

  vtk_hdf_pd.grp_strips_id =
    H5Gopen2(vtk_hdf_pd.vtk_hdf.grp_vtkhdf_id, "/VTKHDF/Strips", H5P_DEFAULT);

  vtk_HDF_check_object(&vtk_hdf_pd.vtk_hdf, vtk_hdf_pd.grp_lines_id);

  vtk_hdf_pd.grp_polygons_id =
    H5Gopen2(vtk_hdf_pd.vtk_hdf.grp_vtkhdf_id, "/VTKHDF/Polygons", H5P_DEFAULT);

  vtk_HDF_check_object(&vtk_hdf_pd.vtk_hdf, vtk_hdf_pd.grp_polygons_id);

  vtk_hdf_pd.grp_steps_id =
    H5Gopen2(vtk_hdf_pd.vtk_hdf.grp_vtkhdf_id, "/VTKHDF/Steps", H5P_DEFAULT);

  vtk_HDF_check_object(&vtk_hdf_pd.vtk_hdf, vtk_hdf_pd.grp_steps_id);

  vtk_hdf_pd.grp_celldata_offsets_id = H5Gopen2(
    vtk_hdf_pd.grp_steps_id, "/VTKHDF/Steps/CellDataOffsets", H5P_DEFAULT);

  vtk_HDF_check_object(&vtk_hdf_pd.vtk_hdf, vtk_hdf_pd.grp_celldata_offsets_id);

  vtk_hdf_pd.grp_pointdata_offsets_id = H5Gopen2(
    vtk_hdf_pd.grp_steps_id, "/VTKHDF/Steps/PointDataOffsets", H5P_DEFAULT);

  vtk_HDF_check_object(&vtk_hdf_pd.vtk_hdf,
                       vtk_hdf_pd.grp_pointdata_offsets_id);

  /*
   * Group: /VTKHDF/Steps
   */
  {

    /*
     * Attribute: /VTKHDF/Steps/NSteps
     */
    {
      int64_t attribute_data = 0;
      const char* attribute_name = "NSteps";
      hid_t attribute_datatype = H5T_STD_I64LE;
      hid_t attribute_group = vtk_hdf_pd.grp_steps_id;

      vtk_HDF_read_scalar_attribute(
        attribute_name, &attribute_data, attribute_group, &vtk_hdf_pd.vtk_hdf);

      attribute_data++;

      vtk_HDF_modify_scalar_attribute(attribute_name,
                                     &attribute_data,
                                     attribute_datatype,
                                     attribute_group,
                                     &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Steps/CellOffsets
     */
    {
      const char* dataset_name = "CellOffsets";
      int64_t dataset_data[] = { 0, 0, 0, 0 };
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_steps_id;
      hsize_t dataset_dims[] = { 1, 4 };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED, 4 };
      hsize_t dataset_chunk_dims[] = { 4000, 4 };

      vtk_HDF_append_chunked_dataset(dataset_name,
                                     dataset_data,
                                     dataset_datatype,
                                     dataset_group,
                                     2,
                                     dataset_dims,
                                     &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Steps/ConnectivityIdOffsets
     */
    {
      const char* dataset_name = "ConnectivityIdOffsets";
      int64_t dataset_data[] = { 0, 0, 0, 0 };
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_steps_id;
      hsize_t dataset_dims[] = { 1, 4 };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000, 4 };

      vtk_HDF_append_chunked_dataset(dataset_name,
                                     dataset_data,
                                     dataset_datatype,
                                     dataset_group,
                                     2,
                                     dataset_dims,
                                     &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Steps/NumberOfParts
     *
     * Description: This is the number of parts at each timestep. This appending
     * function is not set up for MPI so we assume there is only 1 part
     */
    {
      const char* dataset_name = "NumberOfParts";
      int64_t dataset_data[] = { 1 };
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_steps_id;
      hsize_t dataset_dims[] = { 1 };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_append_chunked_dataset(dataset_name,
                                     dataset_data,
                                     dataset_datatype,
                                     dataset_group,
                                     1,
                                     dataset_dims,
                                     &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Steps/PartOffsets
     *
     * Description: PartsOffsets is the index in the partition arrays
     * (NumberOfPoints, NumberOfVertices, ... , NumberOfCells) for the first
     * part of the time step t
     */
    {
      // Count the number of previous points
      const char* dataset_name = "PartOffsets";
      int64_t* dataset_data_previous = NULL;
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_steps_id;
      hsize_t dataset_dims_prev[1] = { 0 };
      hsize_t dataset_max_dims_prev[1] = { 0 };

      vtk_HDF_read_dataset(dataset_name,
                           (void**)&dataset_data_previous,
                           dataset_datatype,
                           dataset_group,
                           1,
                           dataset_dims_prev,
                           dataset_max_dims_prev,
                           &vtk_hdf_pd.vtk_hdf);

      int64_t dataset_data[] = { dataset_data_previous[dataset_dims_prev[0]-1] + 1 };
      hsize_t dataset_dims[] = { 1 };

      free(dataset_data_previous);

      vtk_HDF_append_chunked_dataset(dataset_name,
                                     dataset_data,
                                     dataset_datatype,
                                     dataset_group,
                                     1,
                                     dataset_dims,
                                     &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Steps/PointOffsets
     */ 
    {
      // Count the number of previous points
      const char* read_dataset_name = "NumberOfPoints";
      int64_t* read_dataset_data = NULL;
      hid_t read_dataset_datatype = H5T_STD_I64LE;
      hid_t read_dataset_group = vtk_hdf_pd.vtk_hdf.grp_vtkhdf_id;
      hsize_t read_dataset_dims[1] = { 0 };
      hsize_t read_dataset_max_dims[1] = { 0 };

      vtk_HDF_read_dataset(read_dataset_name,
                           (void**)&read_dataset_data,
                           read_dataset_datatype,
                           read_dataset_group,
                           1,
                           read_dataset_dims,
                           read_dataset_max_dims,
                           &vtk_hdf_pd.vtk_hdf);


      const char* write_dataset_name = "PointOffsets";
      int64_t write_dataset_data[] = { 0 };
      
      // Sum all previous elements to get the offset
      for (hsize_t i = 0; i < read_dataset_dims[0]; i++) {
        write_dataset_data[0] += read_dataset_data[i];
      }
      
      free(read_dataset_data);

      hid_t write_dataset_datatype = H5T_STD_I64LE;
      hid_t write_dataset_group = vtk_hdf_pd.grp_steps_id;
      hsize_t write_dataset_dims[] = { 1 };


      vtk_HDF_append_chunked_dataset(write_dataset_name,
                                     write_dataset_data,
                                     write_dataset_datatype,
                                     write_dataset_group,
                                     1,
                                     write_dataset_dims,
                                     &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Steps/Values
     */
    {
      const char* dataset_name = "Values";
      float dataset_data[] = { time };
      hid_t dataset_datatype = H5T_IEEE_F32LE;
      hid_t dataset_group = vtk_hdf_pd.grp_steps_id;
      hsize_t dataset_dims[] = { 1 };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_append_chunked_dataset(dataset_name,
                                     dataset_data,
                                     dataset_datatype,
                                     dataset_group,
                                     1,
                                     dataset_dims,
                                     &vtk_hdf_pd.vtk_hdf);
    }
  }

  /*
   * Dataset: /VTKHDF/Points
   */
  {
    const char* dataset_name = "Points";
    float* dataset_data = vtk_pd->points;
    hid_t dataset_datatype = H5T_IEEE_F32LE;
    hid_t dataset_group = vtk_hdf_pd.vtk_hdf.grp_vtkhdf_id;
    hsize_t dataset_dims[] = { 0, 3 };
    dataset_dims[0] = vtk_polydata_number_of_points(vtk_pd);
    hsize_t dataset_max_dims[] = { H5S_UNLIMITED, 3 };
    hsize_t dataset_chunk_dims[] = { 4000, 3 };

    vtk_HDF_append_chunked_dataset(dataset_name,
                                   dataset_data,
                                   dataset_datatype,
                                   dataset_group,
                                   2,
                                   dataset_dims,
                                   &vtk_hdf_pd.vtk_hdf);
  }

  /*
   * Dataset: /VTKHDF/NumberOfPoints
   */
  {
    const char* dataset_name = "NumberOfPoints";
    int64_t dataset_data[] = { 0 };
    dataset_data[0] = vtk_polydata_number_of_points(vtk_pd);
    hid_t dataset_datatype = H5T_STD_I64LE;
    hid_t dataset_group = vtk_hdf_pd.vtk_hdf.grp_vtkhdf_id;
    hsize_t dataset_dims[] = { 1 };
    hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
    hsize_t dataset_chunk_dims[] = { 4000 };

    vtk_HDF_append_chunked_dataset(dataset_name,
                                  dataset_data,
                                  dataset_datatype,
                                  dataset_group,
                                  1,
                                  dataset_dims,
                                  &vtk_hdf_pd.vtk_hdf);
  }

  /*
   * Group: /VTKHDF/Lines
   */
  {
    /*
     * Dataset: /VTKHDF/Lines/Connectivity
     */
    {
      const char* dataset_name = "Connectivity";
      int64_t* dataset_data = vtk_pd->lines_connectivity;
      hsize_t dataset_dims[] = { vtk_pd->n_lines_connectivity };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_lines_id;

      vtk_HDF_append_chunked_dataset(dataset_name,
                                     dataset_data,
                                     dataset_datatype,
                                     dataset_group,
                                     1,
                                     dataset_dims,
                                     &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Lines/NumberOfCells
     */
    {
      const char* dataset_name = "NumberOfCells";
      int64_t dataset_data[] = { 0 };
      dataset_data[0] = vtk_polydata_number_of_lines(vtk_pd);
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_lines_id;
      hsize_t dataset_dims[] = { 1 };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_append_chunked_dataset(dataset_name,
                                     dataset_data,
                                     dataset_datatype,
                                     dataset_group,
                                     1,
                                     dataset_dims,
                                     &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Lines/NumberOfConnectivityIds
     */
    {
      const char* dataset_name = "NumberOfConnectivityIds";
      int64_t dataset_data[] = { (int64_t)vtk_pd->n_lines_connectivity };
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_lines_id;
      hsize_t dataset_dims[] = { 1 };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_append_chunked_dataset(dataset_name,
                                     dataset_data,
                                     dataset_datatype,
                                     dataset_group,
                                     1,
                                     dataset_dims,
                                     &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Lines/Offsets
     */
    {
      const char* dataset_name = "Offsets";
      int64_t* dataset_data = vtk_pd->lines_offsets;
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_lines_id;
      hsize_t dataset_dims[] = { vtk_pd->n_lines_offsets };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_append_chunked_dataset(dataset_name,
                                     dataset_data,
                                     dataset_datatype,
                                     dataset_group,
                                     1,
                                     dataset_dims,
                                     &vtk_hdf_pd.vtk_hdf);
    }
  }

  /*
   * Group: /VTKHDF/Polygons
   */
  {
    /*
     * Dataset: /VTKHDF/Polygons/Connectivity
     */
    {
      const char* dataset_name = "Connectivity";
      int64_t* dataset_data = vtk_pd->polygons_connectivity;
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_polygons_id;
      hsize_t dataset_dims[] = { vtk_pd->n_polygons_connectivity };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_append_chunked_dataset(dataset_name,
                                     dataset_data,
                                     dataset_datatype,
                                     dataset_group,
                                     1,
                                     dataset_dims,
                                     &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Polygons/NumberOfCells
     */
    {
      const char* dataset_name = "NumberOfCells";
      int64_t dataset_data[] = { 0 };
      dataset_data[0] = vtk_polydata_number_of_polygons(vtk_pd);
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_polygons_id;
      hsize_t dataset_dims[] = { 1 };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_append_chunked_dataset(dataset_name,
                                     dataset_data,
                                     dataset_datatype,
                                     dataset_group,
                                     1,
                                     dataset_dims,
                                     &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Polygons/NumberOfConnectivityIds
     */
    {
      const char* dataset_name = "NumberOfConnectivityIds";
      int64_t dataset_data[] = { (int64_t)vtk_pd->n_polygons_connectivity };
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_polygons_id;
      hsize_t dataset_dims[] = { 1 };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_append_chunked_dataset(dataset_name,
                                     dataset_data,
                                     dataset_datatype,
                                     dataset_group,
                                     1,
                                     dataset_dims,
                                     &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Polygons/Offsets
     */
    {
      const char* dataset_name = "Offsets";
      int64_t* dataset_data = vtk_pd->polygons_offsets;
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_polygons_id;
      hsize_t dataset_dims[] = { vtk_pd->n_polygons_offsets };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_append_chunked_dataset(dataset_name,
                                     dataset_data,
                                     dataset_datatype,
                                     dataset_group,
                                     1,
                                     dataset_dims,
                                     &vtk_hdf_pd.vtk_hdf);
    }
  }

  /*
   * Group: /VTKHDF/Strips
   */
  {
    /*
     * Dataset: /VTKHDF/Strips/Connectivity
     */
    {
      const char* dataset_name = "Connectivity";
      int64_t* dataset_data = vtk_pd->strips_connectivity;
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_strips_id;
      hsize_t dataset_dims[] = { vtk_pd->n_strips_connectivity };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_append_chunked_dataset(dataset_name,
                                     dataset_data,
                                     dataset_datatype,
                                     dataset_group,
                                     1,
                                     dataset_dims,
                                     &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Strips/NumberOfCells
     */
    {
      const char* dataset_name = "NumberOfCells";
      int64_t dataset_data[] = { 0 };
      dataset_data[0] = vtk_polydata_number_of_strips(vtk_pd);
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_strips_id;
      hsize_t dataset_dims[] = { 1 };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_append_chunked_dataset(dataset_name,
                                     dataset_data,
                                     dataset_datatype,
                                     dataset_group,
                                     1,
                                     dataset_dims,
                                     &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Strips/NumberOfConnectivityIds
     */
    {
      const char* dataset_name = "NumberOfConnectivityIds";
      int64_t dataset_data[] = { (int64_t)vtk_pd->n_strips_connectivity };
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_strips_id;
      hsize_t dataset_dims[] = { 1 };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_append_chunked_dataset(dataset_name,
                                     dataset_data,
                                     dataset_datatype,
                                     dataset_group,
                                     1,
                                     dataset_dims,
                                     &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Strips/Offsets
     */
    {
      const char* dataset_name = "Offsets";
      int64_t* dataset_data = vtk_pd->strips_offsets;
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_strips_id;
      hsize_t dataset_dims[] = { 1 };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_append_chunked_dataset(dataset_name,
                                     dataset_data,
                                     dataset_datatype,
                                     dataset_group,
                                     1,
                                     dataset_dims,
                                     &vtk_hdf_pd.vtk_hdf);
    }
  }

  /*
   * Group: /VTKHDF/Vertices
   */
  {
    /*
     * Dataset: /VTKHDF/Vertices/Connectivity
     */
    {
      const char* dataset_name = "Connectivity";
      int64_t* dataset_data = vtk_pd->vertices_connectivity;
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_vertices_id;
      hsize_t dataset_dims[] = { vtk_pd->n_vertices_connectivity };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_append_chunked_dataset(dataset_name,
                                     dataset_data,
                                     dataset_datatype,
                                     dataset_group,
                                     1,
                                     dataset_dims,
                                     &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Vertices/NumberOfCells
     */
    {
      const char* dataset_name = "NumberOfCells";
      int64_t dataset_data[] = { 0 };
      dataset_data[0] = vtk_polydata_number_of_vertices(vtk_pd);
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_vertices_id;
      hsize_t dataset_dims[] = { 1 };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_append_chunked_dataset(dataset_name,
                                     dataset_data,
                                     dataset_datatype,
                                     dataset_group,
                                     1,
                                     dataset_dims,
                                     &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Vertices/NumberOfConnectivityIds
     */
    {
      const char* dataset_name = "NumberOfConnectivityIds";
      int64_t dataset_data[] = { (int64_t)vtk_pd->n_vertices_connectivity };
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_vertices_id;
      hsize_t dataset_dims[] = { 1 };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_append_chunked_dataset(dataset_name,
                                     dataset_data,
                                     dataset_datatype,
                                     dataset_group,
                                     1,
                                     dataset_dims,
                                     &vtk_hdf_pd.vtk_hdf);
    }

    /*
     * Dataset: /VTKHDF/Vertices/Offsets
     */
    {
      const char* dataset_name = "Offsets";
      int64_t* dataset_data = vtk_pd->vertices_offsets;
      hid_t dataset_datatype = H5T_STD_I64LE;
      hid_t dataset_group = vtk_hdf_pd.grp_vertices_id;
      hsize_t dataset_dims[] = { vtk_pd->n_vertices_offsets };
      hsize_t dataset_max_dims[] = { H5S_UNLIMITED };
      hsize_t dataset_chunk_dims[] = { 4000 };

      vtk_HDF_append_chunked_dataset(dataset_name,
                                     dataset_data,
                                     dataset_datatype,
                                     dataset_group,
                                     1,
                                     dataset_dims,
                                     &vtk_hdf_pd.vtk_hdf);
    }
  }

  return vtk_hdf_pd;
}

} // namespace C
} // namespace io
} // namespace beam