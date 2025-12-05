#include "../../Core/vtkPolyDataSet.h"
#include "vtkHDF.h"

#define COMPRESSION 1
#define COMPRESSION_LEVEL 7
#define CHUNK_SIZE (1 << (4))

typedef struct
{
  /* Parent */
  vtkHDF vtk_hdf;

  hid_t grp_celldata_id;
  hid_t grp_lines_id;
  hid_t grp_vertices_id;
  hid_t grp_polygons_id;
  hid_t grp_strips_id;

  /* Dataspace, datatype, and property-list identifiers */
  hid_t attr_space_id;
  hid_t attr_dtype_id;
  hid_t dset_space_id;
  hid_t dcpl_id;
  hid_t dset_dtype_id;

  /* Dataset and attribute identifiers */
  hid_t dset_id;
  hid_t attr_id;

  /* MPIIO */
  hid_t file_space;
  hid_t mem_space;
  hid_t xfer_plist;

} vtkHDFPolyData;

void
vtk_HDF_polydata_close(vtkHDFPolyData* vtk_hdf_polydata)
{
  //
  if (vtk_hdf_polydata->grp_celldata_id >= 0)
    H5Gclose(vtk_hdf_polydata->grp_celldata_id);

  if (vtk_hdf_polydata->xfer_plist >= 0)
    H5Pclose(vtk_hdf_polydata->xfer_plist);

  // #if _MPI
  if (vtk_hdf_polydata->mem_space >= 0)
    H5Sclose(vtk_hdf_polydata->mem_space);

  if (vtk_hdf_polydata->file_space >= 0)
    H5Sclose(vtk_hdf_polydata->file_space);
  // #endif

  vtk_HDF_close(&vtk_hdf_polydata->vtk_hdf);
}

void
vtk_HDF_polydata_error(vtkHDFPolyData* vtk_hdf_polydata)
{
  //
  vtk_HDF_polydata_close(vtk_hdf_polydata);
  assert(1 == 2);
}

vtkHDFPolyData
vtk_HDF_polydata_init(vtkPolyData_t* data, const char* fname)
{

  vtkHDF vtk_hdf = vtk_HDF_init(fname);
  vtkHDFPolyData vtk_hdf_polydata = {
    .vtk_hdf = vtk_hdf,

    .grp_celldata_id = -1,
    .grp_vertices_id = -1,
    .grp_lines_id = -1,
    .grp_polygons_id = -1,
    .grp_strips_id = -1,

    .attr_space_id = -1,
    .attr_dtype_id = -1,
    .attr_id = -1,

    .dset_space_id = -1,
    .dcpl_id = -1,
    .dset_dtype_id = -1,
    .dset_id = -1,

    .mem_space = -1,
    .file_space = -1,
  };

  /**
  # Type:
  This needs to be a fixed-length string of size 8, value =
  "PolyData"
  */
  // if (pid() == 0)
  {
    const char* type_str = "PolyData";
    // hsize_t scalar_dims = 1;
    vtk_hdf_polydata.attr_space_id = H5Screate(H5S_SCALAR);
    if (vtk_hdf_polydata.attr_space_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    /* Create a fixed-length string datatype of length 8, null-padded, ASCII */
    vtk_hdf_polydata.attr_dtype_id = H5Tcopy(H5T_C_S1);
    if (vtk_hdf_polydata.attr_dtype_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    if (H5Tset_size(vtk_hdf_polydata.attr_dtype_id, (size_t)8) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);
    if (H5Tset_strpad(vtk_hdf_polydata.attr_dtype_id, H5T_STR_NULLPAD) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);
    if (H5Tset_cset(vtk_hdf_polydata.attr_dtype_id, H5T_CSET_ASCII) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.attr_id =
      H5Acreate2(vtk_hdf_polydata.vtk_hdf.grp_vtkhdf_id,
                 "Type",
                 vtk_hdf_polydata.attr_dtype_id,
                 vtk_hdf_polydata.attr_space_id,
                 H5P_DEFAULT,
                 H5P_DEFAULT);
    if (vtk_hdf_polydata.attr_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    /* Write the string (automatically null‐padded up to length 13) */
    if (H5Awrite(vtk_hdf_polydata.attr_id,
                 vtk_hdf_polydata.attr_dtype_id,
                 type_str) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    H5Aclose(vtk_hdf_polydata.attr_id);
    H5Tclose(vtk_hdf_polydata.attr_dtype_id);
    H5Sclose(vtk_hdf_polydata.attr_space_id);
  }

  /**
  # Version
  2-element int64 array {2,4}
  */
  // if (pid() == 0)
  {
    int64_t vers_value[2] = { 2, 0 };
    hsize_t dims_attr[1] = { 2 };
    vtk_hdf_polydata.attr_space_id = H5Screate_simple(1, dims_attr, dims_attr);
    if (vtk_hdf_polydata.attr_space_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.attr_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_polydata.attr_dtype_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.attr_id =
      H5Acreate2(vtk_hdf_polydata.vtk_hdf.grp_vtkhdf_id,
                 "Version",
                 vtk_hdf_polydata.attr_dtype_id,
                 vtk_hdf_polydata.attr_space_id,
                 H5P_DEFAULT,
                 H5P_DEFAULT);
    if (vtk_hdf_polydata.attr_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    if (H5Awrite(vtk_hdf_polydata.attr_id,
                 vtk_hdf_polydata.attr_dtype_id,
                 vers_value) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    H5Aclose(vtk_hdf_polydata.attr_id);
    H5Tclose(vtk_hdf_polydata.attr_dtype_id);
    H5Sclose(vtk_hdf_polydata.attr_space_id);
  }

  /*
    Number Of Points
  */

  {
    int64_t number_of_points[] = { data->points.n_points };

    hsize_t dims_d[1] = { 1 };
    hsize_t maxdims_d[1] = { H5S_UNLIMITED };
    hsize_t chunk_dims[1] = { 1 };

    vtk_hdf_polydata.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_polydata.dset_space_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_polydata.dcpl_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);
    if (H5Pset_chunk(vtk_hdf_polydata.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_polydata.dset_dtype_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_id =
      H5Dcreate2(vtk_hdf_polydata.vtk_hdf.grp_vtkhdf_id,
                 "NumberOfPoints",
                 vtk_hdf_polydata.dset_dtype_id,
                 vtk_hdf_polydata.dset_space_id,
                 H5P_DEFAULT,
                 vtk_hdf_polydata.dcpl_id,
                 H5P_DEFAULT);
    if (vtk_hdf_polydata.dset_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    if (H5Dwrite(vtk_hdf_polydata.dset_id,
                 vtk_hdf_polydata.dset_dtype_id,
                 H5S_ALL,
                 H5S_ALL,
                 H5P_DEFAULT,
                 number_of_points) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    H5Dclose(vtk_hdf_polydata.dset_id);
    H5Tclose(vtk_hdf_polydata.dset_dtype_id);
    H5Pclose(vtk_hdf_polydata.dcpl_id);
    H5Sclose(vtk_hdf_polydata.dset_space_id);
  }

  /*
    Points
  */

  {
    float* points = data->points.data;

    hsize_t dims_d[2] = { (hsize_t)data->points.n_points,
                          (hsize_t)data->points.n_components };
    hsize_t maxdims_d[2] = { H5S_UNLIMITED,
                             (hsize_t)data->points.n_components };
    hsize_t chunk_dims[2] = { 1, (hsize_t)data->points.n_components };

    vtk_hdf_polydata.dset_space_id = H5Screate_simple(2, dims_d, maxdims_d);
    if (vtk_hdf_polydata.dset_space_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_polydata.dcpl_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);
    if (H5Pset_chunk(vtk_hdf_polydata.dcpl_id, 2, chunk_dims) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_dtype_id = H5Tcopy(H5T_IEEE_F32LE);
    if (vtk_hdf_polydata.dset_dtype_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_id =
      H5Dcreate2(vtk_hdf_polydata.vtk_hdf.grp_vtkhdf_id,
                 "Points",
                 vtk_hdf_polydata.dset_dtype_id,
                 vtk_hdf_polydata.dset_space_id,
                 H5P_DEFAULT,
                 vtk_hdf_polydata.dcpl_id,
                 H5P_DEFAULT);
    if (vtk_hdf_polydata.dset_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    if (H5Dwrite(vtk_hdf_polydata.dset_id,
                 vtk_hdf_polydata.dset_dtype_id,
                 H5S_ALL,
                 H5S_ALL,
                 H5P_DEFAULT,
                 points) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    H5Dclose(vtk_hdf_polydata.dset_id);
    H5Tclose(vtk_hdf_polydata.dset_dtype_id);
    H5Pclose(vtk_hdf_polydata.dcpl_id);
    H5Sclose(vtk_hdf_polydata.dset_space_id);
  }

  /*
   * Vertices
   */

  vtk_hdf_polydata.grp_vertices_id =
    H5Gcreate2(vtk_hdf_polydata.vtk_hdf.grp_vtkhdf_id,
               "Vertices",
               H5P_DEFAULT,
               H5P_DEFAULT,
               H5P_DEFAULT);
  if (vtk_hdf_polydata.grp_vertices_id < 0)
    vtk_HDF_polydata_error(&vtk_hdf_polydata);

  /*
   * Vertices / Connectivity
   */

  {
    int64_t *connectivity = data->vertices.connectivity;

    hsize_t dims_d[1] = { data->vertices.n_cells };
    hsize_t maxdims_d[1] = { H5S_UNLIMITED };
    hsize_t chunk_dims[1] = { 1 };

    vtk_hdf_polydata.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_polydata.dset_space_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_polydata.dcpl_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);
    if (H5Pset_chunk(vtk_hdf_polydata.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_polydata.dset_dtype_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_id = H5Dcreate2(vtk_hdf_polydata.grp_vertices_id,
                                          "Connectivity",
                                          vtk_hdf_polydata.dset_dtype_id,
                                          vtk_hdf_polydata.dset_space_id,
                                          H5P_DEFAULT,
                                          vtk_hdf_polydata.dcpl_id,
                                          H5P_DEFAULT);
    if (vtk_hdf_polydata.dset_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    if (H5Dwrite(vtk_hdf_polydata.dset_id,
                 vtk_hdf_polydata.dset_dtype_id,
                 H5S_ALL,
                 H5S_ALL,
                 H5P_DEFAULT,
                 connectivity) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    H5Dclose(vtk_hdf_polydata.dset_id);
    H5Tclose(vtk_hdf_polydata.dset_dtype_id);
    H5Pclose(vtk_hdf_polydata.dcpl_id);
    H5Sclose(vtk_hdf_polydata.dset_space_id);
  }

  /*
   * Vertices / NumberOfCells
   */

  {
    int64_t number_of_cells = data->vertices.n_cells;

    hsize_t dims_d[1] = { 1 };
    hsize_t maxdims_d[1] = { H5S_UNLIMITED };
    hsize_t chunk_dims[1] = { 1 };

    vtk_hdf_polydata.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_polydata.dset_space_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_polydata.dcpl_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);
    if (H5Pset_chunk(vtk_hdf_polydata.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_polydata.dset_dtype_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_id = H5Dcreate2(vtk_hdf_polydata.grp_vertices_id,
                                          "NumberOfCells",
                                          vtk_hdf_polydata.dset_dtype_id,
                                          vtk_hdf_polydata.dset_space_id,
                                          H5P_DEFAULT,
                                          vtk_hdf_polydata.dcpl_id,
                                          H5P_DEFAULT);
    if (vtk_hdf_polydata.dset_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    if (H5Dwrite(vtk_hdf_polydata.dset_id,
                 vtk_hdf_polydata.dset_dtype_id,
                 H5S_ALL,
                 H5S_ALL,
                 H5P_DEFAULT,
                 &number_of_cells) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    H5Dclose(vtk_hdf_polydata.dset_id);
    H5Tclose(vtk_hdf_polydata.dset_dtype_id);
    H5Pclose(vtk_hdf_polydata.dcpl_id);
    H5Sclose(vtk_hdf_polydata.dset_space_id);
  }

  /*
   * Vertices / NumberOfConnectivityIds
   */

  {
    int64_t number_of_connectivity_ids = data->vertices.n_cells;

    hsize_t dims_d[1] = { 1 };
    hsize_t maxdims_d[1] = { H5S_UNLIMITED };
    hsize_t chunk_dims[1] = { 1 };

    vtk_hdf_polydata.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_polydata.dset_space_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_polydata.dcpl_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);
    if (H5Pset_chunk(vtk_hdf_polydata.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_polydata.dset_dtype_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_id = H5Dcreate2(vtk_hdf_polydata.grp_vertices_id,
                                          "NumberOfConnectivityIds",
                                          vtk_hdf_polydata.dset_dtype_id,
                                          vtk_hdf_polydata.dset_space_id,
                                          H5P_DEFAULT,
                                          vtk_hdf_polydata.dcpl_id,
                                          H5P_DEFAULT);
    if (vtk_hdf_polydata.dset_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    if (H5Dwrite(vtk_hdf_polydata.dset_id,
                 vtk_hdf_polydata.dset_dtype_id,
                 H5S_ALL,
                 H5S_ALL,
                 H5P_DEFAULT,
                 &number_of_connectivity_ids) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    H5Dclose(vtk_hdf_polydata.dset_id);
    H5Tclose(vtk_hdf_polydata.dset_dtype_id);
    H5Pclose(vtk_hdf_polydata.dcpl_id);
    H5Sclose(vtk_hdf_polydata.dset_space_id);
  }

  /*
   * Vertices / Offsets
   */

  {
    int64_t *offsets = data->vertices.offsets;

    hsize_t dims_d[1] = { data->vertices.n_cells + 1};
    hsize_t maxdims_d[1] = { H5S_UNLIMITED };
    hsize_t chunk_dims[1] = {1};

    vtk_hdf_polydata.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_polydata.dset_space_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_polydata.dcpl_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);
    if (H5Pset_chunk(vtk_hdf_polydata.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_polydata.dset_dtype_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_id = H5Dcreate2(vtk_hdf_polydata.grp_vertices_id,
                                          "Offsets",
                                          vtk_hdf_polydata.dset_dtype_id,
                                          vtk_hdf_polydata.dset_space_id,
                                          H5P_DEFAULT,
                                          vtk_hdf_polydata.dcpl_id,
                                          H5P_DEFAULT);
    if (vtk_hdf_polydata.dset_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    if (H5Dwrite(vtk_hdf_polydata.dset_id,
                 vtk_hdf_polydata.dset_dtype_id,
                 H5S_ALL,
                 H5S_ALL,
                 H5P_DEFAULT,
                 offsets) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    H5Dclose(vtk_hdf_polydata.dset_id);
    H5Tclose(vtk_hdf_polydata.dset_dtype_id);
    H5Pclose(vtk_hdf_polydata.dcpl_id);
    H5Sclose(vtk_hdf_polydata.dset_space_id);
  }

  /*
   * Lines
   */

  vtk_hdf_polydata.grp_lines_id =
    H5Gcreate2(vtk_hdf_polydata.vtk_hdf.grp_vtkhdf_id,
               "Lines",
               H5P_DEFAULT,
               H5P_DEFAULT,
               H5P_DEFAULT);
  if (vtk_hdf_polydata.grp_lines_id < 0)
    vtk_HDF_polydata_error(&vtk_hdf_polydata);

  /*
   * Lines / Connectivity
   */

  {
    //int64_t connectivity[] = {};
    int64_t *connectivity = NULL;

    hsize_t dims_d[1] = { 0 };
    hsize_t maxdims_d[1] = { H5S_UNLIMITED };
    hsize_t chunk_dims[1] = { 1 };

    vtk_hdf_polydata.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_polydata.dset_space_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_polydata.dcpl_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);
    if (H5Pset_chunk(vtk_hdf_polydata.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_polydata.dset_dtype_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_id = H5Dcreate2(vtk_hdf_polydata.grp_lines_id,
                                          "Connectivity",
                                          vtk_hdf_polydata.dset_dtype_id,
                                          vtk_hdf_polydata.dset_space_id,
                                          H5P_DEFAULT,
                                          vtk_hdf_polydata.dcpl_id,
                                          H5P_DEFAULT);
    if (vtk_hdf_polydata.dset_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    if (H5Dwrite(vtk_hdf_polydata.dset_id,
                 vtk_hdf_polydata.dset_dtype_id,
                 H5S_ALL,
                 H5S_ALL,
                 H5P_DEFAULT,
                 connectivity) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    H5Dclose(vtk_hdf_polydata.dset_id);
    H5Tclose(vtk_hdf_polydata.dset_dtype_id);
    H5Pclose(vtk_hdf_polydata.dcpl_id);
    H5Sclose(vtk_hdf_polydata.dset_space_id);
  }

  /*
   * Lines / NumberOfCells
   */

  {
    int64_t number_of_cells = 0;

    hsize_t dims_d[1] = { 1 };
    hsize_t maxdims_d[1] = { H5S_UNLIMITED };
    hsize_t chunk_dims[1] = { 1 };

    vtk_hdf_polydata.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_polydata.dset_space_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_polydata.dcpl_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);
    if (H5Pset_chunk(vtk_hdf_polydata.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_polydata.dset_dtype_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_id = H5Dcreate2(vtk_hdf_polydata.grp_lines_id,
                                          "NumberOfCells",
                                          vtk_hdf_polydata.dset_dtype_id,
                                          vtk_hdf_polydata.dset_space_id,
                                          H5P_DEFAULT,
                                          vtk_hdf_polydata.dcpl_id,
                                          H5P_DEFAULT);
    if (vtk_hdf_polydata.dset_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    if (H5Dwrite(vtk_hdf_polydata.dset_id,
                 vtk_hdf_polydata.dset_dtype_id,
                 H5S_ALL,
                 H5S_ALL,
                 H5P_DEFAULT,
                 &number_of_cells) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    H5Dclose(vtk_hdf_polydata.dset_id);
    H5Tclose(vtk_hdf_polydata.dset_dtype_id);
    H5Pclose(vtk_hdf_polydata.dcpl_id);
    H5Sclose(vtk_hdf_polydata.dset_space_id);
  }

  /*
   * Lines / NumberOfConnectivityIds
   */

  {
    int64_t number_of_connectivity_ids = 0;

    hsize_t dims_d[1] = { 1 };
    hsize_t maxdims_d[1] = { H5S_UNLIMITED };
    hsize_t chunk_dims[1] = { 1 };

    vtk_hdf_polydata.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_polydata.dset_space_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_polydata.dcpl_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);
    if (H5Pset_chunk(vtk_hdf_polydata.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_polydata.dset_dtype_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_id = H5Dcreate2(vtk_hdf_polydata.grp_lines_id,
                                          "NumberOfConnectivityIds",
                                          vtk_hdf_polydata.dset_dtype_id,
                                          vtk_hdf_polydata.dset_space_id,
                                          H5P_DEFAULT,
                                          vtk_hdf_polydata.dcpl_id,
                                          H5P_DEFAULT);
    if (vtk_hdf_polydata.dset_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    if (H5Dwrite(vtk_hdf_polydata.dset_id,
                 vtk_hdf_polydata.dset_dtype_id,
                 H5S_ALL,
                 H5S_ALL,
                 H5P_DEFAULT,
                 &number_of_connectivity_ids) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    H5Dclose(vtk_hdf_polydata.dset_id);
    H5Tclose(vtk_hdf_polydata.dset_dtype_id);
    H5Pclose(vtk_hdf_polydata.dcpl_id);
    H5Sclose(vtk_hdf_polydata.dset_space_id);
  }

  /*
   * Lines / NumberOfConnectivityIds
   */

  {
    int64_t offsets = 0;

    hsize_t dims_d[1] = { 1 };
    hsize_t maxdims_d[1] = { H5S_UNLIMITED };
    hsize_t chunk_dims[1] = { 1 };

    vtk_hdf_polydata.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_polydata.dset_space_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_polydata.dcpl_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);
    if (H5Pset_chunk(vtk_hdf_polydata.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_polydata.dset_dtype_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_id = H5Dcreate2(vtk_hdf_polydata.grp_lines_id,
                                          "Offsets",
                                          vtk_hdf_polydata.dset_dtype_id,
                                          vtk_hdf_polydata.dset_space_id,
                                          H5P_DEFAULT,
                                          vtk_hdf_polydata.dcpl_id,
                                          H5P_DEFAULT);
    if (vtk_hdf_polydata.dset_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    if (H5Dwrite(vtk_hdf_polydata.dset_id,
                 vtk_hdf_polydata.dset_dtype_id,
                 H5S_ALL,
                 H5S_ALL,
                 H5P_DEFAULT,
                 &offsets) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    H5Dclose(vtk_hdf_polydata.dset_id);
    H5Tclose(vtk_hdf_polydata.dset_dtype_id);
    H5Pclose(vtk_hdf_polydata.dcpl_id);
    H5Sclose(vtk_hdf_polydata.dset_space_id);
  }

  
  /*
   * Polygons
   */

  vtk_hdf_polydata.grp_polygons_id =
    H5Gcreate2(vtk_hdf_polydata.vtk_hdf.grp_vtkhdf_id,
               "Polygons",
               H5P_DEFAULT,
               H5P_DEFAULT,
               H5P_DEFAULT);
  if (vtk_hdf_polydata.grp_polygons_id < 0)
    vtk_HDF_polydata_error(&vtk_hdf_polydata);

  /*
   * Polygons / Connectivity
   */

  {
    //int64_t connectivity[] = {};
    int64_t *connectivity = NULL;

    hsize_t dims_d[1] = { 0 };
    hsize_t maxdims_d[1] = { H5S_UNLIMITED };
    hsize_t chunk_dims[1] = { 1 };

    vtk_hdf_polydata.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_polydata.dset_space_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_polydata.dcpl_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);
    if (H5Pset_chunk(vtk_hdf_polydata.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_polydata.dset_dtype_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_id = H5Dcreate2(vtk_hdf_polydata.grp_polygons_id,
                                          "Connectivity",
                                          vtk_hdf_polydata.dset_dtype_id,
                                          vtk_hdf_polydata.dset_space_id,
                                          H5P_DEFAULT,
                                          vtk_hdf_polydata.dcpl_id,
                                          H5P_DEFAULT);
    if (vtk_hdf_polydata.dset_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    if (H5Dwrite(vtk_hdf_polydata.dset_id,
                 vtk_hdf_polydata.dset_dtype_id,
                 H5S_ALL,
                 H5S_ALL,
                 H5P_DEFAULT,
                 connectivity) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    H5Dclose(vtk_hdf_polydata.dset_id);
    H5Tclose(vtk_hdf_polydata.dset_dtype_id);
    H5Pclose(vtk_hdf_polydata.dcpl_id);
    H5Sclose(vtk_hdf_polydata.dset_space_id);
  }

  /*
   * Polygons / NumberOfCells
   */

  {
    int64_t number_of_cells = 0;

    hsize_t dims_d[1] = { 1 };
    hsize_t maxdims_d[1] = { H5S_UNLIMITED };
    hsize_t chunk_dims[1] = { 1 };

    vtk_hdf_polydata.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_polydata.dset_space_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_polydata.dcpl_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);
    if (H5Pset_chunk(vtk_hdf_polydata.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_polydata.dset_dtype_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_id = H5Dcreate2(vtk_hdf_polydata.grp_polygons_id,
                                          "NumberOfCells",
                                          vtk_hdf_polydata.dset_dtype_id,
                                          vtk_hdf_polydata.dset_space_id,
                                          H5P_DEFAULT,
                                          vtk_hdf_polydata.dcpl_id,
                                          H5P_DEFAULT);
    if (vtk_hdf_polydata.dset_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    if (H5Dwrite(vtk_hdf_polydata.dset_id,
                 vtk_hdf_polydata.dset_dtype_id,
                 H5S_ALL,
                 H5S_ALL,
                 H5P_DEFAULT,
                 &number_of_cells) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    H5Dclose(vtk_hdf_polydata.dset_id);
    H5Tclose(vtk_hdf_polydata.dset_dtype_id);
    H5Pclose(vtk_hdf_polydata.dcpl_id);
    H5Sclose(vtk_hdf_polydata.dset_space_id);
  }

  /*
   * Polygons / NumberOfConnectivityIds
   */

  {
    int64_t number_of_connectivity_ids = 0;

    hsize_t dims_d[1] = { 1 };
    hsize_t maxdims_d[1] = { H5S_UNLIMITED };
    hsize_t chunk_dims[1] = { 1 };

    vtk_hdf_polydata.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_polydata.dset_space_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_polydata.dcpl_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);
    if (H5Pset_chunk(vtk_hdf_polydata.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_polydata.dset_dtype_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_id = H5Dcreate2(vtk_hdf_polydata.grp_polygons_id,
                                          "NumberOfConnectivityIds",
                                          vtk_hdf_polydata.dset_dtype_id,
                                          vtk_hdf_polydata.dset_space_id,
                                          H5P_DEFAULT,
                                          vtk_hdf_polydata.dcpl_id,
                                          H5P_DEFAULT);
    if (vtk_hdf_polydata.dset_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    if (H5Dwrite(vtk_hdf_polydata.dset_id,
                 vtk_hdf_polydata.dset_dtype_id,
                 H5S_ALL,
                 H5S_ALL,
                 H5P_DEFAULT,
                 &number_of_connectivity_ids) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    H5Dclose(vtk_hdf_polydata.dset_id);
    H5Tclose(vtk_hdf_polydata.dset_dtype_id);
    H5Pclose(vtk_hdf_polydata.dcpl_id);
    H5Sclose(vtk_hdf_polydata.dset_space_id);
  }

  /*
   * Polygons / NumberOfConnectivityIds
   */

  {
    int64_t offsets = 0;

    hsize_t dims_d[1] = { 1 };
    hsize_t maxdims_d[1] = { H5S_UNLIMITED };
    hsize_t chunk_dims[1] = { 1 };

    vtk_hdf_polydata.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_polydata.dset_space_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_polydata.dcpl_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);
    if (H5Pset_chunk(vtk_hdf_polydata.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_polydata.dset_dtype_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_id = H5Dcreate2(vtk_hdf_polydata.grp_polygons_id,
                                          "Offsets",
                                          vtk_hdf_polydata.dset_dtype_id,
                                          vtk_hdf_polydata.dset_space_id,
                                          H5P_DEFAULT,
                                          vtk_hdf_polydata.dcpl_id,
                                          H5P_DEFAULT);
    if (vtk_hdf_polydata.dset_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    if (H5Dwrite(vtk_hdf_polydata.dset_id,
                 vtk_hdf_polydata.dset_dtype_id,
                 H5S_ALL,
                 H5S_ALL,
                 H5P_DEFAULT,
                 &offsets) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    H5Dclose(vtk_hdf_polydata.dset_id);
    H5Tclose(vtk_hdf_polydata.dset_dtype_id);
    H5Pclose(vtk_hdf_polydata.dcpl_id);
    H5Sclose(vtk_hdf_polydata.dset_space_id);
  }

  
  /*
   * Strips
   */

  vtk_hdf_polydata.grp_strips_id =
    H5Gcreate2(vtk_hdf_polydata.vtk_hdf.grp_vtkhdf_id,
               "Strips",
               H5P_DEFAULT,
               H5P_DEFAULT,
               H5P_DEFAULT);
  if (vtk_hdf_polydata.grp_strips_id < 0)
    vtk_HDF_polydata_error(&vtk_hdf_polydata);

  /*
   * Strips / Connectivity
   */

  {
    //int64_t connectivity[] = {};
    int64_t *connectivity = NULL;

    hsize_t dims_d[1] = { 0 };
    hsize_t maxdims_d[1] = { H5S_UNLIMITED };
    hsize_t chunk_dims[1] = { 1 };

    vtk_hdf_polydata.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_polydata.dset_space_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_polydata.dcpl_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);
    if (H5Pset_chunk(vtk_hdf_polydata.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_polydata.dset_dtype_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_id = H5Dcreate2(vtk_hdf_polydata.grp_strips_id,
                                          "Connectivity",
                                          vtk_hdf_polydata.dset_dtype_id,
                                          vtk_hdf_polydata.dset_space_id,
                                          H5P_DEFAULT,
                                          vtk_hdf_polydata.dcpl_id,
                                          H5P_DEFAULT);
    if (vtk_hdf_polydata.dset_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    if (H5Dwrite(vtk_hdf_polydata.dset_id,
                 vtk_hdf_polydata.dset_dtype_id,
                 H5S_ALL,
                 H5S_ALL,
                 H5P_DEFAULT,
                 connectivity) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    H5Dclose(vtk_hdf_polydata.dset_id);
    H5Tclose(vtk_hdf_polydata.dset_dtype_id);
    H5Pclose(vtk_hdf_polydata.dcpl_id);
    H5Sclose(vtk_hdf_polydata.dset_space_id);
  }

  /*
   * Strips / NumberOfCells
   */

  {
    int64_t number_of_cells = 0;

    hsize_t dims_d[1] = { 1 };
    hsize_t maxdims_d[1] = { H5S_UNLIMITED };
    hsize_t chunk_dims[1] = { 1 };

    vtk_hdf_polydata.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_polydata.dset_space_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_polydata.dcpl_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);
    if (H5Pset_chunk(vtk_hdf_polydata.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_polydata.dset_dtype_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_id = H5Dcreate2(vtk_hdf_polydata.grp_strips_id,
                                          "NumberOfCells",
                                          vtk_hdf_polydata.dset_dtype_id,
                                          vtk_hdf_polydata.dset_space_id,
                                          H5P_DEFAULT,
                                          vtk_hdf_polydata.dcpl_id,
                                          H5P_DEFAULT);
    if (vtk_hdf_polydata.dset_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    if (H5Dwrite(vtk_hdf_polydata.dset_id,
                 vtk_hdf_polydata.dset_dtype_id,
                 H5S_ALL,
                 H5S_ALL,
                 H5P_DEFAULT,
                 &number_of_cells) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    H5Dclose(vtk_hdf_polydata.dset_id);
    H5Tclose(vtk_hdf_polydata.dset_dtype_id);
    H5Pclose(vtk_hdf_polydata.dcpl_id);
    H5Sclose(vtk_hdf_polydata.dset_space_id);
  }

  /*
   * Strips / NumberOfConnectivityIds
   */

  {
    int64_t number_of_connectivity_ids = 0;

    hsize_t dims_d[1] = { 1 };
    hsize_t maxdims_d[1] = { H5S_UNLIMITED };
    hsize_t chunk_dims[1] = { 1 };

    vtk_hdf_polydata.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_polydata.dset_space_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_polydata.dcpl_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);
    if (H5Pset_chunk(vtk_hdf_polydata.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_polydata.dset_dtype_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_id = H5Dcreate2(vtk_hdf_polydata.grp_strips_id,
                                          "NumberOfConnectivityIds",
                                          vtk_hdf_polydata.dset_dtype_id,
                                          vtk_hdf_polydata.dset_space_id,
                                          H5P_DEFAULT,
                                          vtk_hdf_polydata.dcpl_id,
                                          H5P_DEFAULT);
    if (vtk_hdf_polydata.dset_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    if (H5Dwrite(vtk_hdf_polydata.dset_id,
                 vtk_hdf_polydata.dset_dtype_id,
                 H5S_ALL,
                 H5S_ALL,
                 H5P_DEFAULT,
                 &number_of_connectivity_ids) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    H5Dclose(vtk_hdf_polydata.dset_id);
    H5Tclose(vtk_hdf_polydata.dset_dtype_id);
    H5Pclose(vtk_hdf_polydata.dcpl_id);
    H5Sclose(vtk_hdf_polydata.dset_space_id);
  }

  /*
   * Strips / NumberOfConnectivityIds
   */

  {
    int64_t offsets = 0;

    hsize_t dims_d[1] = { 1 };
    hsize_t maxdims_d[1] = { H5S_UNLIMITED };
    hsize_t chunk_dims[1] = { 1 };

    vtk_hdf_polydata.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_polydata.dset_space_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_polydata.dcpl_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);
    if (H5Pset_chunk(vtk_hdf_polydata.dcpl_id, 1, chunk_dims) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_polydata.dset_dtype_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    vtk_hdf_polydata.dset_id = H5Dcreate2(vtk_hdf_polydata.grp_strips_id,
                                          "Offsets",
                                          vtk_hdf_polydata.dset_dtype_id,
                                          vtk_hdf_polydata.dset_space_id,
                                          H5P_DEFAULT,
                                          vtk_hdf_polydata.dcpl_id,
                                          H5P_DEFAULT);
    if (vtk_hdf_polydata.dset_id < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    if (H5Dwrite(vtk_hdf_polydata.dset_id,
                 vtk_hdf_polydata.dset_dtype_id,
                 H5S_ALL,
                 H5S_ALL,
                 H5P_DEFAULT,
                 &offsets) < 0)
      vtk_HDF_polydata_error(&vtk_hdf_polydata);

    H5Dclose(vtk_hdf_polydata.dset_id);
    H5Tclose(vtk_hdf_polydata.dset_dtype_id);
    H5Pclose(vtk_hdf_polydata.dcpl_id);
    H5Sclose(vtk_hdf_polydata.dset_space_id);
  }

  

  return vtk_hdf_polydata;
}