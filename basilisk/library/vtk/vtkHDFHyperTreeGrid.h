#include "vtkHDF.h"
#include "vtkHDFHyperTreeGridData.h"

#include <float.h>
#include <math.h>

#define COMPRESSION 1
#define COMPRESSION_LEVEL 7

/**
 * @brief This struct holds various IDs needed by the HDF5 library to read and
 * write HDF5 files for our particular HyperTreeGrid/PHyperTreeGrid schema.
 */
typedef struct
{
  /* Parent */
  vtkHDF vtk_hdf;

  hid_t grp_celldata_id;
  // hid_t grp_steps_id;

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

} vtkHDFHyperTreeGrid;

/**
 * @brief This function closes all objects (if any are still open).
 *
 * @param vtk_hdf_htg Struct holding the object ids for our file
 *
 * @memberof vtkHDFHyperTreeGrid
 */
void
vtk_HDF_hypertreegrid_close(vtkHDFHyperTreeGrid* vtk_hdf_htg)
{
  if (vtk_hdf_htg->grp_celldata_id >= 0)
    H5Gclose(vtk_hdf_htg->grp_celldata_id);

  if (vtk_hdf_htg->xfer_plist >= 0)
    H5Pclose(vtk_hdf_htg->xfer_plist);

  if (vtk_hdf_htg->mem_space >= 0)
    H5Sclose(vtk_hdf_htg->mem_space);

  if (vtk_hdf_htg->file_space >= 0)
    H5Sclose(vtk_hdf_htg->file_space);

  vtk_HDF_close(&vtk_hdf_htg->vtk_hdf);
}

/**
 * @brief This function is called on error, and immediately closes any open
 * objects, so that information already written is saved correctly. Then an
 * error is raised to alert the user.
 *
 * @param vtk_hdf_htg Struct holding the object ids for our file
 *
 * @memberof vtkHDFHyperTreeGrid
 */
void
vtk_HDF_hypertreegrid_error(vtkHDFHyperTreeGrid* vtk_hdf_htg)
{
  vtk_HDF_hypertreegrid_close(vtk_hdf_htg);
  perror("Error in writing vtkhdf file.");
  abort();
}

macro
ON_ERROR_VTK_HDF_HTG(herr_t result, vtkHDFHyperTreeGrid* vtk_hdf_htg)
{
  if (result <= H5I_INVALID_HID) {
    vtk_HDF_hypertreegrid_error(vtk_hdf_htg);
  }
}

macro
ON_ERROR_OBJ_ID_VTK_HDF_HTG(hid_t obj_id, vtkHDFHyperTreeGrid* vtk_hdf_htg)
{
  if (obj_id < 0) {
    vtk_HDF_hypertreegrid_error(vtk_hdf_htg);
  }
}

/**
 * @brief Write new attribute
 *
 * @param attribute_name The name of the attrbiute
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the attribute in
 * @param dims The dimensions of the data
 * @param vtk_hdf_htg Pointer or reference to vtkHDFHyperTreeGrid object/struct
 *
 * @memberof vtkHDFHyperTreeGrids
 */
void
vtk_HDF_hypertreegrid_write_attribute(const char* attribute_name,
                                      const void* data,
                                      const hid_t dtype_id,
                                      const hid_t group_id,
                                      const hsize_t dims[],
                                      vtkHDFHyperTreeGrid* vtk_hdf_htg)
{
  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the attribute space
  vtk_hdf_htg->attr_space_id = H5Screate_simple(1, dims, dims);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->attr_space_id, vtk_hdf_htg);

  // Set the datatype
  vtk_hdf_htg->attr_dtype_id = H5Tcopy(dtype_id);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->attr_dtype_id, vtk_hdf_htg);

  // Create the attribute
  vtk_hdf_htg->attr_id = H5Acreate2(group_id,
                                    attribute_name,
                                    vtk_hdf_htg->attr_dtype_id,
                                    vtk_hdf_htg->attr_space_id,
                                    H5P_DEFAULT,
                                    H5P_DEFAULT);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->attr_id, vtk_hdf_htg);

  // Write the attribute
  result = H5Awrite(vtk_hdf_htg->attr_id, vtk_hdf_htg->attr_dtype_id, data);
  ON_ERROR_VTK_HDF_HTG(result, vtk_hdf_htg);

  // Close the open objects
  H5Aclose(vtk_hdf_htg->attr_id);
  H5Tclose(vtk_hdf_htg->attr_dtype_id);
  H5Sclose(vtk_hdf_htg->attr_space_id);
}

/**
 * @brief Write new unchunked dataset
 *
 * @param dataset_name The name of the dataset
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the data in
 * @param rank
 * @param dims The dimensions of the dataset
 * @param vtk_hdf_htg Pointer or reference to vtkHDFHyperTreeGrid object/struct
 *
 * @memberof vtkHDFHyperTreeGrids
 */
void
vtk_HDF_hypertreegrid_write_dataset(const char* dataset_name,
                                    const void* data,
                                    const hid_t dtype_id,
                                    const hid_t group_id,
                                    const int rank,
                                    const hsize_t dims[],
                                    vtkHDFHyperTreeGrid* vtk_hdf_htg)
{
  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the data space with dimensions and maximum dimensions
  vtk_hdf_htg->dset_space_id = H5Screate_simple(rank, dims, dims);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->dset_space_id, vtk_hdf_htg);

  // Create a property list for this dataset
  vtk_hdf_htg->dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->dcpl_id, vtk_hdf_htg);

  // Set the datatype
  vtk_hdf_htg->dset_dtype_id = H5Tcopy(dtype_id);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->dset_dtype_id, vtk_hdf_htg);

  // Create the dataset
  vtk_hdf_htg->dset_id = H5Dcreate2(group_id,
                                    dataset_name,
                                    vtk_hdf_htg->dset_dtype_id,
                                    vtk_hdf_htg->dset_space_id,
                                    H5P_DEFAULT,
                                    vtk_hdf_htg->dcpl_id,
                                    H5P_DEFAULT);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->dset_id, vtk_hdf_htg);

  result = H5Dwrite(vtk_hdf_htg->dset_id,
                    vtk_hdf_htg->dset_dtype_id,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    data);
  ON_ERROR_VTK_HDF_HTG(result, vtk_hdf_htg);

  // Close opened objects
  H5Dclose(vtk_hdf_htg->dset_id);
  H5Tclose(vtk_hdf_htg->dset_dtype_id);
  H5Pclose(vtk_hdf_htg->dcpl_id);
  H5Sclose(vtk_hdf_htg->dset_space_id);
}

/**
 * @brief Write new chunked dataset
 *
 * @param dataset_name The name of the dataset
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the data in
 * @param rank
 * @param dims The dimensions of the dataset
 * @param max_dims The maximum dimensions of the dataset
 * @param chunk_dims The chunk dimensions
 * @param vtk_hdf_htg Pointer or reference to vtkHDFHyperTreeGrid object/struct
 */
void
vtk_HDF_hypertreegrid_write_chunked_dataset(const char* dataset_name,
                                            const void* data,
                                            const hid_t dtype_id,
                                            const hid_t group_id,
                                            const int rank,
                                            const hsize_t dims[],
                                            const hsize_t max_dims[],
                                            const hsize_t chunk_dims[],
                                            vtkHDFHyperTreeGrid* vtk_hdf_htg)
{

  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the data space with dimensions and maximum dimensions
  vtk_hdf_htg->dset_space_id = H5Screate_simple(rank, dims, max_dims);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->dset_space_id, vtk_hdf_htg);

  // Create a property list for this dataset
  vtk_hdf_htg->dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->dcpl_id, vtk_hdf_htg);

  // Set the chunk size for the dataset
  result = H5Pset_chunk(vtk_hdf_htg->dcpl_id, rank, chunk_dims);
  ON_ERROR_VTK_HDF_HTG(result, vtk_hdf_htg);

  // Set the datatype
  vtk_hdf_htg->dset_dtype_id = H5Tcopy(dtype_id);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->dset_dtype_id, vtk_hdf_htg);

  // Create the dataset
  vtk_hdf_htg->dset_id = H5Dcreate2(group_id,
                                    dataset_name,
                                    vtk_hdf_htg->dset_dtype_id,
                                    vtk_hdf_htg->dset_space_id,
                                    H5P_DEFAULT,
                                    vtk_hdf_htg->dcpl_id,
                                    H5P_DEFAULT);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->dset_id, vtk_hdf_htg);

  result = H5Dwrite(vtk_hdf_htg->dset_id,
                    vtk_hdf_htg->dset_dtype_id,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    data);
  ON_ERROR_VTK_HDF_HTG(result, vtk_hdf_htg);

  // Close opened objects
  H5Dclose(vtk_hdf_htg->dset_id);
  H5Tclose(vtk_hdf_htg->dset_dtype_id);
  H5Pclose(vtk_hdf_htg->dcpl_id);
  H5Sclose(vtk_hdf_htg->dset_space_id);
}

/**
 * @brief Write new chunked and compressed dataset
 *
 * @param dataset_name The name of the dataset
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the data in
 * @param rank
 * @param dims The dimensions of the dataset
 * @param max_dims The maximum dimensions of the dataset
 * @param chunk_dims The chunk dimensions
 * @param compression_level The level of compression
 * @param vtk_hdf_htg Pointer or reference to vtkHDFHyperTreeGrid object/struct
 *
 * @memberof vtkHDFHyperTreeGrid
 */
void
vtk_HDF_hypertreegrid_write_compressed_dataset(const char* dataset_name,
                                               const void* data,
                                               const hid_t dtype_id,
                                               const hid_t group_id,
                                               const int rank,
                                               const hsize_t dims[],
                                               const hsize_t max_dims[],
                                               const hsize_t chunk_dims[],
                                               unsigned int compression_level,
                                               vtkHDFHyperTreeGrid* vtk_hdf_htg)
{

  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the data space with dimensions and maximum dimensions
  vtk_hdf_htg->dset_space_id = H5Screate_simple(rank, dims, max_dims);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->dset_space_id, vtk_hdf_htg);

  // Create a property list for this dataset
  vtk_hdf_htg->dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->dcpl_id, vtk_hdf_htg);

  // Set the chunk size for the dataset
  result = H5Pset_chunk(vtk_hdf_htg->dcpl_id, rank, chunk_dims);
  ON_ERROR_VTK_HDF_HTG(result, vtk_hdf_htg);

  // Set compression on the dataset
  result = H5Pset_deflate(vtk_hdf_htg->dcpl_id, compression_level);
  ON_ERROR_VTK_HDF_HTG(result, vtk_hdf_htg);

  // Set the datatype
  vtk_hdf_htg->dset_dtype_id = H5Tcopy(dtype_id);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->dset_dtype_id, vtk_hdf_htg);

  // Create the dataset
  vtk_hdf_htg->dset_id = H5Dcreate2(group_id,
                                    dataset_name,
                                    vtk_hdf_htg->dset_dtype_id,
                                    vtk_hdf_htg->dset_space_id,
                                    H5P_DEFAULT,
                                    vtk_hdf_htg->dcpl_id,
                                    H5P_DEFAULT);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->dset_id, vtk_hdf_htg);

  result = H5Dwrite(vtk_hdf_htg->dset_id,
                    vtk_hdf_htg->dset_dtype_id,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    data);
  ON_ERROR_VTK_HDF_HTG(result, vtk_hdf_htg);

  // Close opened objects
  H5Dclose(vtk_hdf_htg->dset_id);
  H5Tclose(vtk_hdf_htg->dset_dtype_id);
  H5Pclose(vtk_hdf_htg->dcpl_id);
  H5Sclose(vtk_hdf_htg->dset_space_id);
}

/**
 * @brief Write a collectively written dataset with MPI-IO
 *
 * @param dataset_name The name of the dataset
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the data in
 * @param rank
 * @param dims The dimensions of the dataset
 * @param max_dims The maximum dimensions of the dataset
 * @param local_size The size of the sub-array this process writes into
 * @param local_offset The position of the sub-array this process writes into
 * @param vtk_hdf_htg Pointer or reference to vtkHDFHyperTreeGrid object/struct
 *
 * @memberof vtkHDFHyperTreeGrid
 */
void
vtk_HDF_hypertreegrid_collective_write_dataset(const char* dataset_name,
                                               const void* data,
                                               const hid_t dtype_id,
                                               const hid_t group_id,
                                               const int rank,
                                               const hsize_t dims[],
                                               const hsize_t local_size[],
                                               const hsize_t local_offset[],
                                               vtkHDFHyperTreeGrid* vtk_hdf_htg)
{

  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the data space with dimensions and maximum dimensions
  vtk_hdf_htg->dset_space_id = H5Screate_simple(rank, dims, dims);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->dset_space_id, vtk_hdf_htg);

  // Set the datatype
  vtk_hdf_htg->dset_dtype_id = H5Tcopy(dtype_id);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->dset_dtype_id, vtk_hdf_htg);

  // Create the dataset
  vtk_hdf_htg->dset_id = H5Dcreate2(group_id,
                                    dataset_name,
                                    vtk_hdf_htg->dset_dtype_id,
                                    vtk_hdf_htg->dset_space_id,
                                    H5P_DEFAULT,
                                    H5P_DEFAULT,
                                    H5P_DEFAULT);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->dset_id, vtk_hdf_htg);

  // Create a property list for MPI-IO transfer
  vtk_hdf_htg->xfer_plist = H5Pcreate(H5P_DATASET_XFER);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->xfer_plist, vtk_hdf_htg);

  // Set MPIO_COLLECTIVE writing
  result = H5Pset_dxpl_mpio(vtk_hdf_htg->xfer_plist, H5FD_MPIO_COLLECTIVE);
  ON_ERROR_VTK_HDF_HTG(result, vtk_hdf_htg);

  // Get the dataset space
  vtk_hdf_htg->file_space = H5Dget_space(vtk_hdf_htg->dset_id);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->file_space, vtk_hdf_htg);

  // Select a hyperslab for our process
  result = H5Sselect_hyperslab(vtk_hdf_htg->file_space,
                               H5S_SELECT_SET,
                               local_offset,
                               NULL,
                               local_size,
                               NULL);
  ON_ERROR_VTK_HDF_HTG(result, vtk_hdf_htg);

  // Create a memory dataspace for our process
  vtk_hdf_htg->mem_space = H5Screate_simple(rank, local_size, NULL);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->mem_space, vtk_hdf_htg);

  // Actually do the collective write to the file
  result = H5Dwrite(
    vtk_hdf_htg->dset_id,       /* dataset handle */
    vtk_hdf_htg->dset_dtype_id, /* H5T_IEEE_F64LE */
    vtk_hdf_htg->mem_space,     /* memory dataspace [local_nx] */
    vtk_hdf_htg->file_space,    /* file dataspace with hyperslab selected */
    vtk_hdf_htg->xfer_plist,    /* collective MPI‐IO transfer property */
    data                        /* pointer to local data */
  );
  ON_ERROR_VTK_HDF_HTG(result, vtk_hdf_htg);

  // Close opened objects
  H5Pclose(vtk_hdf_htg->xfer_plist);
  vtk_hdf_htg->xfer_plist = H5I_INVALID_HID;
  H5Sclose(vtk_hdf_htg->mem_space);
  vtk_hdf_htg->mem_space = H5I_INVALID_HID;
  H5Sclose(vtk_hdf_htg->file_space);
  vtk_hdf_htg->file_space = H5I_INVALID_HID;
  H5Tclose(vtk_hdf_htg->dset_dtype_id);
  H5Sclose(vtk_hdf_htg->dset_space_id);
  H5Dclose(vtk_hdf_htg->dset_id);
}

/**
 * @brief Write a chunked dataset collectively using MPI-IO
 *
 * @param dataset_name The name of the dataset
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the data in
 * @param dims The dimensions of the dataset
 * @param rank
 * @param max_dims The maximum dimensions of the dataset
 * @param chunk_dims The chunk dimensions
 * @param local_size The size of the sub-array this process writes into
 * @param local_offset The position of the sub-array this process writes into
 * @param vtk_hdf_htg Pointer or reference to vtkHDFHyperTreeGrid object/struct
 *
 * @memberof vtkHDFHyperTreeGrid
 */
void
vtk_HDF_hypertreegrid_collective_write_chunked_dataset(
  const char* dataset_name,
  const void* data,
  const hid_t dtype_id,
  const hid_t group_id,
  const int rank,
  const hsize_t dims[],
  const hsize_t max_dims[],
  const hsize_t chunk_dims[],
  const hsize_t local_size[],
  const hsize_t local_offset[],
  vtkHDFHyperTreeGrid* vtk_hdf_htg)
{

  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the data space with dimensions and maximum dimensions
  vtk_hdf_htg->dset_space_id = H5Screate_simple(rank, dims, max_dims);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->dset_space_id, vtk_hdf_htg);

  // Create a property list for this dataset
  vtk_hdf_htg->dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->dcpl_id, vtk_hdf_htg);

  // Set the chunk size for the dataset
  result = H5Pset_chunk(vtk_hdf_htg->dcpl_id, rank, chunk_dims);
  ON_ERROR_VTK_HDF_HTG(result, vtk_hdf_htg);

  // Set the datatype
  vtk_hdf_htg->dset_dtype_id = H5Tcopy(dtype_id);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->dset_dtype_id, vtk_hdf_htg);

  // Create the dataset
  vtk_hdf_htg->dset_id = H5Dcreate2(group_id,
                                    dataset_name,
                                    vtk_hdf_htg->dset_dtype_id,
                                    vtk_hdf_htg->dset_space_id,
                                    H5P_DEFAULT,
                                    vtk_hdf_htg->dcpl_id,
                                    H5P_DEFAULT);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->dset_id, vtk_hdf_htg);

  // Create a property list for MPI-IO transfer
  vtk_hdf_htg->xfer_plist = H5Pcreate(H5P_DATASET_XFER);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->xfer_plist, vtk_hdf_htg);

  // Set MPIO_COLLECTIVE writing
  result = H5Pset_dxpl_mpio(vtk_hdf_htg->xfer_plist, H5FD_MPIO_COLLECTIVE);
  ON_ERROR_VTK_HDF_HTG(result, vtk_hdf_htg);

  // Get the dataset space
  vtk_hdf_htg->file_space = H5Dget_space(vtk_hdf_htg->dset_id);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->file_space, vtk_hdf_htg);

  // Select a hyperslab for our process
  result = H5Sselect_hyperslab(vtk_hdf_htg->file_space,
                               H5S_SELECT_SET,
                               local_offset,
                               NULL,
                               local_size,
                               NULL);
  ON_ERROR_VTK_HDF_HTG(result, vtk_hdf_htg);

  // Create a memory dataspace for our process
  vtk_hdf_htg->mem_space = H5Screate_simple(rank, local_size, NULL);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->mem_space, vtk_hdf_htg);

  // Actually do the collective write to the file
  result = H5Dwrite(
    vtk_hdf_htg->dset_id,       /* dataset handle */
    vtk_hdf_htg->dset_dtype_id, /* H5T_IEEE_F64LE */
    vtk_hdf_htg->mem_space,     /* memory dataspace */
    vtk_hdf_htg->file_space,    /* file dataspace with hyperslab selected */
    vtk_hdf_htg->xfer_plist,    /* collective MPI‐IO transfer property */
    data                        /* pointer to local data */
  );
  ON_ERROR_VTK_HDF_HTG(result, vtk_hdf_htg);

  // Close opened objects
  H5Pclose(vtk_hdf_htg->xfer_plist);
  vtk_hdf_htg->xfer_plist = H5I_INVALID_HID;
  H5Sclose(vtk_hdf_htg->mem_space);
  vtk_hdf_htg->mem_space = H5I_INVALID_HID;
  H5Sclose(vtk_hdf_htg->file_space);
  vtk_hdf_htg->file_space = H5I_INVALID_HID;
  H5Tclose(vtk_hdf_htg->dset_dtype_id);
  H5Pclose(vtk_hdf_htg->dcpl_id);
  H5Sclose(vtk_hdf_htg->dset_space_id);
  H5Dclose(vtk_hdf_htg->dset_id);
}

/**
 * @brief Write a collectively written dataset using compression with MPI-IO
 *
 * @param dataset_name The name of the dataset
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the data in
 * @param rank
 * @param dims The dimensions of the dataset
 * @param max_dims The maximum dimensions of the dataset
 * @param chunk_dims The chunk dimensions
 * @param local_size The size of the sub-array this process writes into
 * @param local_offset The position of the sub-array this process writes into
 * @param compression_level The level of compression
 * @param vtk_hdf_htg Pointer or reference to vtkHDFHyperTreeGrid object/struct
 *
 * @memberof vtkHDFHyperTreeGrid
 */
void
vtk_HDF_hypertreegrid_collective_write_compressed_dataset(
  const char* dataset_name,
  const void* data,
  const hid_t dtype_id,
  const hid_t group_id,
  const int rank,
  const hsize_t dims[],
  const hsize_t max_dims[],
  const hsize_t chunk_dims[],
  const hsize_t local_size[],
  const hsize_t local_offset[],
  unsigned int compression_level,
  vtkHDFHyperTreeGrid* vtk_hdf_htg)
{
  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the data space with dimensions and maximum dimensions
  vtk_hdf_htg->dset_space_id = H5Screate_simple(rank, dims, max_dims);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->dset_space_id, vtk_hdf_htg);

  // Create a property list for this dataset
  vtk_hdf_htg->dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->dcpl_id, vtk_hdf_htg);

  // Set the chunk size for the dataset
  result = H5Pset_chunk(vtk_hdf_htg->dcpl_id, rank, chunk_dims);
  ON_ERROR_VTK_HDF_HTG(result, vtk_hdf_htg);

  // Set compression on the dataset
  result = H5Pset_deflate(vtk_hdf_htg->dcpl_id, compression_level);
  ON_ERROR_VTK_HDF_HTG(result, vtk_hdf_htg);

  // Set the datatype
  vtk_hdf_htg->dset_dtype_id = H5Tcopy(dtype_id);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->dset_dtype_id, vtk_hdf_htg);

  // Create the dataset
  vtk_hdf_htg->dset_id = H5Dcreate2(group_id,
                                    dataset_name,
                                    vtk_hdf_htg->dset_dtype_id,
                                    vtk_hdf_htg->dset_space_id,
                                    H5P_DEFAULT,
                                    vtk_hdf_htg->dcpl_id,
                                    H5P_DEFAULT);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->dset_id, vtk_hdf_htg);

  // Create a property list for MPI-IO transfer
  vtk_hdf_htg->xfer_plist = H5Pcreate(H5P_DATASET_XFER);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->xfer_plist, vtk_hdf_htg);

  // Set MPIO_COLLECTIVE writing
  result = H5Pset_dxpl_mpio(vtk_hdf_htg->xfer_plist, H5FD_MPIO_COLLECTIVE);
  ON_ERROR_VTK_HDF_HTG(result, vtk_hdf_htg);

  // Get the dataset space
  vtk_hdf_htg->file_space = H5Dget_space(vtk_hdf_htg->dset_id);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->file_space, vtk_hdf_htg);

  // Select a hyperslab for our process
  result = H5Sselect_hyperslab(vtk_hdf_htg->file_space,
                               H5S_SELECT_SET,
                               local_offset,
                               NULL,
                               local_size,
                               NULL);
  ON_ERROR_VTK_HDF_HTG(result, vtk_hdf_htg);

  // Create a memory dataspace for our process
  vtk_hdf_htg->mem_space = H5Screate_simple(rank, local_size, NULL);
  ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg->mem_space, vtk_hdf_htg);

  // Actually do the collective write to the file
  result = H5Dwrite(
    vtk_hdf_htg->dset_id,       /* dataset handle */
    vtk_hdf_htg->dset_dtype_id, /* H5T_IEEE_F64LE */
    vtk_hdf_htg->mem_space,     /* memory dataspace [local_nx] */
    vtk_hdf_htg->file_space,    /* file dataspace with hyperslab selected */
    vtk_hdf_htg->xfer_plist,    /* collective MPI‐IO transfer property */
    data                        /* pointer to local data */
  );
  ON_ERROR_VTK_HDF_HTG(result, vtk_hdf_htg);

  // Close opened objects
  H5Pclose(vtk_hdf_htg->xfer_plist);
  vtk_hdf_htg->xfer_plist = H5I_INVALID_HID;
  H5Sclose(vtk_hdf_htg->mem_space);
  vtk_hdf_htg->mem_space = H5I_INVALID_HID;
  H5Sclose(vtk_hdf_htg->file_space);
  vtk_hdf_htg->file_space = H5I_INVALID_HID;
  H5Tclose(vtk_hdf_htg->dset_dtype_id);
  H5Pclose(vtk_hdf_htg->dcpl_id);
  H5Sclose(vtk_hdf_htg->dset_space_id);
  H5Dclose(vtk_hdf_htg->dset_id);
}

/**
 * @brief This function does the actual writing of the
 * HyperTreeGrid/PHyperTreeGrid in the HDF5 format.
 *
 * @param scalar_list List of any basilisk scalar fields
 * @param vector_list List of any basilisk vector fields
 * @param fname The filename to write the vtkhdf (HDF5) file to
 *
 * @memberof vtkHDFHyperTreeGrid
 *
 * The HyperTreeGrid in [VTK](https://vtk.org/) is a
 * [class](https://vtk.org/doc/nightly/html/classvtkHyperTreeGrid.html#details)
 * that is intended exactly for representation of tree-based grids in VTK. It is
 * similar to Basilisk, although the major difference is that while the
 * foreach_cell() iterator in Basilisk uses [depth-first
 * search](https://en.wikipedia.org/wiki/Depth-first_search), the HyperTreeGrid
 * format requires us to describe the tree using [breadth-first
 * search](https://en.wikipedia.org/wiki/Breadth-first_search).
 *
 * Some other differences worth noting is that the HyperTreeGrid also supports
 * multiple trees in a regular grid, even though basilisk does not support this.
 * Thus, there are some seemingly extraneous fields included that we do not make
 * use of, though are neccesary if we wanted to describe such a "forest" of
 * trees.
 *
 * It is also worth noting that the HDF5 format itself is somewhat generic (like
 * XML), even if VTK imposes a specific schema for HyperTreeGrid formats. In
 * HDF5, we have two basic building blocks: Dataset and Group objects. These are
 * self-descriptive. In addition, each may have "attributes" attached to them,
 * which themselves are small pieces of data used to help interpret and
 * understand the Group and Dataset objects. In addition, it is worthwhile to
 * understand HDF5 chunking, since this has some implications for writing
 * appendable datasets (H5S_UNLIMITED) or when using compression features of
 * HDF5.
 *
 * In regards to what collection of Groups, Datasets, and Attributes in HDF5
 * define a HyperTreeGrid, the only official documentation on the schema from
 * Kitware is found
 * [here](https://web.archive.org/web/20250804033704/https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html#hypertreegrid).
 * Note that even though the format supports the creation of temporal datasets,
 * we write one separate file per time step. This is to avoid the risk of
 * corruption. In addition, this writer was made for version 2.4 (of the file
 * format). The format
 *
 * In summary, our schema in the case of a non-temporal (single-timestep)
 * HyperTreeGrid is as follows
 *
 * - **Group**: "/VTKHDF"
 *
 *  - **Attribute**: "Version"
 *    - **Datatype**: `int64`
 *    - **Dimension**: `{2}`
 *    - **Description**: This is the VTKHDF file 'version' we use: The first
 * entry is the major version and the second entry represents the minor version.
 *
 *  - **Attribute**: "Type"
 *    - **Datatype**: `char`
 *    - **Dimension**: `{13}`
 *    - **Description**: This tells VTK what format we are using (HyperTreeGrid)
 * as opposed to ImageGrid, UnstructuredGrid, etc.
 *
 *  - **Attribute**: "BranchFactor"
 *    - **Datatype**: `int64`
 *    - **Dimension**: `{1}`
 *    - **Description**: When a tree cell is refined, it should have n^d
 * children, where d is the dimension and n is the branch factor. For basilisk,
 * this is always 2.
 *
 *  - **Attribute**: "Dimensions"
 *    - **Datatype**: `int64`
 *    - **Dimension**: `{3}`
 *    - **Description**: The number of coordinates describing the trees in each
 * dimension; this should match/describe the number of entries in the
 * XCoordinate, YCoordinate, and ZCoordinate datasets, when the tree is not a
 * temporal tree or PHyperTree.
 *
 *  - **Attribute**: "TransposedRootIndexing"
 *    - **Datatype**: `int64`
 *    - **Dimension**: `{1}`
 *    - **Value**: `{0}`
 *    - **Description**: I am unsure what the significance of this one is, but
 * we do not seem to need this feature.
 *
 *  - **Dataset**: "DepthPerTree"
 *    - **Datatype**: `int64`
 *    - **Dimension**: `{1}`
 *    - **Description**: The maximum depth of the tree. In the case of basilisk
 * we only have one.
 *
 *  - **Dataset**: "Descriptors"
 *    - **Datatype**: `char`
 *    - **Description**: This is the description of the actual structure of the
 * tree by breadth-first search. \sa @ref hdf_get_descriptors
 *
 *  - **Dataset**: "DescriptorsSize"
 *    - **Datatype**: `int64`
 *    - **Dimension**: `{1}`
 *    - **Description**: The size of the descriptors dataset. This only becomes
 * useful when there are a forest of trees (or parallel copies of the same
 * tree), so that the VTK reader knows where the descriptors of each tree starts
 * and ends in the concatenation of all descriptors.
 *
 *  - **Dataset**: "NumberOfCells"
 *    - **Datatype**: `int64`
 *    - **Dimension**: {1}
 *    - **Description**: Cell count across the entire tree, for each tree. Since
 * we have only one tree in basilisk, this has dimension 1.
 *
 *  - **Dataset**: "NumberOfCellsPerDepth"
 *    - **Datatype**: `int64`
 *    - **Dimension**: {NumberOfDepths}
 *    - **Description**: Cell count on each depth level of the tree.
 *
 *  - **Dataset**: "NumberOfDepths"
 *    - **Datatype**: `int64`
 *    - **Dimension**: `{1}`
 *    - **Description**: I am not sure what the significance of this is, but we
 * just write the same value(s) as DepthPerTree, which is the maximum depth of
 * the tree.
 *
 *  - **Dataset**: "NumberOfTrees"
 *    - **Datatype**: `int64`
 *    - **Dimension**: `{}`
 *    - **Value**: `{1}`
 *    - **Description**: The total number of trees. In basilisk we have just
 * one.
 *
 *  - **Dataset**: "Mask"
 *    - **Datatype**: `char`
 *    - **Dimension**: `{}`
 *    - **Description**: This field is entirely optional. If it is included, VTK
 * will hide or show certain cells. The format is nearly identical to that of
 * the descriptors dataset.
 *
 *  - **Dataset**: "TreeIds"
 *    - **Datatype**: `int64`
 *    - **Dimension**: `{1}`
 *    - **Value**: `{0}`
 *    - **Description**: There is only one tree in basilisk, and the ID of the
 * first tree in a HyperTreeGrid should be zero. If we had more than one tree,
 * the placement/ordering of the trees into the regular grid (forest) is
 * determined by the IDs given here.
 *
 *  - **Dataset**: "XCoordinates"
 *    - **Datatype**: `double`
 *    - **Dimension**: `{2}`
 *    - **Description**: The x coordinate of two corners of our tree.
 *
 *  - **Dataset**: "YCoordinates"
 *    - **Datatype**: `double`
 *    - **Dimension**: `{2} or {1}`
 *    - **Description**: The y coordinate of two corners of our tree. If the
 * tree is really a binary tree, just write 0 as the single entry.
 *
 *  - **Dataset**: "ZCoordinates"
 *    - **Datatype**: `double`
 *    - **Dimension**: `{2}` or `{1}`
 *    - **Description**: The z coordinate of two corners of our tree. If the
 * tree is really a quad tree, just write 0 as the single entry.
 *
 *  - **Group**: "/VTKHDF/CellData"
 *
 *    - **Dataset**: "MyCellCenteredScalarField"
 *      - **DataType**: `float` or `double`
 *      - **Dimension**: `{}`
 *    - **Dataset**: "MyCellCenteredVectorField"
 *      - **DataType**: `float` or `double`
 *      - **Dimension**: `{,d}`
 *
 * If done correctly, we should be able to run `h5dump -H` on our file and get
 * something like
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  HDF5 "myfile.vtkhdf" {
 *  GROUP "/" {
 *     GROUP "VTKHDF" {
 *        ATTRIBUTE "BranchFactor" {
 *           DATATYPE  H5T_STD_I64LE
 *           DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
 *        }
 *        ATTRIBUTE "Dimensions" {
 *           DATATYPE  H5T_STD_I64LE
 *           DATASPACE  SIMPLE { ( 3 ) / ( 3 ) }
 *        }
 *        ATTRIBUTE "TransposedRootIndexing" {
 *           DATATYPE  H5T_STD_I64LE
 *           DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
 *        }
 *        ATTRIBUTE "Type" {
 *           DATATYPE  H5T_STRING {
 *              STRSIZE 13;
 *              STRPAD H5T_STR_NULLPAD;
 *              CSET H5T_CSET_ASCII;
 *              CTYPE H5T_C_S1;
 *           }
 *           DATASPACE  SCALAR
 *        }
 *        ATTRIBUTE "Version" {
 *           DATATYPE  H5T_STD_I64LE
 *           DATASPACE  SIMPLE { ( 2 ) / ( 2 ) }
 *        }
 *        GROUP "CellData" {
 *           DATASET "p" {
 *              DATATYPE  H5T_IEEE_F32LE
 *              DATASPACE  SIMPLE { ( 87381 ) / ( 87381 ) }
 *           }
 *           DATASET "u" {
 *              DATATYPE  H5T_IEEE_F32LE
 *              DATASPACE  SIMPLE { ( 87381, 2 ) / ( 87381, 2 ) }
 *           }
 *        }
 *        DATASET "DepthPerTree" {
 *           DATATYPE  H5T_STD_I64LE
 *           DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
 *        }
 *        DATASET "Descriptors" {
 *           DATATYPE  H5T_STD_U8LE
 *           DATASPACE  SIMPLE { ( 2731 ) / ( 2731 ) }
 *        }
 *        DATASET "DescriptorsSize" {
 *           DATATYPE  H5T_STD_I64LE
 *           DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
 *        }
 *        DATASET "NumberOfCells" {
 *           DATATYPE  H5T_STD_I64LE
 *           DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
 *        }
 *        DATASET "NumberOfCellsPerTreeDepth" {
 *           DATATYPE  H5T_STD_I64LE
 *           DATASPACE  SIMPLE { ( 9 ) / ( 9 ) }
 *        }
 *        DATASET "NumberOfDepths" {
 *           DATATYPE  H5T_STD_I64LE
 *           DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
 *        }
 *        DATASET "NumberOfTrees" {
 *           DATATYPE  H5T_STD_I64LE
 *           DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
 *        }
 *        DATASET "TreeIds" {
 *           DATATYPE  H5T_STD_I64LE
 *           DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
 *        }
 *        DATASET "XCoordinates" {
 *           DATATYPE  H5T_IEEE_F64LE
 *           DATASPACE  SIMPLE { ( 2 ) / ( 2 ) }
 *        }
 *        DATASET "YCoordinates" {
 *           DATATYPE  H5T_IEEE_F64LE
 *           DATASPACE  SIMPLE { ( 2 ) / ( 2 ) }
 *        }
 *        DATASET "ZCoordinates" {
 *           DATATYPE  H5T_IEEE_F64LE
 *           DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
 *        }
 *     }
 *  }
 *  }
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * For the non-temporal parallel "PHyperTreeGrid", the schema is exactly
 * identical, although we now write multiple trees over the same X/Y/Z
 * coordinates and (for the most part) VTK will stitch the trees from each
 * process together in a coherent way. This requires some MPI coordination on
 * our part to know the offsets into each dataset--though otherwise, the HDF5
 * library handles all other details of MPIIO coordination.
 *
 *
 *
 */
vtkHDFHyperTreeGrid
vtk_HDF_hypertreegrid_init(scalar* scalar_list,
                           vector* vector_list,
                           const char* fname)
{

  // Create the vtkHDF struct
  vtkHDF vtk_hdf = vtk_HDF_init(fname);

  // Initialize our object ids with -1 or H5I_INVALID_HID and the vtkHDF struct
  vtkHDFHyperTreeGrid vtk_hdf_htg = {
    .vtk_hdf = vtk_hdf,
    .grp_celldata_id = H5I_INVALID_HID,
    .attr_space_id = H5I_INVALID_HID,
    .attr_dtype_id = H5I_INVALID_HID,
    .attr_id = H5I_INVALID_HID,
    .dset_space_id = H5I_INVALID_HID,
    .dcpl_id = H5I_INVALID_HID,
    .dset_dtype_id = H5I_INVALID_HID,
    .dset_id = H5I_INVALID_HID,
    .file_space = H5I_INVALID_HID,
    .mem_space = H5I_INVALID_HID,
    .xfer_plist = H5I_INVALID_HID,
  };

  /*
   * Attribute: /VTKHDF/BranchFactor
   * Datatype: H5T_NATIVE_INT64 / int64_t
   * Dimension: {1}
   *
   * This attribute describes the number of children of a refined cell. This is
   * \f(n^d\f) where \f(n\f) is the branch factor and \f(d\f) is the dimension
   * of the tree. In basilisk \f(n=2\f) always.
   */
  {
    // Value for the BranchFactor attrbiute
    const int64_t bf_value = 2;

    // Dimensions of the attribute
    hsize_t dims_attr[1] = { 1 };

    vtk_HDF_hypertreegrid_write_attribute(
      "BranchFactor",                    /* attribute_name */
      &bf_value,                         /* data */
      H5T_NATIVE_INT64,                  /* dtype_id */
      vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, /* group_id */
      dims_attr,                         /* dims */
      &vtk_hdf_htg                       /* vtkHDFHyperTreeGrid object/struct */
    );
  }

  /*
   * Attribute: /VTKHDF/Dimensions
   * Datatype: H5T_NATIVE_INT64 / int64_t
   * Dimension: {3}
   *
   * This attribute describes the number of coordinates in each cartesian X,Y,Z
   * direction. For basilisk, each will never be greater than 2, since we can
   * only have 1 tree.
   */
  {
    // Set the value of this attribute, depending on the compile-time definition
    // of dimension in basilisk
#if dimension == 1
    const int64_t dims_value[3] = { 2, 1, 1 };
#elif dimension == 2
    const int64_t dims_value[3] = { 2, 2, 1 };
#else
    const int64_t dims_value[3] = { 2, 2, 2 };
#endif
    // Dimensions of the attribute
    const hsize_t dims_attr[1] = { 3 };

    vtk_HDF_hypertreegrid_write_attribute(
      "Dimensions",                      /* attribute_name */
      dims_value,                        /* data */
      H5T_NATIVE_INT64,                  /* dtype_id */
      vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, /* group_id */
      dims_attr,                         /* dims */
      &vtk_hdf_htg                       /* vtkHDFHyperTreeGrid object/struct */
    );
  }

  /*
   * Attribute: /VTKHDF/TransposedRootIndexing
   * Datatype: H5T_NATIVE_INT64 / int64_t
   * Dimension: {1}
   */
  {
    // Value for the TransposedRootIndexing attribute
    const int64_t tri_value = 0;

    // Dimensions for the attribute
    const hsize_t dims_attr[1] = { 1 };

    vtk_HDF_hypertreegrid_write_attribute(
      "TransposedRootIndexing",          /* attribute_name */
      &tri_value,                        /* data */
      H5T_NATIVE_INT64,                  /* dtype_id */
      vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, /* group_id */
      dims_attr,                         /* dims */
      &vtk_hdf_htg                       /* vtkHDFHyperTreeGrid object/struct */
    );
  }

  /*
   * Attribute: /VTKHDF/Type
   * Datatype: H5T_C_S1 / char
   * Dimension: {1}
   *
   * This should be "HyperTreeGrid" to distinguish from other supported types
   * (e.g. "ImageData"). Note that this is written into a scalar data space
   * (rather than an data array with dimensions) using a string type with (at
   * most) length 13
   */
  {
    // Store the result of HDF5 functions that do not return an identifier
    herr_t result = 0;

    // Value of the Type attribute
    const char* type_str = "HyperTreeGrid";

    // Create the attribute space
    vtk_hdf_htg.attr_space_id = H5Screate(H5S_SCALAR);
    ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg.attr_space_id, &vtk_hdf_htg);

    // Create a fixed-length string datatype of length 13, null-padded, ASCII
    vtk_hdf_htg.attr_dtype_id = H5Tcopy(H5T_C_S1);
    ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg.attr_dtype_id, &vtk_hdf_htg);

    // Set the size
    result = H5Tset_size(vtk_hdf_htg.attr_dtype_id, (size_t)13);
    ON_ERROR_VTK_HDF_HTG(result, &vtk_hdf_htg);

    // Pad the string
    result = H5Tset_strpad(vtk_hdf_htg.attr_dtype_id, H5T_STR_NULLPAD);
    ON_ERROR_VTK_HDF_HTG(result, &vtk_hdf_htg);

    // Set the string value
    result = H5Tset_cset(vtk_hdf_htg.attr_dtype_id, H5T_CSET_ASCII);
    ON_ERROR_VTK_HDF_HTG(result, &vtk_hdf_htg);

    // Create the attribute
    vtk_hdf_htg.attr_id = H5Acreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id,
                                     "Type",
                                     vtk_hdf_htg.attr_dtype_id,
                                     vtk_hdf_htg.attr_space_id,
                                     H5P_DEFAULT,
                                     H5P_DEFAULT);
    ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg.attr_id, &vtk_hdf_htg);

    /* Write the string (automatically null‐padded up to length 13) */
    result = H5Awrite(vtk_hdf_htg.attr_id, vtk_hdf_htg.attr_dtype_id, type_str);
    ON_ERROR_VTK_HDF_HTG(result, &vtk_hdf_htg);

    // Close the open objects
    H5Aclose(vtk_hdf_htg.attr_id);
    H5Tclose(vtk_hdf_htg.attr_dtype_id);
    H5Sclose(vtk_hdf_htg.attr_space_id);
  }

  /*
   * Attribute: /VTKHDF/Version
   * Datatype: H5T_NATIVE_INT64 / int64_t
   * Dimension: {2}
   */
  {
    // Value of the VTKHDF version attribute
    const int64_t vers_value[2] = { 2, 4 };

    // Dimensions of the attribute
    const hsize_t dims_attr[1] = { 2 };

    vtk_HDF_hypertreegrid_write_attribute(
      "Version",                         /* attribute_name */
      vers_value,                        /* data */
      H5T_NATIVE_INT64,                  /* dtype_id */
      vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, /* group_id */
      dims_attr,                         /* dims */
      &vtk_hdf_htg                       /* vtkHDFHyperTreeGrid object/struct */
    );
  }

  /*
   * Create the vtkHDFHyperTreeGridData object: This object contains a lot of
   * local-view data that we will write after this. See @ref
   * vtk_hdf_hypertreegrid_data_init in vtkHDFHyperTreeGridData.h
   */
  vtkHDFHyperTreeGridData* vtk_hdf_htg_data = vtk_hdf_hypertreegrid_data_init();

  /*
   * Group: /VTKHDF/CellData
   *
   * Here we place all of our grid, tree, and field data. This group does not
   * require any attributes.
   */
  {
    // Create the new CellData group inside group VTKHDF
    vtk_hdf_htg.grp_celldata_id = H5Gcreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id,
                                             "CellData",
                                             H5P_DEFAULT,
                                             H5P_DEFAULT,
                                             H5P_DEFAULT);
    ON_ERROR_OBJ_ID_VTK_HDF_HTG(vtk_hdf_htg.grp_celldata_id, &vtk_hdf_htg);
  }

#if MPI_SINGLE_FILE
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  /*
   * Dataset: /VTKHDF/DepthPerTree
   * Datatype: H5T_NATIVE_INT64 / int64_t
   * Dimension: {1*npe()}
   *
   * The maximum depth of our single tree. In the MPI_SINGLE_FILE case, we must
   * write this once per proc according to the maximum (actual) depth of the
   * tree on that proc.
   */
  {
#if MPI_SINGLE_FILE
    hsize_t local_size[] = { 1 };
    hsize_t global_size[] = { npe() };
    hsize_t local_offset[] = { pid() };

    vtk_HDF_hypertreegrid_collective_write_dataset(
      "DepthPerTree",                    /* dataset_name */
      &vtk_hdf_htg_data->depth_per_tree, /* pointer to data */
      H5T_NATIVE_INT64,                  /* datatype */
      vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, /* group id */
      1,                                 /* rank */
      global_size,                       /* dimensions of dataset */
      local_size,                        /* dimensions of local size */
      local_offset, /* where in the global dataset to write */
      &vtk_hdf_htg  /* pointer to vtkHDFHyperTreeGrid object */
    );

#else
    hsize_t local_size[] = { 1 };

    vtk_HDF_hypertreegrid_write_dataset(
      "DepthPerTree",                    /* dataset_name */
      &vtk_hdf_htg_data->depth_per_tree, /* pointer to data */
      H5T_NATIVE_INT64,                  /* datatype */
      vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, /* group id */
      1,                                 /* rank */
      local_size,                        /*dimensions of dataset */
      &vtk_hdf_htg /* pointer to vtkHDFHyperTreeGrid object */
    );
#endif
  }

  /*
   * Dataset: /VTKHDF/Descriptors
   * Datatype: H5T_STD_U8LE / UInt8_t
   * Dimension: {global_size}
   *
   * The [breadth-first
   * search](https://en.wikipedia.org/wiki/Breadth-first_search) description of
   * the tree. Each bit represents a node of the tree, starting with the root.
   * If a given node is refined, then a '1' is written. If a given node is not
   * refined, then a '0' is written instead. For example, the binary tree
   *
   * [Level 0]                      o
   *                               / \
   *                              /   \
   * [Level 1]                   o     o
   *                            / \
   *                           /   \
   * [Level 2]                o     o
   *
   * would be written as "1 10 00" or "11000000" formatted as a byte. The actual
   * values of this are computed in @ref hdf_get_descriptors() in @ref
   * vtkHDFHyperTreeGridData.h and the actual breadth-first traversal (in the
   * correct order) is performed by the @ref foreach_cell_bfs() macro in @ref
   * foreach_cell_bfs.h.
   *
   * In the MPI_SINGLE_FILE case, we write this descriptor array for the tree
   * local to each proc.
   *
   * \sa @ref hdf_get_descriptors
   * \sa @ref foreach_cell_bdf
   */
  {
#if MPI_SINGLE_FILE
    // Calculate local and global sizes for collective MPI-IO operations
    hsize_t local_size = (hsize_t)vtk_hdf_htg_data->descriptors_size;
    hsize_t global_size = local_size;

    // MPI_Allreduce: Aggregate local sizes to compute total global dataset size
    // across all ranks
    MPI_Allreduce(&local_size,
                  &global_size,
                  1,
                  MPI_UNSIGNED_LONG_LONG,
                  MPI_SUM,
                  MPI_COMM_WORLD);

    // MPI_Exscan: Calculate exclusive prefix sum to determine each rank's
    // offset in the global dataset This ensures each MPI rank writes to a
    // non-overlapping region of the HDF5 data array
    hsize_t local_offset = 0;

    MPI_Exscan(&local_size,
               &local_offset,
               1,
               MPI_UNSIGNED_LONG_LONG,
               MPI_SUM,
               MPI_COMM_WORLD);
    if (pid() == 0)
      local_offset = 0;

    vtk_HDF_hypertreegrid_collective_write_dataset(
      "Descriptors",                     /* dataset_name */
      vtk_hdf_htg_data->descriptors,     /* pointer to data */
      H5T_STD_U8LE,                      /* datatype */
      vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, /* group id */
      1,                                 /* rank */

      &global_size,  /* dimensions of dataset */
      &local_size,   /* dimensions of local size */
      &local_offset, /* where in the global dataset to write */
      &vtk_hdf_htg   /* pointer to vtkHDFHyperTreeGrid object */
    );

#else
    hsize_t local_size[] = { (hsize_t)vtk_hdf_htg_data->descriptors_size };

    vtk_HDF_hypertreegrid_write_dataset(
      "Descriptors",                     /* dataset_name */
      vtk_hdf_htg_data->descriptors,     /* pointer to data */
      H5T_STD_U8LE,                      /* datatype */
      vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, /* group id */
      1,                                 /* rank */

      local_size,  /*dimensions of dataset */
      &vtk_hdf_htg /* pointer to vtkHDFHyperTreeGrid object */
    );
#endif
  }

  /*
   * Dataset: /VTKHDF/DescriptorsSize
   * Datatype: H5T_NATIVE_INT64 / int64_t
   * Dimension: {1*npe()}
   *
   * This is the number of descriptors (the individual bits, rather than the
   * number of bytes they pack into) of the tree on each process. Since we have
   * only one tree in basilisk, the dimensions of the dataset is always
   * `{1*npe()}`.
   */
  {
#if MPI_SINGLE_FILE
    hsize_t local_size[] = { 1 };
    hsize_t global_size[] = { (hsize_t)npe() };
    hsize_t local_offset[] = { (hsize_t)pid() };

    vtk_HDF_hypertreegrid_collective_write_dataset(
      "DescriptorsSize",                 /* dataset_name */
      &vtk_hdf_htg_data->n_descriptors,  /* pointer to data */
      H5T_NATIVE_INT64,                  /* datatype */
      vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, /* group id */
      1,                                 /* rank */

      global_size,  /* dimensions of dataset */
      local_size,   /* dimensions of local size */
      local_offset, /* where in the global dataset to write */
      &vtk_hdf_htg  /* pointer to vtkHDFHyperTreeGrid object */
    );
#else
    hsize_t local_size[] = { 1 };

    vtk_HDF_hypertreegrid_write_dataset(
      "DescriptorsSize",                 /* dataset_name */
      &vtk_hdf_htg_data->n_descriptors,  /* pointer to data */
      H5T_NATIVE_INT64,                  /* datatype */
      vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, /* group id */
      1,                                 /* rank */
      local_size,                        /*dimensions of dataset */
      &vtk_hdf_htg /* pointer to vtkHDFHyperTreeGrid object */
    );
#endif
  }

  /*
   * Dataset: /VTKHDF/Mask
   * Datatype: H5T_STD_U8LE / UInt8_t
   * Dimension: {global_size}
   */
  if (vtk_hdf_htg_data->has_mask) {
#if MPI_SINGLE_FILE

    hsize_t local_size = (hsize_t)vtk_hdf_htg_data->mask_size;
    hsize_t global_size = local_size;
    MPI_Allreduce(&local_size,
                  &global_size,
                  1,
                  MPI_UNSIGNED_LONG_LONG,
                  MPI_SUM,
                  MPI_COMM_WORLD);

    hsize_t local_offset = 0;
    MPI_Exscan(&local_size,
               &local_offset,
               1,
               MPI_UNSIGNED_LONG_LONG,
               MPI_SUM,
               MPI_COMM_WORLD);
    if (pid() == 0)
      local_offset = 0;

    // hsize_t chunk_size = (hsize_t)(global_size / npe());

    vtk_HDF_hypertreegrid_collective_write_dataset(
      "Mask",                            /* dataset_name */
      vtk_hdf_htg_data->mask,            /* pointer to data */
      H5T_STD_U8LE,                      /* datatype */
      vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, /* group id */
      1,                                 /* rank */
      &global_size,                      /* dimensions of dataset */
      &local_size,                       /* dimensions of local size */
      &local_offset, /* where in the global dataset to write */
      &vtk_hdf_htg   /* pointer to vtkHDFHyperTreeGrid object */
    );
#else
    hsize_t local_size[] = { (hsize_t)vtk_hdf_htg_data->mask_size };

    vtk_HDF_hypertreegrid_write_dataset(
      "Mask",                            /* dataset_name */
      vtk_hdf_htg_data->mask,            /* pointer to data */
      H5T_STD_U8LE,                      /* datatype */
      vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, /* group id */
      1,                                 /* rank */
      local_size,                        /* dimensions of dataset */
      &vtk_hdf_htg /* pointer to vtkHDFHyperTreeGrid object */
    );
#endif
  }

  /*
   * Dataset: /VTKHDF/NumberOfCells
   * Datatype: H5T_NATIVE_INT64 / int64_t
   * Dimension: {1*npe()}
   *
   * Here we write the number of cells in the tree.
   */
  {
#if MPI_SINGLE_FILE
    hsize_t local_size[] = { 1 };
    hsize_t global_size[] = { (hsize_t)npe() };
    hsize_t local_offset[] = { (hsize_t)pid() };

    vtk_HDF_hypertreegrid_collective_write_dataset(
      "NumberOfCells",                    /* dataset_name */
      &vtk_hdf_htg_data->number_of_cells, /* pointer to data */
      H5T_NATIVE_INT64,                   /* datatype */
      vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id,  /* group id */
      1,                                  /* rank */
      global_size,                        /* dimensions of dataset */
      local_size,                         /* dimensions of local size */
      local_offset, /* where in the global dataset to write */
      &vtk_hdf_htg  /* pointer to vtkHDFHyperTreeGrid object */
    );
#else
    hsize_t local_size[] = { 1 };

    vtk_HDF_hypertreegrid_write_dataset(
      "NumberOfCells",                    /* dataset_name */
      &vtk_hdf_htg_data->number_of_cells, /* pointer to data */
      H5T_NATIVE_INT64,                   /* datatype */
      vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id,  /* group id */
      1,                                  /* rank */
      local_size,                         /*dimensions of dataset */
      &vtk_hdf_htg /* pointer to vtkHDFHyperTreeGrid object */
    );
#endif
  }

  /*
   * Dataset: /VTKHDF/NumberOfCellsPerTreeDepth
   * Datatype: H5T_NATIVE_INT64 / int64_t
   * Dimension: {global_size}
   *
   * Here we write the number of cells at each level/depth in the tree. For
   * example the tree
   *
   * [Level 0]                      o
   *                               / \
   *                              /   \
   * [Level 1]                   o     o
   *                            / \
   *                           /   \
   * [Level 2]                o     o
   *
   * would write the array `int64_t depth_per_tree = {1 2 2}`.
   */
  {
#if MPI_SINGLE_FILE
    hsize_t local_size = (hsize_t)vtk_hdf_htg_data->depth_per_tree;
    hsize_t global_size = local_size;
    MPI_Allreduce(&local_size,
                  &global_size,
                  1,
                  MPI_UNSIGNED_LONG_LONG,
                  MPI_SUM,
                  MPI_COMM_WORLD);

    hsize_t local_offset = 0;
    MPI_Exscan(&local_size,
               &local_offset,
               1,
               MPI_UNSIGNED_LONG_LONG,
               MPI_SUM,
               MPI_COMM_WORLD);
    if (pid() == 0)
      local_offset = 0;

    vtk_HDF_hypertreegrid_collective_write_dataset(
      "NumberOfCellsPerTreeDepth",                      /* dataset_name */
      vtk_hdf_htg_data->number_of_cells_per_tree_depth, /* pointer to data */
      H5T_NATIVE_INT64,                                 /* datatype */
      vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id,                /* group id */
      1,                                                /*rank*/
      &global_size,  /* dimensions of dataset */
      &local_size,   /* dimensions of local size */
      &local_offset, /* where in the global dataset to write */
      &vtk_hdf_htg   /* pointer to vtkHDFHyperTreeGrid object */
    );
#else
    hsize_t local_size[] = { (hsize_t)vtk_hdf_htg_data->depth_per_tree };

    vtk_HDF_hypertreegrid_write_dataset(
      "NumberOfCellsPerTreeDepth",                      /* dataset_name */
      vtk_hdf_htg_data->number_of_cells_per_tree_depth, /* pointer to data */
      H5T_NATIVE_INT64,                                 /* datatype */
      vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id,                /* group id */
      1,                                                /* rank */
      local_size,  /* dimensions of local size */
      &vtk_hdf_htg /* pointer to vtkHDFHyperTreeGrid object */
    );
#endif
  }

  /*
   * Dataset: /VTKHDF/XCoordinates
   * Datatype: H5T_IEEE_F64LE / double
   * Dimension: {2 * npe()}
   *
   * In the MPI_SINGLE_FILE case we must (redudantly) write this for each
   * proccess.
   */
  {
#if MPI_SINGLE_FILE
    hsize_t local_size[] = { (hsize_t)vtk_hdf_htg_data->n_x };
    hsize_t global_size[] = { (hsize_t)npe() * local_size[0] };
    hsize_t local_offset[] = { (hsize_t)pid() * local_size[0] };
    vtk_HDF_hypertreegrid_collective_write_dataset(
      "XCoordinates",                    /* dataset_name */
      vtk_hdf_htg_data->x,               /* pointer to data */
      H5T_IEEE_F64LE,                    /* datatype */
      vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, /* group id */
      1,                                 /* rank */
      global_size,                       /* dimensions of dataset */
      local_size,                        /* dimensions of local size */
      local_offset, /* where in the global dataset to write */
      &vtk_hdf_htg  /* pointer to vtkHDFHyperTreeGrid object */
    );
#else
    hsize_t local_size[] = { (hsize_t)vtk_hdf_htg_data->n_x };
    vtk_HDF_hypertreegrid_write_dataset(
      "XCoordinates",                    /* dataset_name */
      vtk_hdf_htg_data->x,               /* pointer to data */
      H5T_IEEE_F64LE,                    /* datatype */
      vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, /* group id */
      1,                                 /* rank */
      local_size,                        /* dimensions of local size */
      &vtk_hdf_htg /* pointer to vtkHDFHyperTreeGrid object */
    );
#endif
  }

  /*
   * Dataset: /VTKHDF/YCoordinates
   * Datatype: H5T_IEEE_F64LE / double
   * Dimension: {(2 or 1) * npe()}
   *
   * In the MPI_SINGLE_FILE case we must (redudantly) write this for each
   * proccess.
   */
  {
#if MPI_SINGLE_FILE
    hsize_t local_size[] = { (hsize_t)vtk_hdf_htg_data->n_y };
    hsize_t global_size[] = { (hsize_t)npe() * local_size[0] };
    hsize_t local_offset[] = { (hsize_t)pid() * local_size[0] };
    vtk_HDF_hypertreegrid_collective_write_dataset(
      "YCoordinates",                    /* dataset_name */
      vtk_hdf_htg_data->y,               /* pointer to data */
      H5T_IEEE_F64LE,                    /* datatype */
      vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, /* group id */
      1,                                 /* rank */
      global_size,                       /* dimensions of dataset */
      local_size,                        /* dimensions of local size */
      local_offset, /* where in the global dataset to write */
      &vtk_hdf_htg  /* pointer to vtkHDFHyperTreeGrid object */
    );
#else
    hsize_t local_size[] = { (hsize_t)vtk_hdf_htg_data->n_y };
    vtk_HDF_hypertreegrid_write_dataset(
      "YCoordinates",                    /* dataset_name */
      vtk_hdf_htg_data->y,               /* pointer to data */
      H5T_IEEE_F64LE,                    /* datatype */
      vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, /* group id */
      1,                                 /* rank */
      local_size,                        /* dimensions of local size */
      &vtk_hdf_htg /* pointer to vtkHDFHyperTreeGrid object */
    );
#endif
  }

  /*
   * Dataset: /VTKHDF/ZCoordinates
   * Datatype: H5T_IEEE_F64LE / double
   * Dimension: {(2 or 1) * npe()}
   *
   * In the MPI_SINGLE_FILE case we must (redudantly) write this for each
   * proccess.
   */
  {
#if MPI_SINGLE_FILE
    hsize_t local_size[] = { (hsize_t)vtk_hdf_htg_data->n_z };
    hsize_t global_size[] = { (hsize_t)npe() * local_size[0] };
    hsize_t local_offset[] = { (hsize_t)pid() * local_size[0] };
    vtk_HDF_hypertreegrid_collective_write_dataset(
      "ZCoordinates",                    /* dataset_name */
      vtk_hdf_htg_data->z,               /* pointer to data */
      H5T_IEEE_F64LE,                    /* datatype */
      vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, /* group id */
      1,                                 /* rank */
      global_size,                       /* dimensions of dataset */
      local_size,                        /* dimensions of local size */
      local_offset, /* where in the global dataset to write */
      &vtk_hdf_htg  /* pointer to vtkHDFHyperTreeGrid object */
    );
#else
    hsize_t local_size[] = { (hsize_t)vtk_hdf_htg_data->n_z };
    vtk_HDF_hypertreegrid_write_dataset(
      "ZCoordinates",                    /* dataset_name */
      vtk_hdf_htg_data->z,               /* pointer to data */
      H5T_IEEE_F64LE,                    /* datatype */
      vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, /* group id */
      1,                                 /* rank */
      local_size,                        /* dimensions of local size */
      &vtk_hdf_htg /* pointer to vtkHDFHyperTreeGrid object */
    );
#endif
  }

  /*
   * Dataset: /VTKHDF/NumberOfDepths
   * Datatype: H5T_NATIVE_INT64 / int64_t
   * Dimension: {1 * npe()}
   *
   * The actual maximum depth per tree.
   */
  {
#if MPI_SINGLE_FILE
    hsize_t local_size[] = { 1 };
    hsize_t global_size[] = { (hsize_t)npe() };
    hsize_t local_offset[] = { (hsize_t)pid() };
    vtk_HDF_hypertreegrid_collective_write_dataset(
      "NumberOfDepths",                  /* dataset_name */
      &vtk_hdf_htg_data->depth_per_tree, /* pointer to data */
      H5T_NATIVE_INT64,                  /* datatype */
      vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, /* group id */
      1,                                 /* rank */
      global_size,                       /* dimensions of dataset */
      local_size,                        /* dimensions of local size */
      local_offset, /* where in the global dataset to write */
      &vtk_hdf_htg  /* pointer to vtkHDFHyperTreeGrid object */
    );
#else
    hsize_t local_size[] = { 1 };
    vtk_HDF_hypertreegrid_write_dataset(
      "NumberOfDepths",                  /* dataset_name */
      &vtk_hdf_htg_data->depth_per_tree, /* pointer to data */
      H5T_NATIVE_INT64,                  /* datatype */
      vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, /* group id */
      1,                                 /* rank */
      local_size,                        /* dimensions of local size */
      &vtk_hdf_htg /* pointer to vtkHDFHyperTreeGrid object */
    );
#endif
  }

  /*
   * Dataset: /VTKHDF/NumberOfTrees
   * Datatype: H5T_NATIVE_INT64 / int64_t
   * Dimension: {1 * npe()}
   *
   * The number of trees. Since basilisk can only ever have one tree, we create
   * a dataset with just 1.
   */
  {
#if MPI_SINGLE_FILE
    hsize_t local_size[] = { 1 };
    hsize_t global_size[] = { (hsize_t)npe() };
    hsize_t local_offset[] = { (hsize_t)pid() };
    vtk_HDF_hypertreegrid_collective_write_dataset(
      "NumberOfTrees",                    /* dataset_name */
      &vtk_hdf_htg_data->number_of_trees, /* pointer to data */
      H5T_NATIVE_INT64,                   /* datatype */
      vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id,  /* group id */
      1,                                  /* rank */
      global_size,                        /* dimensions of dataset */
      local_size,                         /* dimensions of local size */
      local_offset, /* where in the global dataset to write */
      &vtk_hdf_htg  /* pointer to vtkHDFHyperTreeGrid object */
    );
#else
    hsize_t local_size[] = { 1 };
    vtk_HDF_hypertreegrid_write_dataset(
      "NumberOfTrees",                    /* dataset_name */
      &vtk_hdf_htg_data->number_of_trees, /* pointer to data */
      H5T_NATIVE_INT64,                   /* datatype */
      vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id,  /* group id */
      1,                                  /* rank */
      local_size,                         /* dimensions of local size */
      &vtk_hdf_htg /* pointer to vtkHDFHyperTreeGrid object */
    );
#endif
  }

  /*
   * Dataset: /VTKHDF/TreeIds
   * Datatype: H5T_NATIVE_INT64 / int64_t
   * Dimension: {1 * npe()}
   *
   * The TreeIds of each tree. Again, since basilisk has only one tree, the
   * TreeId is always zero. However, the ordering in the case of a
   * forest-of-trees has implications in informing VTK which tree belongs to
   * which cell of the regular grid described in X/Y/Z coordinates.
   */
  {
#if MPI_SINGLE_FILE
    hsize_t local_size[] = { 1 };
    hsize_t global_size[] = { (hsize_t)npe() };
    hsize_t local_offset[] = { (hsize_t)pid() };
    vtk_HDF_hypertreegrid_collective_write_dataset(
      "TreeIds",                         /* dataset_name */
      &vtk_hdf_htg_data->tree_ids,       /* pointer to data */
      H5T_NATIVE_INT64,                  /* datatype */
      vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, /* group id */
      1,                                 /* rank */
      global_size,                       /* dimensions of dataset */
      local_size,                        /* dimensions of local size */
      local_offset, /* where in the global dataset to write */
      &vtk_hdf_htg  /* pointer to vtkHDFHyperTreeGrid object */
    );
#else
    hsize_t local_size[] = { 1 };
    vtk_HDF_hypertreegrid_write_dataset(
      "TreeIds",                         /* dataset_name */
      &vtk_hdf_htg_data->tree_ids,       /* pointer to data */
      H5T_NATIVE_INT64,                  /* datatype */
      vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, /* group id */
      1,                                 /* rank */
      local_size,                        /* dimensions of local size */
      &vtk_hdf_htg /* pointer to vtkHDFHyperTreeGrid object */
    );
#endif
  }

  /*
   * Dataset: /VTKHDF/CellData/*
   * Datatype: H5T_IEEE_F32LE / float
   * Dimension: {global_size} (Scalars) or {global_size,dimension} (Vectors)
   *
   * In the remaining part of the function, we write both scalar fields and
   * vector fields defined on the tree. These must be written in breadth-first
   * search order--the same order as the descriptors.
   *
   * In the MPI_SINGLE_FILE case, we write a PHyperTreeGrid, which neccesitates
   * a short MPI exchange to determine offsets. Since the fields are all
   * described by the same tree, we only need to do this once for all of the
   * vector and scalar fields.
   *
   * In addition, we will use this information to choose a chunk size. This
   * chunk size (measured in the number of elements) should be sized such that
   * we minimize the instances in which ranks/processes write into the same
   * chunk. Thus we set the chunk size to be the average of the local number of
   * cells.
   */
  {
#if MPI_SINGLE_FILE
    hsize_t scalar_local_size = (hsize_t)vtk_hdf_htg_data->number_of_cells;
    hsize_t scalar_global_size = scalar_local_size;
    MPI_Allreduce(&scalar_local_size,
                  &scalar_global_size,
                  1,
                  MPI_UNSIGNED_LONG_LONG,
                  MPI_SUM,
                  MPI_COMM_WORLD);
    hsize_t scalar_local_offset = 0;
    MPI_Exscan(&scalar_local_size,
               &scalar_local_offset,
               1,
               MPI_UNSIGNED_LONG_LONG,
               MPI_SUM,
               MPI_COMM_WORLD);
    if (pid() == 0) {
      scalar_local_offset = 0;
    }

    // hsize_t scalar_target_chunk_size = (1 << 20) / (hsize_t) sizeof(float);
    // // 1 MB hsize_t scalar_chunk_size = scalar_global_size <
    // scalar_target_chunk_size ? scalar_global_size : scalar_target_chunk_size;

    hsize_t scalar_chunk_size = scalar_global_size / (hsize_t)npe();
#else
    hsize_t scalar_local_size = (hsize_t)vtk_hdf_htg_data->number_of_cells;
    hsize_t scalar_global_size = scalar_local_size;
    hsize_t scalar_chunk_size = scalar_global_size;
#endif

    /*
     * Scalars
     */
    for (scalar s in scalar_list) {

      // Value for the Scalar dataset
      float* s_data = malloc(scalar_local_size * sizeof(float));

      // Malloc error handling
      if (!s_data) {
        perror("malloc(s_data)");
        exit(1);
      }

      // Copy the data from the tree into s_data
      size_t si = 0;
      foreach_cell_BFS()
      {
        bool write = false;
        if (is_leaf(cell)) {
          if (is_local(cell)) {
            write = true;
          } else {
            foreach_neighbor(1)
            {
              if (is_local(cell))
                write = true;
            }
          }
        }

        s_data[si++] = write ? (float)val(s) : 0.;
      }

#if COMPRESSION && MPI_SINGLE_FILE
      hsize_t dims[] = { scalar_global_size };
      hsize_t max_dims[] = { scalar_global_size };
      hsize_t chunk_dims[] = { scalar_chunk_size };
      hsize_t local_size[] = { scalar_local_size };
      hsize_t local_offset[] = { scalar_local_offset };

      vtk_HDF_hypertreegrid_collective_write_compressed_dataset(
        s.name,                      /* dataset_name */
        s_data,                      /* data */
        H5T_IEEE_F32LE,              /* dtype_id */
        vtk_hdf_htg.grp_celldata_id, /* group_id */
        1,                           /* rank */
        dims,                        /* dims[] */
        max_dims,                    /* max_dims[] */
        chunk_dims,                  /* chunk_dims[] */
        local_size,                  /* local_size[] */
        local_offset,                /* local_offset[] */
        COMPRESSION_LEVEL,           /* compression_level */
        &vtk_hdf_htg                 /* vtkHDFHyperTreeGrid */
      );
#elif COMPRESSION && !MPI_SINGLE_FILE
      hsize_t dims[] = { scalar_local_size };
      hsize_t max_dims[] = { scalar_local_size };
      hsize_t chunk_dims[] = { scalar_chunk_size };
      vtk_HDF_hypertreegrid_write_compressed_dataset(
        s.name,                      /* dataset_name */
        s_data,                      /* data */
        H5T_IEEE_F32LE,              /* dtype_id */
        vtk_hdf_htg.grp_celldata_id, /* group_id */
        1,                           /* rank */
        dims,                        /* dims[] */
        max_dims,                    /* max_dims[] */
        chunk_dims,                  /* chunk_dims[] */
        COMPRESSION_LEVEL,           /* compression_level */
        &vtk_hdf_htg                 /* vtkHDFHyperTreeGrid */
      );
#elif !COMPRESSION && MPI_SINGLE_FILE
      hsize_t dims[] = { scalar_global_size };
      hsize_t max_dims[] = { scalar_global_size };
      hsize_t local_size[] = { scalar_local_size };
      hsize_t local_offset[] = { scalar_local_offset };
      vtk_HDF_hypertreegrid_collective_write_dataset(
        s.name,                      /* dataset_name */
        s_data,                      /* data */
        H5T_IEEE_F32LE,              /* dtype_id */
        vtk_hdf_htg.grp_celldata_id, /* group_id */
        1,                           /* rank */
        dims,                        /* dims[] */
        max_dims,                    /* max_dims[] */
        local_size,                  /* local_size[] */
        local_offset,                /* local_offset[] */
        &vtk_hdf_htg                 /* vtkHDFHyperTreeGrid */
      );
#else // !COMPRESSION && !MPI_SINGLE_FILE
      hsize_t dims[] = { scalar_local_size };
      hsize_t max_dims[] = { scalar_local_size };
      vtk_HDF_hypertreegrid_collective_write_dataset(
        s.name,                      /* dataset_name */
        s_data,                      /* data */
        H5T_IEEE_F32LE,              /* dtype_id */
        vtk_hdf_htg.grp_celldata_id, /* group_id */
        1,                           /* rank */
        dims,                        /* dims[] */
        max_dims,                    /* max_dims[] */
        &vtk_hdf_htg                 /* vtkHDFHyperTreeGrid */
      );
#endif
      // Free the copied data
      free(s_data);

    } /* end of “for (scalar s in scalar_list)” */

    /*
     * Vectors
     */
    for (vector v in vector_list) {

      // Obtain the name of the vector
      char* vector_name;
      size_t trunc_len = (size_t)(strlen(v.x.name) - 2);
      vector_name = malloc((trunc_len + 1) * sizeof(char));
      strncpy(vector_name, v.x.name, trunc_len);
      vector_name[trunc_len] = '\0';

      // Create the BFS-ordered array
      float* v_data = malloc(scalar_local_size * dimension * sizeof(float));

      // Malloc error handling
      if (!v_data) {
        perror("malloc(v_data)");
        exit(1);
      }

      int vi = 0;
#if dimension == 1
      foreach_cell_BFS()
      {
        bool write = false;
        if (is_leaf(cell)) {
          if (is_local(cell)) {
            write = true;
          } else {
            foreach_neighbor(1)
            {
              if (is_local(cell))
                write = true;
            }
          }
        }
        v_data[vi * dimension + 0] = write ? (float)val(v.x) : 0.;
        vi++;
      }
#elif dimension == 2
      foreach_cell_BFS()
      {
        bool write = false;
        if (is_leaf(cell)) {
          if (is_local(cell)) {
            write = true;
          } else {
            foreach_neighbor(1)
            {
              if (is_local(cell))
                write = true;
            }
          }
        }
        v_data[vi * dimension + 0] = write ? (float)val(v.x) : 0.;
        v_data[vi * dimension + 1] = write ? (float)val(v.y) : 0.;
        vi++;
      }
#else // dimension == 3
      foreach_cell_BFS()
      {
        bool write = false;
        if (is_leaf(cell)) {
          if (is_local(cell)) {
            write = true;
          } else {
            foreach_neighbor(1)
            {
              if (is_local(cell))
                write = true;
            }
          }
        }
        v_data[vi * dimension + 0] = write ? (float)val(v.x) : 0.;
        v_data[vi * dimension + 1] = write ? (float)val(v.y) : 0.;
        v_data[vi * dimension + 2] = write ? (float)val(v.z) : 0.;
        vi++;
      }
#endif

#if COMPRESSION && MPI_SINGLE_FILE
      hsize_t dims[2] = { scalar_global_size, dimension };
      hsize_t max_dims[2] = { scalar_global_size, dimension };
      hsize_t chunk_dims[2] = { scalar_chunk_size, dimension };
      hsize_t local_size[2] = { scalar_local_size, dimension };
      hsize_t local_offset[2] = { scalar_local_offset, 0 };
      vtk_HDF_hypertreegrid_collective_write_compressed_dataset(
        vector_name,                 /* dataset_name */
        v_data,                      /* data */
        H5T_IEEE_F32LE,              /* dtype_id */
        vtk_hdf_htg.grp_celldata_id, /* group_id */
        2,                           /* rank */
        dims,                        /* dims[] */
        max_dims,                    /* max_dims[] */
        chunk_dims,                  /* chunk_dims[] */
        local_size,                  /* local_size[] */
        local_offset,                /* local_offset[] */
        COMPRESSION_LEVEL,           /* compression_level */
        &vtk_hdf_htg                 /* vtkHDFHyperTreeGrid */
      );
#elif COMPRESSION && !MPI_SINGLE_FILE
      hsize_t dims[] = { scalar_global_size, dimension };
      hsize_t max_dims[] = { scalar_global_size, dimension };
      hsize_t chunk_dims[] = { scalar_chunk_size, dimension };
      vtk_HDF_hypertreegrid_write_compressed_dataset(
        vector_name,                 /* dataset_name */
        v_data,                      /* data */
        H5T_IEEE_F32LE,              /* dtype_id */
        vtk_hdf_htg.grp_celldata_id, /* group_id */
        2,                           /* rank */
        dims,                        /* dims[] */
        max_dims,                    /* max_dims[] */
        chunk_dims,                  /* chunk_dims[] */
        COMPRESSION_LEVEL,           /* compression_level */
        &vtk_hdf_htg                 /* vtkHDFHyperTreeGrid */
      );
#elif !COMPRESSION && MPI_SINGLE_FILE
      hsize_t dims[2] = { scalar_global_size, dimension };
      hsize_t max_dims[2] = { scalar_global_size, dimension };
      hsize_t chunk_dims[2] = { scalar_chunk_size, dimension };
      hsize_t local_size[2] = { scalar_local_size, dimension };
      hsize_t local_offset[2] = { scalar_local_offset, 0 };
      vtk_HDF_hypertreegrid_collective_write_dataset(
        vector_name,                 /* dataset_name */
        v_data,                      /* data */
        H5T_IEEE_F32LE,              /* dtype_id */
        vtk_hdf_htg.grp_celldata_id, /* group_id */
        2,                           /* rank */
        dims,                        /* dims[] */
        max_dims,                    /* max_dims[] */
        local_size,                  /* local_size[] */
        local_offset,                /* local_offset[] */
        &vtk_hdf_htg                 /* vtkHDFHyperTreeGrid */
      );
#else // !COMPRESSION && !MPI_SINGLE_FILE
      hsize_t dims[] = { scalar_global_size, dimension };
      hsize_t max_dims[] = { scalar_global_size, dimension };
      vtk_HDF_hypertreegrid_collective_write_dataset(
        vector_name,                 /* dataset_name */
        v_data,                      /* data */
        H5T_IEEE_F32LE,              /* dtype_id */
        vtk_hdf_htg.grp_celldata_id, /* group_id */
        2,                           /* rank */
        dims,                        /* dims[] */
        max_dims,                    /* max_dims[] */
        &vtk_hdf_htg                 /* vtkHDFHyperTreeGrid */
      );
#endif
      // Free the copied data
      free(v_data);
    } /* end of “for (vector v in vector_list)” */
  }

  vtk_hdf_hypertreegrid_data_free(vtk_hdf_htg_data);

  return vtk_hdf_htg;
}
