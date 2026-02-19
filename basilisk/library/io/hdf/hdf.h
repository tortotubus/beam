#pragma once

#include <hdf5.h>
#include <mpi.h>

/**
 * @brief Parent class with many helper functions for writing VTKHDF files
 */
typedef struct {
  hid_t file_id; /**< File id */
  hid_t fapl;    /**< File access property list  */

  hid_t attr_id;       /**< Attribute id */
  hid_t attr_space_id; /**< Attribute dataspace id */
  hid_t attr_dtype_id; /**< Attribute datatype id */

  hid_t dset_id;       /**< Dataset id */
  hid_t dset_space_id; /**< Dataset dataspace id */
  hid_t dset_dtype_id; /**< Dataset datatype id */

  hid_t dcpl_id; /**< Dataset property-list id */

  hid_t file_space; /**< MPI-IO file space */
  hid_t mem_space;  /**< MPI-IO memory space */
  hid_t xfer_plist; /**< MPI-IO  */
} HDF;

/**
 * @brief Close a HDF file
 *
 * @param hdf Pointer or reference to HDF object/struct
 *
 * @memberof HDF
 */
void HDF_close (HDF* hdf) {
  // if (hdf->grp_vtkhdf_id >= 0)
  //   H5Gclose (hdf->grp_vtkhdf_id);

  if (hdf->xfer_plist >= 0)
    H5Pclose (hdf->xfer_plist);

  if (hdf->mem_space >= 0)
    H5Sclose (hdf->mem_space);

  if (hdf->file_space >= 0)
    H5Sclose (hdf->file_space);

  if (hdf->fapl >= 0)
    H5Pclose (hdf->fapl);

  if (hdf->file_id >= 0)
    H5Fclose (hdf->file_id);
}

/**
 * @brief Function called on any error which immediately closes the file to
 * avoid corruptions
 *
 * @param hdf Pointer or reference to HDF object/struct
 *
 * @memberof HDF
 */
void HDF_error (HDF* hdf) {
  HDF_close (hdf);
}

/**
 * @brief Check for error of some operation that returns herr_t type
 *
 * @param hdf Pointer or reference to HDF object/struct
 * @param result Result of a HDF5 operation to check
 *
 * @memberof HDF
 */
void HDF_check_result (HDF* hdf, herr_t result) {
  if (result < 0) {
    HDF_error (hdf);
  }
}

/**
 * @brief Check for error of some operation that returns an object id
 *
 * @param hdf Pointer or reference to HDF object/struct
 * @param object_id Object id to check
 *
 * @memberof HDF
 */
void HDF_check_object (HDF* hdf, hid_t object_id) {
  if (object_id <= H5I_INVALID_HID) {
    HDF_error (hdf);
  }
}

/**
 * @brief Initialize the HDF struct by writing a new file
 *
 * @memberof HDF
 */
HDF HDF_init (const char* fname, bool overwrite) {
  HDF hdf = {
    .file_id = H5I_INVALID_HID, /**< File id */
    .fapl = H5I_INVALID_HID,    /**< File access property list  */
    // .grp_vtkhdf_id = H5I_INVALID_HID, /**< VTKHDF group id */
    .attr_id = H5I_INVALID_HID,       /**< Attribute id */
    .attr_space_id = H5I_INVALID_HID, /**< Attribute dataspace id */
    .attr_dtype_id = H5I_INVALID_HID, /**< Attribute datatype id */
    .dset_id = H5I_INVALID_HID,       /**< Dataset id */
    .dset_space_id = H5I_INVALID_HID, /**< Dataset dataspace id */
    .dset_dtype_id = H5I_INVALID_HID, /**< Dataset datatype id */
    .dcpl_id = H5I_INVALID_HID,       /**< Dataset property-list id */
    .file_space = H5I_INVALID_HID,    /**< MPI-IO file space */
    .mem_space = H5I_INVALID_HID,     /**< MPI-IO memory space */
    .xfer_plist = H5I_INVALID_HID,    /**< MPI-IO  */
  };

  // Create the actual file on disk
  hdf.file_id = H5Fcreate (fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  HDF_check_object (&hdf, hdf.file_id);

  // Create the VTKHDF group
  // hdf.grp_vtkhdf_id = H5Gcreate2 (
  //   hdf.file_id, "/VTKHDF", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  // HDF_check_object (&hdf, hdf.grp_vtkhdf_id);

  return hdf;
}

HDF HDF_init_MPIIO (const char* fname, bool overwrite) {
  HDF hdf = {
    .file_id = H5I_INVALID_HID, /**< File id */
    .fapl = H5I_INVALID_HID,    /**< File access property list  */
    // .grp_vtkhdf_id = H5I_INVALID_HID, /**< VTKHDF group id */
    .attr_id = H5I_INVALID_HID,       /**< Attribute id */
    .attr_space_id = H5I_INVALID_HID, /**< Attribute dataspace id */
    .attr_dtype_id = H5I_INVALID_HID, /**< Attribute datatype id */
    .dset_id = H5I_INVALID_HID,       /**< Dataset id */
    .dset_space_id = H5I_INVALID_HID, /**< Dataset dataspace id */
    .dset_dtype_id = H5I_INVALID_HID, /**< Dataset datatype id */
    .dcpl_id = H5I_INVALID_HID,       /**< Dataset property-list id */
    .file_space = H5I_INVALID_HID,    /**< MPI-IO file space */
    .mem_space = H5I_INVALID_HID,     /**< MPI-IO memory space */
    .xfer_plist = H5I_INVALID_HID,    /**< MPI-IO  */
  };

  // Create a file access property list
  hdf.fapl = H5Pcreate (H5P_FILE_ACCESS);
  if (hdf.fapl < 0)
    HDF_error (&hdf);

  // If we use MPI_SINGLE_FILE, we set all procs to have access and write into
  // the same file
  if (H5Pset_fapl_mpio (hdf.fapl, MPI_COMM_WORLD, MPI_INFO_NULL) < 0)
    HDF_error (&hdf);

  // Create the acutal file on disk
  hdf.file_id = H5Fcreate (fname, H5F_ACC_TRUNC, H5P_DEFAULT, hdf.fapl);
  if (hdf.file_id < 0)
    HDF_error (&hdf);

  // Create the VTKHDF group
  // hdf.grp_vtkhdf_id = H5Gcreate2 (
  //   hdf.file_id, "/VTKHDF", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  // if (hdf.grp_vtkhdf_id < 0)
  //   HDF_error (&hdf);

  return hdf;
}

/**
 * @brief Initialize the HDF struct by opening an existing file
 *
 * @memberof HDF
 */
HDF HDF_open (const char* fname) {
  HDF hdf = {
    .file_id = H5I_INVALID_HID, /**< File id */
    .fapl = H5I_INVALID_HID,    /**< File access property list  */
    // .grp_vtkhdf_id = H5I_INVALID_HID, /**< VTKHDF group id */
    .attr_id = H5I_INVALID_HID,       /**< Attribute id */
    .attr_space_id = H5I_INVALID_HID, /**< Attribute dataspace id */
    .attr_dtype_id = H5I_INVALID_HID, /**< Attribute datatype id */
    .dset_id = H5I_INVALID_HID,       /**< Dataset id */
    .dset_space_id = H5I_INVALID_HID, /**< Dataset dataspace id */
    .dset_dtype_id = H5I_INVALID_HID, /**< Dataset datatype id */
    .dcpl_id = H5I_INVALID_HID,       /**< Dataset property-list id */
    .file_space = H5I_INVALID_HID,    /**< MPI-IO file space */
    .mem_space = H5I_INVALID_HID,     /**< MPI-IO memory space */
    .xfer_plist = H5I_INVALID_HID,    /**< MPI-IO  */
  };

  // Open the actual file on disk
  hdf.file_id = H5Fopen (fname, H5F_ACC_RDWR, H5P_DEFAULT);
  HDF_check_object (&hdf, hdf.file_id);

  // Open the VTKHDF group
  // hdf.grp_vtkhdf_id = H5Gopen (hdf.file_id, "/VTKHDF", H5P_DEFAULT);
  // HDF_check_object (&hdf, hdf.grp_vtkhdf_id);

  return hdf;
}

/**
 * @brief Initialize the HDF struct by opening an existing file
 *
 * @memberof HDF
 */
HDF HDF_open_MPIIO (const char* fname) {
  HDF hdf = {
    .file_id = H5I_INVALID_HID, /**< File id */
    .fapl = H5I_INVALID_HID,    /**< File access property list  */
    // .grp_vtkhdf_id = H5I_INVALID_HID, /**< VTKHDF group id */
    .attr_id = H5I_INVALID_HID,       /**< Attribute id */
    .attr_space_id = H5I_INVALID_HID, /**< Attribute dataspace id */
    .attr_dtype_id = H5I_INVALID_HID, /**< Attribute datatype id */
    .dset_id = H5I_INVALID_HID,       /**< Dataset id */
    .dset_space_id = H5I_INVALID_HID, /**< Dataset dataspace id */
    .dset_dtype_id = H5I_INVALID_HID, /**< Dataset datatype id */
    .dcpl_id = H5I_INVALID_HID,       /**< Dataset property-list id */
    .file_space = H5I_INVALID_HID,    /**< MPI-IO file space */
    .mem_space = H5I_INVALID_HID,     /**< MPI-IO memory space */
    .xfer_plist = H5I_INVALID_HID,    /**< MPI-IO  */
  };

  herr_t result = 0;
  // Create a file access property list
  hdf.fapl = H5Pcreate (H5P_FILE_ACCESS);
  HDF_check_object (&hdf, hdf.fapl);

  // If we use MPI_SINGLE_FILE, we set all procs to have access and write into
  // the same file
  result = H5Pset_fapl_mpio (hdf.fapl, MPI_COMM_WORLD, MPI_INFO_NULL);
  HDF_check_result (&hdf, result);

  // Open the actual file on disk
  hdf.file_id = H5Fopen (fname, H5F_ACC_RDWR, hdf.fapl);
  HDF_check_object (&hdf, hdf.file_id);

  // Open the VTKHDF group
  // hdf.grp_vtkhdf_id = H5Gopen (hdf.file_id, "/VTKHDF", H5P_DEFAULT);
  // HDF_check_object (&hdf, hdf.grp_vtkhdf_id);

  return hdf;
}

/**
 * @brief Write new attribute
 *
 * @param attribute_name The name of the attrbiute
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the attribute in
 * @param dims The dimensions of the data
 * @param hdf Pointer or reference to HDF object/struct
 *
 * @memberof HDF
 */
void HDF_write_attribute (const char* attribute_name,
                          const void* data,
                          const hid_t dtype_id,
                          const hid_t group_id,
                          const hsize_t dims[],
                          HDF* hdf) {
  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the attribute space
  hdf->attr_space_id = H5Screate_simple (1, dims, dims);
  HDF_check_object (hdf, hdf->attr_space_id);

  // Set the datatype
  hdf->attr_dtype_id = H5Tcopy (dtype_id);
  HDF_check_object (hdf, hdf->attr_dtype_id);

  // Create the attribute
  hdf->attr_id = H5Acreate2 (group_id,
                             attribute_name,
                             hdf->attr_dtype_id,
                             hdf->attr_space_id,
                             H5P_DEFAULT,
                             H5P_DEFAULT);
  HDF_check_object (hdf, hdf->attr_id);

  // Write the attribute
  result = H5Awrite (hdf->attr_id, hdf->attr_dtype_id, data);
  HDF_check_object (hdf, result);

  // Close the open objects
  H5Aclose (hdf->attr_id);
  H5Tclose (hdf->attr_dtype_id);
  H5Sclose (hdf->attr_space_id);
}

/**
 * @brief Write new scalar attribute
 *
 * @param attribute_name The name of the attrbiute
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the attribute in
 * @param hdf Pointer or reference to HDF object/struct
 *
 * @memberof HDF
 */

void HDF_write_scalar_attribute (const char* attribute_name,
                                 const void* data,
                                 const hid_t dtype_id,
                                 const hid_t group_id,
                                 HDF* hdf) {
  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the attribute space
  hdf->attr_space_id = H5Screate (H5S_SCALAR);
  HDF_check_object (hdf, hdf->attr_space_id);

  // Set the datatype
  hdf->attr_dtype_id = H5Tcopy (dtype_id);
  HDF_check_object (hdf, hdf->attr_dtype_id);

  // Create the attribute
  hdf->attr_id = H5Acreate2 (group_id,
                             attribute_name,
                             hdf->attr_dtype_id,
                             hdf->attr_space_id,
                             H5P_DEFAULT,
                             H5P_DEFAULT);
  HDF_check_object (hdf, hdf->attr_id);

  // Write the attribute
  result = H5Awrite (hdf->attr_id, hdf->attr_dtype_id, data);
  HDF_check_object (hdf, result);

  // Close the open objects
  H5Aclose (hdf->attr_id);
  H5Tclose (hdf->attr_dtype_id);
  H5Sclose (hdf->attr_space_id);
}

/**
 * @brief Modify existing scalar attribute
 *
 * @param attribute_name The name of the attrbiute
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the attribute in
 * @param hdf Pointer or reference to HDF object/struct
 *
 * @memberof HDF
 */
void HDF_modify_scalar_attribute (const char* attribute_name,
                                  const void* data,
                                  const hid_t dtype_id,
                                  const hid_t group_id,
                                  HDF* hdf) {
  herr_t result = 0;

  // 1. Open existing attribute
  hdf->attr_id = H5Aopen (group_id, attribute_name, H5P_DEFAULT);
  HDF_check_object (hdf, hdf->attr_id);

  // 2. (Optional but nice) get its dataspace and verify it is scalar
  hdf->attr_space_id = H5Aget_space (hdf->attr_id);
  HDF_check_object (hdf, hdf->attr_space_id);

  int ndims = H5Sget_simple_extent_ndims (hdf->attr_space_id);
  HDF_check_result (hdf, ndims);
  if (ndims != 0) {
    perror ("HDF_modify_scalar_attribute: attribute is not scalar.\n");
    exit (1);
  }

  // 3. (Optional) get file datatype if you want to sanity-check or inspect it
  hdf->attr_dtype_id = H5Aget_type (hdf->attr_id);
  HDF_check_object (hdf, hdf->attr_dtype_id);
  // You *could* compare mem_dtype_id vs attr_dtype_id here if desired.

  // 4. Overwrite the attribute value in-place
  //    mem_dtype_id is the in-memory type (often H5T_NATIVE_*)
  result = H5Awrite (hdf->attr_id, dtype_id, data);
  HDF_check_result (hdf, result);

  // 5. Close only what we actually opened
  H5Tclose (hdf->attr_dtype_id); // from H5Aget_type
  H5Sclose (hdf->attr_space_id); // from H5Aget_space
  H5Aclose (hdf->attr_id);       // from H5Aopen
}

/**
 * @brief Read scalar attribute
 *
 * @param attribute_name The name of the attrbiute
 * @param data Pointer to the array of data
 * @param group_id The HDF5 group to write the attribute in
 * @param hdf Pointer or reference to HDF object/struct
 *
 * @memberof HDF
 */
void HDF_read_scalar_attribute (const char* attribute_name,
                                void* data,
                                const hid_t group_id,
                                HDF* hdf) {

  herr_t result;

  // Open the attribute
  hdf->attr_id = H5Aopen (group_id, attribute_name, H5P_DEFAULT);
  HDF_check_object (hdf, hdf->attr_id);

  // Get the datatype
  hdf->attr_dtype_id = H5Aget_type (hdf->attr_id);
  HDF_check_object (hdf, hdf->attr_dtype_id);

  // Get the dataspace
  hdf->attr_space_id = H5Aget_space (hdf->attr_id);
  HDF_check_object (hdf, hdf->attr_space_id);

  // Read the scalar value into data pointer
  result = H5Aread (hdf->attr_id, hdf->attr_dtype_id, data);
  HDF_check_result (hdf, result);

  H5Aclose (hdf->attr_id);
  H5Tclose (hdf->attr_dtype_id);
  H5Sclose (hdf->attr_space_id);
}

/**
 * @brief Write new attribute
 *
 * @param type_name The name of the attrbiute
 * @param group_id The HDF5 group to write the attribute in
 * @param hdf Pointer or reference to HDF object/struct
 *
 * @memberof HDF
 */
void HDF_write_type_attribute (const char* type_name,
                               const hid_t group_id,
                               HDF* hdf) {
  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the attribute space
  hdf->attr_space_id = H5Screate (H5S_SCALAR);
  HDF_check_object (hdf, hdf->attr_space_id);

  // Set the datatype
  hdf->attr_dtype_id = H5Tcopy (H5T_C_S1);
  HDF_check_object (hdf, hdf->attr_dtype_id);

  // Set the size
  result = H5Tset_size (hdf->attr_dtype_id, strlen (type_name));
  HDF_check_result (hdf, result);

  // Pad the string
  result = H5Tset_strpad (hdf->attr_dtype_id, H5T_STR_NULLTERM);
  HDF_check_result (hdf, result);

  // Set the string value
  result = H5Tset_cset (hdf->attr_dtype_id, H5T_CSET_ASCII);
  HDF_check_result (hdf, result);

  // Create the attribute
  hdf->attr_id = H5Acreate2 (group_id,
                             "Type",
                             hdf->attr_dtype_id,
                             hdf->attr_space_id,
                             H5P_DEFAULT,
                             H5P_DEFAULT);
  HDF_check_object (hdf, hdf->attr_id);

  // Write the attribute
  result = H5Awrite (hdf->attr_id, hdf->attr_dtype_id, type_name);
  HDF_check_object (hdf, result);

  // Close the open objects
  H5Aclose (hdf->attr_id);
  H5Tclose (hdf->attr_dtype_id);
  H5Sclose (hdf->attr_space_id);
}

/**
 * @brief Read dataset
 *
 * @param dataset_name The name of the dataset
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the data in
 * @param rank The rank of the dataset (e.g. the size of the dims array)
 * @param dims Pointer to an array of hsize_t which will afterwards contain the
 * dims of the dataset
 * @param max_dims Pointer to an array of hsize_t which will afterwards contain
 * the maximum dimensions of the dataset
 * @param hdf Pointer or reference to HDF object/struct
 *
 * @memberof HDF
 */
void HDF_read_dataset (const char* dataset_name,
                       void** data,
                       const hid_t dtype_id,
                       const hid_t group_id,
                       const int rank,
                       hsize_t dims[],
                       hsize_t max_dims[],
                       HDF* hdf) {
  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Open the dataset
  hdf->dset_id = H5Dopen2 (group_id, dataset_name, H5P_DEFAULT);
  HDF_check_object (hdf, hdf->dset_id);

  // Dataspace
  hdf->dset_space_id = H5Dget_space (hdf->dset_id);
  HDF_check_object (hdf, hdf->dset_space_id);

  // Rank
  int file_rank = H5Sget_simple_extent_ndims (hdf->dset_space_id);

  if (file_rank != rank) {
    perror ("HDF_read_dataset : Rank mismatch.\n");
    exit (1);
  }

  // Dims / max_dims
  result = H5Sget_simple_extent_dims (hdf->dset_space_id, dims, max_dims);

  // Compute the number of total elements
  size_t n_elems = 1;
  for (int i = 0; i < rank; ++i) {
    n_elems *= (size_t) dims[i];
  }

  // Size of one element (from mem_type_id)
  size_t type_size = (size_t) H5Tget_size (dtype_id);

  // Allocate the buffer
  *data = malloc (n_elems * type_size);

  if (!*data) {
    perror ("malloc failed\n");
    exit (1);
  }

  result =
    H5Dread (hdf->dset_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, *data);
  HDF_check_result (hdf, result);

  // Close opened objects
  H5Dclose (hdf->dset_id);
  H5Sclose (hdf->dset_space_id);
}

/**
 * @brief Read dataset
 *
 * @param dataset_name The name of the dataset
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the data in
 * @param rank The rank of the dataset (e.g. the size of the dims array)
 * @param dims Pointer to an array of hsize_t which will afterwards contain the
 * dims of the dataset
 * @param max_dims Pointer to an array of hsize_t which will afterwards contain
 * the maximum dimensions of the dataset
 * @param hdf Pointer or reference to HDF object/struct
 *
 * @memberof HDF
 */
void HDF_read_dataset_dims (const char* dataset_name,
                            const hid_t dtype_id,
                            const hid_t group_id,
                            const int rank,
                            hsize_t dims[],
                            hsize_t max_dims[],
                            HDF* hdf) {
  // Store the result of HDF5 functions that do not return an identifier
  herr_t result;

  // Open the dataset
  hdf->dset_id = H5Dopen2 (group_id, dataset_name, H5P_DEFAULT);
  HDF_check_object (hdf, hdf->dset_id);

  // Dataspace
  hdf->dset_space_id = H5Dget_space (hdf->dset_id);
  HDF_check_object (hdf, hdf->dset_space_id);

  // Rank
  int file_rank = H5Sget_simple_extent_ndims (hdf->dset_space_id);

  if (file_rank != rank) {
    perror ("HDF_read_dataset_dims: Rank mismatch.\n");
    exit (1);
  }

  // Dims / max_dims
  result = H5Sget_simple_extent_dims (hdf->dset_space_id, dims, max_dims);
  HDF_check_result (hdf, result);

  // Close opened objects
  H5Dclose (hdf->dset_id);
  H5Sclose (hdf->dset_space_id);
}

/**
 * @brief Write new unchunked dataset
 *
 * @param dataset_name The name of the dataset
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the data in
 * @param rank The rank of the dataset (e.g. the size of the dims array)
 * @param dims The dimensions of the dataset
 * @param hdf Pointer or reference to HDF object/struct
 *
 * @memberof HDF
 */
void HDF_write_dataset (const char* dataset_name,
                        const void* data,
                        const hid_t dtype_id,
                        const hid_t group_id,
                        const int rank,
                        const hsize_t dims[],
                        HDF* hdf) {
  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the data space with dimensions and maximum dimensions
  hdf->dset_space_id = H5Screate_simple (rank, dims, dims);
  HDF_check_object (hdf, hdf->dset_space_id);

  // Create a property list for this dataset
  hdf->dcpl_id = H5Pcreate (H5P_DATASET_CREATE);
  HDF_check_object (hdf, hdf->dcpl_id);

  // Set the datatype
  hdf->dset_dtype_id = H5Tcopy (dtype_id);
  HDF_check_object (hdf, hdf->dset_dtype_id);

  // Create the dataset
  hdf->dset_id = H5Dcreate2 (group_id,
                             dataset_name,
                             hdf->dset_dtype_id,
                             hdf->dset_space_id,
                             H5P_DEFAULT,
                             hdf->dcpl_id,
                             H5P_DEFAULT);
  HDF_check_object (hdf, hdf->dset_id);

  result = H5Dwrite (
    hdf->dset_id, hdf->dset_dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  HDF_check_result (hdf, result);

  // Close opened objects
  H5Dclose (hdf->dset_id);
  H5Tclose (hdf->dset_dtype_id);
  H5Pclose (hdf->dcpl_id);
  H5Sclose (hdf->dset_space_id);
}

/**
 * @brief Write new chunked dataset
 *
 * @param dataset_name The name of the dataset
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the data in
 * @param rank The rank of the dataset (e.g. the size of the dims array)
 * @param dims The dimensions of the dataset
 * @param max_dims The maximum dimensions of the dataset
 * @param chunk_dims The chunk dimensions
 * @param hdf Pointer or reference to HDF object/struct
 */
void HDF_write_chunked_dataset (const char* dataset_name,
                                const void* data,
                                const hid_t dtype_id,
                                const hid_t group_id,
                                const int rank,
                                const hsize_t dims[],
                                const hsize_t max_dims[],
                                const hsize_t chunk_dims[],
                                HDF* hdf) {

  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the data space with dimensions and maximum dimensions
  hdf->dset_space_id = H5Screate_simple (rank, dims, max_dims);
  HDF_check_object (hdf, hdf->dset_space_id);

  // Create a property list for this dataset
  hdf->dcpl_id = H5Pcreate (H5P_DATASET_CREATE);
  HDF_check_object (hdf, hdf->dcpl_id);

  // Set the chunk size for the dataset
  result = H5Pset_chunk (hdf->dcpl_id, rank, chunk_dims);
  HDF_check_result (hdf, result);

  // Set the datatype
  hdf->dset_dtype_id = H5Tcopy (dtype_id);
  HDF_check_object (hdf, hdf->dset_dtype_id);

  // Create the dataset
  hdf->dset_id = H5Dcreate2 (group_id,
                             dataset_name,
                             hdf->dset_dtype_id,
                             hdf->dset_space_id,
                             H5P_DEFAULT,
                             hdf->dcpl_id,
                             H5P_DEFAULT);
  HDF_check_object (hdf, hdf->dset_id);

  result = H5Dwrite (
    hdf->dset_id, hdf->dset_dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  HDF_check_result (hdf, result);

  // Close opened objects
  H5Dclose (hdf->dset_id);
  H5Tclose (hdf->dset_dtype_id);
  H5Pclose (hdf->dcpl_id);
  H5Sclose (hdf->dset_space_id);
}

/**
 * @brief Append existing chunked dataset
 *
 * @param dataset_name The name of the dataset
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the data in
 * @param rank The rank of the dataset (e.g. the size of the dims array)
 * @param dims The dimensions of the dataset to append
 * @param hdf Pointer or reference to HDF object/struct
 */

void HDF_append_chunked_dataset (const char* dataset_name,
                                 const void* data,
                                 const hid_t dtype_id,
                                 const hid_t group_id,
                                 const int rank,
                                 const hsize_t dims[],
                                 HDF* hdf) {
  herr_t result = 0;

  // Open existing dataset
  hdf->dset_id = H5Dopen2 (group_id, dataset_name, H5P_DEFAULT);
  HDF_check_object (hdf, hdf->dset_id);

  // Get current dataspace and its dimensions
  hdf->dset_space_id = H5Dget_space (hdf->dset_id);
  HDF_check_object (hdf, hdf->dset_space_id);

  hsize_t cur_dims[H5S_MAX_RANK];
  hsize_t cur_max_dims[H5S_MAX_RANK];
  int file_rank =
    H5Sget_simple_extent_dims (hdf->dset_space_id, cur_dims, cur_max_dims);
  HDF_check_result (hdf, file_rank);

  // Check that the rank should match
  if (file_rank != rank) {
    perror ("Rank mismatch when appending chunked dataset.\n");
    exit (1);
  }

  // We assume we append along dimension 0, and that all other dimensions stay
  // the same.
  hsize_t new_dims[H5S_MAX_RANK];
  for (int i = 0; i < rank; ++i) {
    new_dims[i] = cur_dims[i];
  }
  new_dims[0] = cur_dims[0] + dims[0];

  // Ensure other dims are unchanged
  for (int i = 1; i < rank; ++i) {
    if (dims[i] != cur_dims[i]) {
      perror ("Append dims must match existing dims in non-unlimited axes.\n");
      exit (1);
    }
  }

  // Extend the dataset to the new size
  result = H5Dset_extent (hdf->dset_id, new_dims);
  HDF_check_result (hdf, result);

  // After extending, the old dataspace is invalid; close and re-open
  H5Sclose (hdf->dset_space_id);
  hdf->dset_space_id = H5Dget_space (hdf->dset_id);
  HDF_check_object (hdf, hdf->dset_space_id);

  // 4. Select hyperslab corresponding to the newly appended region
  hsize_t start[H5S_MAX_RANK];
  hsize_t count[H5S_MAX_RANK];

  start[0] = cur_dims[0]; // start right after old end
  count[0] = dims[0];     // append this many along dim 0

  for (int i = 1; i < rank; ++i) {
    start[i] = 0;
    count[i] = dims[i]; // whole extent in other dims
  }

  result = H5Sselect_hyperslab (hdf->dset_space_id,
                                H5S_SELECT_SET,
                                start,
                                NULL, // stride
                                count,
                                NULL); // block

  HDF_check_result (hdf, result);

  // Create memspace for the data block we’re appending
  hdf->mem_space = H5Screate_simple (rank, dims, NULL);
  HDF_check_object (hdf, hdf->mem_space);

  // Write the new chunk into the extended portion
  result = H5Dwrite (hdf->dset_id,
                     dtype_id,           // in-memory type
                     hdf->mem_space,     // mem dataspace selection
                     hdf->dset_space_id, // file dataspace selection
                     H5P_DEFAULT,        // xfer plist
                     data);
  HDF_check_result (hdf, result);

  // Close everything
  H5Sclose (hdf->mem_space);
  H5Sclose (hdf->dset_space_id);
  H5Dclose (hdf->dset_id);
}

/**
 * @brief Write new chunked and compressed dataset
 *
 * @param dataset_name The name of the dataset
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the data in
 * @param rank The rank of the dataset (e.g. the size of the dims array)
 * @param dims The dimensions of the dataset
 * @param max_dims The maximum dimensions of the dataset
 * @param chunk_dims The chunk dimensions
 * @param compression_level The level of compression
 * @param hdf Pointer or reference to HDF object/struct
 *
 * @memberof HDF
 */
void HDF_write_compressed_dataset (const char* dataset_name,
                                   const void* data,
                                   const hid_t dtype_id,
                                   const hid_t group_id,
                                   const int rank,
                                   const hsize_t dims[],
                                   const hsize_t max_dims[],
                                   const hsize_t chunk_dims[],
                                   unsigned int compression_level,
                                   HDF* hdf) {

  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the data space with dimensions and maximum dimensions
  hdf->dset_space_id = H5Screate_simple (rank, dims, max_dims);
  HDF_check_object (hdf, hdf->dset_space_id);

  // Create a property list for this dataset
  hdf->dcpl_id = H5Pcreate (H5P_DATASET_CREATE);
  HDF_check_object (hdf, hdf->dcpl_id);

  // Set the chunk size for the dataset
  result = H5Pset_chunk (hdf->dcpl_id, rank, chunk_dims);
  HDF_check_result (hdf, result);

  // Set compression on the dataset
  result = H5Pset_deflate (hdf->dcpl_id, compression_level);
  HDF_check_result (hdf, result);

  // Set the datatype
  hdf->dset_dtype_id = H5Tcopy (dtype_id);
  HDF_check_object (hdf, hdf->dset_dtype_id);

  // Create the dataset
  hdf->dset_id = H5Dcreate2 (group_id,
                             dataset_name,
                             hdf->dset_dtype_id,
                             hdf->dset_space_id,
                             H5P_DEFAULT,
                             hdf->dcpl_id,
                             H5P_DEFAULT);
  HDF_check_object (hdf, hdf->dset_id);

  // Write the dataset
  result = H5Dwrite (
    hdf->dset_id, hdf->dset_dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  HDF_check_result (hdf, result);

  // Close opened objects
  H5Dclose (hdf->dset_id);
  H5Tclose (hdf->dset_dtype_id);
  H5Pclose (hdf->dcpl_id);
  H5Sclose (hdf->dset_space_id);
}

/**
 * @brief Write a collectively written dataset with MPI-IO
 *
 * @param dataset_name The name of the dataset
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the data in
 * @param rank The rank of the dataset (e.g. the size of the dims array)
 * @param dims The dimensions of the dataset
 * @param local_size The size of the sub-array this process writes into
 * @param local_offset The position of the sub-array this process writes into
 * @param hdf Pointer or reference to HDF object/struct
 *
 * @memberof HDF
 */
void HDF_collective_write_dataset (const char* dataset_name,
                                   const void* data,
                                   const hid_t dtype_id,
                                   const hid_t group_id,
                                   const int rank,
                                   const hsize_t dims[],
                                   const hsize_t local_size[],
                                   const hsize_t local_offset[],
                                   HDF* hdf) {

  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the data space with dimensions and maximum dimensions
  hdf->dset_space_id = H5Screate_simple (rank, dims, dims);
  HDF_check_object (hdf, hdf->dset_space_id);

  // Set the datatype
  hdf->dset_dtype_id = H5Tcopy (dtype_id);
  HDF_check_object (hdf, hdf->dset_dtype_id);

  // Create the dataset
  hdf->dset_id = H5Dcreate2 (group_id,
                             dataset_name,
                             hdf->dset_dtype_id,
                             hdf->dset_space_id,
                             H5P_DEFAULT,
                             H5P_DEFAULT,
                             H5P_DEFAULT);
  HDF_check_object (hdf, hdf->dset_id);

  // Create a property list for MPI-IO transfer
  hdf->xfer_plist = H5Pcreate (H5P_DATASET_XFER);
  HDF_check_object (hdf, hdf->xfer_plist);

  // Set MPIO_COLLECTIVE writing
  result = H5Pset_dxpl_mpio (hdf->xfer_plist, H5FD_MPIO_COLLECTIVE);
  HDF_check_result (hdf, result);

  // Get the dataset space
  hdf->file_space = H5Dget_space (hdf->dset_id);
  HDF_check_object (hdf, hdf->file_space);

  // Select a hyperslab for our process
  result = H5Sselect_hyperslab (
    hdf->file_space, H5S_SELECT_SET, local_offset, NULL, local_size, NULL);
  HDF_check_result (hdf, result);

  // Create a memory dataspace for our process
  hdf->mem_space = H5Screate_simple (rank, local_size, NULL);
  HDF_check_object (hdf, hdf->mem_space);

  // Actually do the collective write to the file
  result =
    H5Dwrite (hdf->dset_id,       /* dataset handle */
              hdf->dset_dtype_id, /* H5T_IEEE_F64LE */
              hdf->mem_space,     /* memory dataspace [local_nx] */
              hdf->file_space,    /* file dataspace with hyperslab selected */
              hdf->xfer_plist,    /* collective MPI‐IO transfer property */
              data                /* pointer to local data */
    );
  HDF_check_result (hdf, result);

  // Close opened objects
  H5Pclose (hdf->xfer_plist);
  hdf->xfer_plist = H5I_INVALID_HID;
  H5Sclose (hdf->mem_space);
  hdf->mem_space = H5I_INVALID_HID;
  H5Sclose (hdf->file_space);
  hdf->file_space = H5I_INVALID_HID;
  H5Tclose (hdf->dset_dtype_id);
  H5Sclose (hdf->dset_space_id);
  H5Dclose (hdf->dset_id);
}

/**
 * @brief Write a chunked dataset collectively using MPI-IO
 *
 * @param dataset_name The name of the dataset
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the data in
 * @param dims The dimensions of the dataset
 * @param rank The rank of the dataset (e.g. the size of the dims array)
 * @param max_dims The maximum dimensions of the dataset
 * @param chunk_dims The chunk dimensions
 * @param local_size The size of the sub-array this process writes into
 * @param local_offset The position of the sub-array this process writes into
 * @param hdf Pointer or reference to HDF object/struct
 *
 * @memberof HDF
 */
void HDF_collective_write_chunked_dataset (const char* dataset_name,
                                           const void* data,
                                           const hid_t dtype_id,
                                           const hid_t group_id,
                                           const int rank,
                                           const hsize_t dims[],
                                           const hsize_t max_dims[],
                                           const hsize_t chunk_dims[],
                                           const hsize_t local_size[],
                                           const hsize_t local_offset[],
                                           HDF* hdf) {

  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the data space with dimensions and maximum dimensions
  hdf->dset_space_id = H5Screate_simple (rank, dims, max_dims);
  HDF_check_object (hdf, hdf->dset_space_id);

  // Create a property list for this dataset
  hdf->dcpl_id = H5Pcreate (H5P_DATASET_CREATE);
  HDF_check_object (hdf, hdf->dcpl_id);

  // Set the chunk size for the dataset
  result = H5Pset_chunk (hdf->dcpl_id, rank, chunk_dims);
  HDF_check_result (hdf, result);

  // Set the datatype
  hdf->dset_dtype_id = H5Tcopy (dtype_id);
  HDF_check_object (hdf, hdf->dset_dtype_id);

  // Create the dataset
  hdf->dset_id = H5Dcreate2 (group_id,
                             dataset_name,
                             hdf->dset_dtype_id,
                             hdf->dset_space_id,
                             H5P_DEFAULT,
                             hdf->dcpl_id,
                             H5P_DEFAULT);
  HDF_check_object (hdf, hdf->dset_id);

  // Create a property list for MPI-IO transfer
  hdf->xfer_plist = H5Pcreate (H5P_DATASET_XFER);
  HDF_check_object (hdf, hdf->xfer_plist);

  // Set MPIO_COLLECTIVE writing
  result = H5Pset_dxpl_mpio (hdf->xfer_plist, H5FD_MPIO_COLLECTIVE);
  HDF_check_result (hdf, result);

  // Get the dataset space
  hdf->file_space = H5Dget_space (hdf->dset_id);
  HDF_check_object (hdf, hdf->file_space);

  // Select a hyperslab for our process
  result = H5Sselect_hyperslab (
    hdf->file_space, H5S_SELECT_SET, local_offset, NULL, local_size, NULL);
  HDF_check_result (hdf, result);

  // Create a memory dataspace for our process
  hdf->mem_space = H5Screate_simple (rank, local_size, NULL);
  HDF_check_object (hdf, hdf->mem_space);

  // Actually do the collective write to the file
  result =
    H5Dwrite (hdf->dset_id,       /* dataset handle */
              hdf->dset_dtype_id, /* H5T_IEEE_F64LE */
              hdf->mem_space,     /* memory dataspace */
              hdf->file_space,    /* file dataspace with hyperslab selected */
              hdf->xfer_plist,    /* collective MPI‐IO transfer property */
              data                /* pointer to local data */
    );
  HDF_check_result (hdf, result);

  // Close opened objects
  H5Pclose (hdf->xfer_plist);
  hdf->xfer_plist = H5I_INVALID_HID;
  H5Sclose (hdf->mem_space);
  hdf->mem_space = H5I_INVALID_HID;
  H5Sclose (hdf->file_space);
  hdf->file_space = H5I_INVALID_HID;
  H5Tclose (hdf->dset_dtype_id);
  H5Pclose (hdf->dcpl_id);
  H5Sclose (hdf->dset_space_id);
  H5Dclose (hdf->dset_id);
}

/**
 * @brief Collectiviely append existing chunked dataset
 *
 * @param dataset_name The name of the dataset
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the data in
 * @param rank The rank of the dataset (e.g. the size of the dims array)
 * @param dims The dimensions of the dataset to append
 * @param local_size The size of the sub-array this process writes into
 * @param local_offset The position of the sub-array this process writes into
 * @param hdf Pointer or reference to HDF object/struct
 */
void HDF_collective_append_chunked_dataset (
  const char* dataset_name,
  const void* data,
  const hid_t dtype_id,
  const hid_t group_id,
  const int rank,
  const hsize_t dims[],         // global size of appended block
  const hsize_t local_size[],   // local tile within appended block
  const hsize_t local_offset[], // local offset within appended block
  HDF* hdf) {
  herr_t result = 0;

  // 1. Open existing dataset collectively
  hdf->dset_id = H5Dopen2 (group_id, dataset_name, H5P_DEFAULT);
  HDF_check_object (hdf, hdf->dset_id);

  // 2. Get current global dims
  hdf->dset_space_id = H5Dget_space (hdf->dset_id);
  HDF_check_object (hdf, hdf->dset_space_id);

  hsize_t cur_dims[H5S_MAX_RANK];
  hsize_t cur_max_dims[H5S_MAX_RANK];
  int file_rank =
    H5Sget_simple_extent_dims (hdf->dset_space_id, cur_dims, cur_max_dims);
  if (file_rank < 0)
    HDF_error (hdf);
  if (file_rank != rank) {
    fprintf (stderr,
             "HDF_collective_append_chunked_dataset: rank mismatch "
             "(file_rank=%d, expected=%d)\n",
             file_rank,
             rank);
    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  // 3. Compute new global dims: append along dim 0 only
  hsize_t new_dims[H5S_MAX_RANK];
  for (int i = 0; i < rank; ++i)
    new_dims[i] = cur_dims[i];

  new_dims[0] = cur_dims[0] + dims[0];

  // Sanity: other dims must match existing ones
  for (int i = 1; i < rank; ++i) {
    if (dims[i] != cur_dims[i]) {
      fprintf (stderr,
               "Append dims mismatch in dim %d: existing=%llu, appended=%llu\n",
               i,
               (unsigned long long) cur_dims[i],
               (unsigned long long) dims[i]);
      MPI_Abort (MPI_COMM_WORLD, 1);
    }
  }

  // 4. Extend dataset (collective)
  result = H5Dset_extent (hdf->dset_id, new_dims);
  HDF_check_result (hdf, result);

  // Old dataspace no longer valid
  H5Sclose (hdf->dset_space_id);
  hdf->dset_space_id = H5I_INVALID_HID;

  // 5. Get new file dataspace for selection
  hdf->file_space = H5Dget_space (hdf->dset_id);
  HDF_check_object (hdf, hdf->file_space);

  // 6. Select hyperslab corresponding to *this rank's tile* in appended region
  hsize_t start[H5S_MAX_RANK];
  hsize_t count[H5S_MAX_RANK];

  // We append along dim 0: new slab starts at cur_dims[0]
  // Contract: local_offset[0] is measured *within* the appended slab,
  // so we shift it by cur_dims[0] in the global file space.
  start[0] = cur_dims[0] + local_offset[0];
  count[0] = local_size[0];

  for (int i = 1; i < rank; ++i) {
    start[i] = local_offset[i];
    count[i] = local_size[i];
  }

  result = H5Sselect_hyperslab (hdf->file_space,
                                H5S_SELECT_SET,
                                start,
                                NULL, // stride
                                count,
                                NULL); // block
  HDF_check_result (hdf, result);

  // 7. Create memspace matching this rank's local buffer
  //    (shape of 'data' as seen locally)
  hsize_t mem_dims[H5S_MAX_RANK];
  for (int i = 0; i < rank; ++i)
    mem_dims[i] = local_size[i];

  hdf->mem_space = H5Screate_simple (rank, mem_dims, NULL);
  HDF_check_object (hdf, hdf->mem_space);

  // 8. Create transfer plist and enable collective MPI-IO
  hdf->xfer_plist = H5Pcreate (H5P_DATASET_XFER);
  HDF_check_object (hdf, hdf->xfer_plist);

  result = H5Pset_dxpl_mpio (hdf->xfer_plist, H5FD_MPIO_COLLECTIVE);
  HDF_check_result (hdf, result);

  // 9. Collective write into the appended region
  result =
    H5Dwrite (hdf->dset_id,
              dtype_id,        // in-memory type (e.g. H5T_NATIVE_DOUBLE)
              hdf->mem_space,  // local mem dataspace
              hdf->file_space, // file dataspace (with hyperslab selected)
              hdf->xfer_plist, // collective property list
              data);           // local buffer
  HDF_check_result (hdf, result);

  // 10. Close everything we opened
  H5Pclose (hdf->xfer_plist);
  hdf->xfer_plist = H5I_INVALID_HID;

  H5Sclose (hdf->mem_space);
  hdf->mem_space = H5I_INVALID_HID;

  H5Sclose (hdf->file_space);
  hdf->file_space = H5I_INVALID_HID;

  H5Dclose (hdf->dset_id);
  hdf->dset_id = H5I_INVALID_HID;
}

/**
 * @brief Write a collectively written dataset using compression with MPI-IO
 *
 * @param dataset_name The name of the dataset
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the data in
 * @param rank The rank of the dataset (e.g. the size of the dims array)
 * @param dims The dimensions of the dataset
 * @param max_dims The maximum dimensions of the dataset
 * @param chunk_dims The chunk dimensions
 * @param local_size The size of the sub-array this process writes into
 * @param local_offset The position of the sub-array this process writes into
 * @param compression_level The level of compression
 * @param hdf Pointer or reference to HDF object/struct
 *
 * @memberof HDF
 */

void HDF_collective_write_compressed_dataset (const char* dataset_name,
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
                                              HDF* hdf) {
  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the data space with dimensions and maximum dimensions
  hdf->dset_space_id = H5Screate_simple (rank, dims, max_dims);
  HDF_check_object (hdf, hdf->dset_space_id);

  // Create a property list for this dataset
  hdf->dcpl_id = H5Pcreate (H5P_DATASET_CREATE);
  HDF_check_object (hdf, hdf->dcpl_id);

  // Set the chunk size for the dataset
  result = H5Pset_chunk (hdf->dcpl_id, rank, chunk_dims);
  HDF_check_result (hdf, result);

  // Set compression on the dataset
  result = H5Pset_deflate (hdf->dcpl_id, compression_level);
  HDF_check_result (hdf, result);

  // Set the datatype
  hdf->dset_dtype_id = H5Tcopy (dtype_id);
  HDF_check_object (hdf, hdf->dset_dtype_id);

  // Create the dataset
  hdf->dset_id = H5Dcreate2 (group_id,
                             dataset_name,
                             hdf->dset_dtype_id,
                             hdf->dset_space_id,
                             H5P_DEFAULT,
                             hdf->dcpl_id,
                             H5P_DEFAULT);
  HDF_check_object (hdf, hdf->dset_id);

  // Create a property list for MPI-IO transfer
  hdf->xfer_plist = H5Pcreate (H5P_DATASET_XFER);
  HDF_check_object (hdf, hdf->xfer_plist);

  // Set MPIO_COLLECTIVE writing
  result = H5Pset_dxpl_mpio (hdf->xfer_plist, H5FD_MPIO_COLLECTIVE);
  HDF_check_result (hdf, result);

  // Get the dataset space
  hdf->file_space = H5Dget_space (hdf->dset_id);
  HDF_check_object (hdf, hdf->file_space);

  // Select a hyperslab for our process
  result = H5Sselect_hyperslab (
    hdf->file_space, H5S_SELECT_SET, local_offset, NULL, local_size, NULL);
  HDF_check_result (hdf, result);

  // Create a memory dataspace for our process
  hdf->mem_space = H5Screate_simple (rank, local_size, NULL);
  HDF_check_object (hdf, hdf->mem_space);

  // Actually do the collective write to the file
  result =
    H5Dwrite (hdf->dset_id,       /* dataset handle */
              hdf->dset_dtype_id, /* H5T_IEEE_F64LE */
              hdf->mem_space,     /* memory dataspace [local_nx] */
              hdf->file_space,    /* file dataspace with hyperslab selected */
              hdf->xfer_plist,    /* collective MPI‐IO transfer property */
              data                /* pointer to local data */
    );
  HDF_check_result (hdf, result);

  // Close opened objects
  H5Pclose (hdf->xfer_plist);
  hdf->xfer_plist = H5I_INVALID_HID;
  H5Sclose (hdf->mem_space);
  hdf->mem_space = H5I_INVALID_HID;
  H5Sclose (hdf->file_space);
  hdf->file_space = H5I_INVALID_HID;
  H5Tclose (hdf->dset_dtype_id);
  H5Pclose (hdf->dcpl_id);
  H5Sclose (hdf->dset_space_id);
  H5Dclose (hdf->dset_id);
}
