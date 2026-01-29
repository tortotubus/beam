#pragma once

#include <hdf5.h>
#include <mpi.h>

/**
 * @brief vtkHDF file type
 *
 * @memberof vtkHDF
 */
typedef enum {
  IMAGEDATA,
  POLYDATA,
  UNSTRUCTUREDGRID,
  HYPERTREEGRID,
  OVERLAPPINGAMR,
  PARTITIONEDDATASETCOLLECTION,
  MULTIBLOCKDATASET
} vtkHDFType;

/**
 * @brief Major and minor version of the vtkHDF file format
 *
 * @memberof vtkHDF
 */
typedef struct {
  int major;
  int minor;
} vtkHDFVersion;

/**
 * @brief Parent class with many helper functions for writing VTKHDF files
 */
typedef struct {
  hid_t file_id;       /**< File id */
  hid_t fapl;          /**< File access property list  */
  hid_t grp_vtkhdf_id; /**< VTKHDF group id */

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
} vtkHDF;

/**
 * @brief Close a vtkHDF file
 *
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 *
 * @memberof vtkHDF
 */
void vtk_HDF_close (vtkHDF* vtk_hdf) {
  if (vtk_hdf->grp_vtkhdf_id >= 0)
    H5Gclose (vtk_hdf->grp_vtkhdf_id);

  if (vtk_hdf->xfer_plist >= 0)
    H5Pclose (vtk_hdf->xfer_plist);

  if (vtk_hdf->mem_space >= 0)
    H5Sclose (vtk_hdf->mem_space);

  if (vtk_hdf->file_space >= 0)
    H5Sclose (vtk_hdf->file_space);

  if (vtk_hdf->fapl >= 0)
    H5Pclose (vtk_hdf->fapl);

  if (vtk_hdf->file_id >= 0)
    H5Fclose (vtk_hdf->file_id);
}

/**
 * @brief Function called on any error which immediately closes the file to
 * avoid corruptions
 *
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 *
 * @memberof vtkHDF
 */
void vtk_HDF_error (vtkHDF* vtk_hdf) {
  vtk_HDF_close (vtk_hdf);
}

/**
 * @brief Check for error of some operation that returns herr_t type
 *
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 * @param result Result of a HDF5 operation to check
 *
 * @memberof vtkHDF
 */
void vtk_HDF_check_result (vtkHDF* vtk_hdf, herr_t result) {
  if (result < 0) {
    vtk_HDF_error (vtk_hdf);
  }
}

/**
 * @brief Check for error of some operation that returns an object id
 *
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 * @param object_id Object id to check
 *
 * @memberof vtkHDF
 */
void vtk_HDF_check_object (vtkHDF* vtk_hdf, hid_t object_id) {
  if (object_id <= H5I_INVALID_HID) {
    vtk_HDF_error (vtk_hdf);
  }
}

/**
 * @brief Initialize the vtkHDF struct by writing a new file
 *
 * @memberof vtkHDF
 */
vtkHDF vtk_HDF_init (const char* fname, bool overwrite) {
  vtkHDF vtk_hdf = {
    .file_id = H5I_INVALID_HID,       /**< File id */
    .fapl = H5I_INVALID_HID,          /**< File access property list  */
    .grp_vtkhdf_id = H5I_INVALID_HID, /**< VTKHDF group id */
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
  vtk_hdf.file_id = H5Fcreate (fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  vtk_HDF_check_object (&vtk_hdf, vtk_hdf.file_id);

  // Create the VTKHDF group
  vtk_hdf.grp_vtkhdf_id = H5Gcreate2 (
    vtk_hdf.file_id, "/VTKHDF", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  vtk_HDF_check_object (&vtk_hdf, vtk_hdf.grp_vtkhdf_id);

  return vtk_hdf;
}

vtkHDF vtk_HDF_init_MPIIO (const char* fname, bool overwrite) {
  vtkHDF vtk_hdf = {
    .file_id = H5I_INVALID_HID,       /**< File id */
    .fapl = H5I_INVALID_HID,          /**< File access property list  */
    .grp_vtkhdf_id = H5I_INVALID_HID, /**< VTKHDF group id */
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
  vtk_hdf.fapl = H5Pcreate (H5P_FILE_ACCESS);
  if (vtk_hdf.fapl < 0)
    vtk_HDF_error (&vtk_hdf);

  // If we use MPI_SINGLE_FILE, we set all procs to have access and write into
  // the same file
  if (H5Pset_fapl_mpio (vtk_hdf.fapl, MPI_COMM_WORLD, MPI_INFO_NULL) < 0)
    vtk_HDF_error (&vtk_hdf);

  // Create the acutal file on disk
  vtk_hdf.file_id = H5Fcreate (fname, H5F_ACC_TRUNC, H5P_DEFAULT, vtk_hdf.fapl);
  if (vtk_hdf.file_id < 0)
    vtk_HDF_error (&vtk_hdf);

  // Create the VTKHDF group
  vtk_hdf.grp_vtkhdf_id = H5Gcreate2 (
    vtk_hdf.file_id, "/VTKHDF", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (vtk_hdf.grp_vtkhdf_id < 0)
    vtk_HDF_error (&vtk_hdf);

  return vtk_hdf;
}

/**
 * @brief Initialize the vtkHDF struct by opening an existing file
 *
 * @memberof vtkHDF
 */
vtkHDF vtk_HDF_open (const char* fname) {
  vtkHDF vtk_hdf = {
    .file_id = H5I_INVALID_HID,       /**< File id */
    .fapl = H5I_INVALID_HID,          /**< File access property list  */
    .grp_vtkhdf_id = H5I_INVALID_HID, /**< VTKHDF group id */
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
  vtk_hdf.file_id = H5Fopen (fname, H5F_ACC_RDWR, H5P_DEFAULT);
  vtk_HDF_check_object (&vtk_hdf, vtk_hdf.file_id);

  // Open the VTKHDF group
  vtk_hdf.grp_vtkhdf_id = H5Gopen (vtk_hdf.file_id, "/VTKHDF", H5P_DEFAULT);
  vtk_HDF_check_object (&vtk_hdf, vtk_hdf.grp_vtkhdf_id);

  return vtk_hdf;
}

/**
 * @brief Initialize the vtkHDF struct by opening an existing file
 *
 * @memberof vtkHDF
 */
vtkHDF vtk_HDF_open_MPIIO (const char* fname) {
  vtkHDF vtk_hdf = {
    .file_id = H5I_INVALID_HID,       /**< File id */
    .fapl = H5I_INVALID_HID,          /**< File access property list  */
    .grp_vtkhdf_id = H5I_INVALID_HID, /**< VTKHDF group id */
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
  vtk_hdf.fapl = H5Pcreate (H5P_FILE_ACCESS);
  vtk_HDF_check_object (&vtk_hdf, vtk_hdf.fapl);

  // If we use MPI_SINGLE_FILE, we set all procs to have access and write into
  // the same file
  result = H5Pset_fapl_mpio (vtk_hdf.fapl, MPI_COMM_WORLD, MPI_INFO_NULL);
  vtk_HDF_check_result (&vtk_hdf, result);

  // Open the actual file on disk
  vtk_hdf.file_id = H5Fopen (fname, H5F_ACC_RDWR, vtk_hdf.fapl);
  vtk_HDF_check_object (&vtk_hdf, vtk_hdf.file_id);

  // Open the VTKHDF group
  vtk_hdf.grp_vtkhdf_id = H5Gopen (vtk_hdf.file_id, "/VTKHDF", H5P_DEFAULT);
  vtk_HDF_check_object (&vtk_hdf, vtk_hdf.grp_vtkhdf_id);

  return vtk_hdf;
}

/**
 * @brief Write new attribute
 *
 * @param attribute_name The name of the attrbiute
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the attribute in
 * @param dims The dimensions of the data
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 *
 * @memberof vtkHDF
 */
void vtk_HDF_write_attribute (const char* attribute_name,
                              const void* data,
                              const hid_t dtype_id,
                              const hid_t group_id,
                              const hsize_t dims[],
                              vtkHDF* vtk_hdf) {
  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the attribute space
  vtk_hdf->attr_space_id = H5Screate_simple (1, dims, dims);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->attr_space_id);

  // Set the datatype
  vtk_hdf->attr_dtype_id = H5Tcopy (dtype_id);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->attr_dtype_id);

  // Create the attribute
  vtk_hdf->attr_id = H5Acreate2 (group_id,
                                 attribute_name,
                                 vtk_hdf->attr_dtype_id,
                                 vtk_hdf->attr_space_id,
                                 H5P_DEFAULT,
                                 H5P_DEFAULT);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->attr_id);

  // Write the attribute
  result = H5Awrite (vtk_hdf->attr_id, vtk_hdf->attr_dtype_id, data);
  vtk_HDF_check_object (vtk_hdf, result);

  // Close the open objects
  H5Aclose (vtk_hdf->attr_id);
  H5Tclose (vtk_hdf->attr_dtype_id);
  H5Sclose (vtk_hdf->attr_space_id);
}

/**
 * @brief Write new scalar attribute
 *
 * @param attribute_name The name of the attrbiute
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the attribute in
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 *
 * @memberof vtkHDF
 */

void vtk_HDF_write_scalar_attribute (const char* attribute_name,
                                     const void* data,
                                     const hid_t dtype_id,
                                     const hid_t group_id,
                                     vtkHDF* vtk_hdf) {
  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the attribute space
  vtk_hdf->attr_space_id = H5Screate (H5S_SCALAR);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->attr_space_id);

  // Set the datatype
  vtk_hdf->attr_dtype_id = H5Tcopy (dtype_id);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->attr_dtype_id);

  // Create the attribute
  vtk_hdf->attr_id = H5Acreate2 (group_id,
                                 attribute_name,
                                 vtk_hdf->attr_dtype_id,
                                 vtk_hdf->attr_space_id,
                                 H5P_DEFAULT,
                                 H5P_DEFAULT);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->attr_id);

  // Write the attribute
  result = H5Awrite (vtk_hdf->attr_id, vtk_hdf->attr_dtype_id, data);
  vtk_HDF_check_object (vtk_hdf, result);

  // Close the open objects
  H5Aclose (vtk_hdf->attr_id);
  H5Tclose (vtk_hdf->attr_dtype_id);
  H5Sclose (vtk_hdf->attr_space_id);
}

/**
 * @brief Modify existing scalar attribute
 *
 * @param attribute_name The name of the attrbiute
 * @param data Pointer to the array of data
 * @param dtype_id The HDF5 datatype ID
 * @param group_id The HDF5 group to write the attribute in
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 *
 * @memberof vtkHDF
 */
void vtk_HDF_modify_scalar_attribute (const char* attribute_name,
                                      const void* data,
                                      const hid_t dtype_id,
                                      const hid_t group_id,
                                      vtkHDF* vtk_hdf) {
  herr_t result = 0;

  // 1. Open existing attribute
  vtk_hdf->attr_id = H5Aopen (group_id, attribute_name, H5P_DEFAULT);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->attr_id);

  // 2. (Optional but nice) get its dataspace and verify it is scalar
  vtk_hdf->attr_space_id = H5Aget_space (vtk_hdf->attr_id);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->attr_space_id);

  int ndims = H5Sget_simple_extent_ndims (vtk_hdf->attr_space_id);
  vtk_HDF_check_result (vtk_hdf, ndims);
  if (ndims != 0) {
    perror ("vtk_HDF_modify_scalar_attribute: attribute is not scalar.\n");
    exit (1);
  }

  // 3. (Optional) get file datatype if you want to sanity-check or inspect it
  vtk_hdf->attr_dtype_id = H5Aget_type (vtk_hdf->attr_id);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->attr_dtype_id);
  // You *could* compare mem_dtype_id vs attr_dtype_id here if desired.

  // 4. Overwrite the attribute value in-place
  //    mem_dtype_id is the in-memory type (often H5T_NATIVE_*)
  result = H5Awrite (vtk_hdf->attr_id, dtype_id, data);
  vtk_HDF_check_result (vtk_hdf, result);

  // 5. Close only what we actually opened
  H5Tclose (vtk_hdf->attr_dtype_id); // from H5Aget_type
  H5Sclose (vtk_hdf->attr_space_id); // from H5Aget_space
  H5Aclose (vtk_hdf->attr_id);       // from H5Aopen
}

/**
 * @brief Read scalar attribute
 *
 * @param attribute_name The name of the attrbiute
 * @param data Pointer to the array of data
 * @param group_id The HDF5 group to write the attribute in
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 *
 * @memberof vtkHDF
 */
void vtk_HDF_read_scalar_attribute (const char* attribute_name,
                                    void* data,
                                    const hid_t group_id,
                                    vtkHDF* vtk_hdf) {

  herr_t result;

  // Open the attribute
  vtk_hdf->attr_id = H5Aopen (group_id, attribute_name, H5P_DEFAULT);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->attr_id);

  // Get the datatype
  vtk_hdf->attr_dtype_id = H5Aget_type (vtk_hdf->attr_id);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->attr_dtype_id);

  // Get the dataspace
  vtk_hdf->attr_space_id = H5Aget_space (vtk_hdf->attr_id);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->attr_space_id);

  // Read the scalar value into data pointer
  result = H5Aread (vtk_hdf->attr_id, vtk_hdf->attr_dtype_id, data);
  vtk_HDF_check_result (vtk_hdf, result);

  H5Aclose (vtk_hdf->attr_id);
  H5Tclose (vtk_hdf->attr_dtype_id);
  H5Sclose (vtk_hdf->attr_space_id);
}

/**
 * @brief Write new attribute
 *
 * @param type_name The name of the attrbiute
 * @param group_id The HDF5 group to write the attribute in
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 *
 * @memberof vtkHDF
 */
void vtk_HDF_write_type_attribute (const char* type_name,
                                   const hid_t group_id,
                                   vtkHDF* vtk_hdf) {
  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the attribute space
  vtk_hdf->attr_space_id = H5Screate (H5S_SCALAR);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->attr_space_id);

  // Set the datatype
  vtk_hdf->attr_dtype_id = H5Tcopy (H5T_C_S1);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->attr_dtype_id);

  // Set the size
  result = H5Tset_size (vtk_hdf->attr_dtype_id, strlen (type_name));
  vtk_HDF_check_result (vtk_hdf, result);

  // Pad the string
  result = H5Tset_strpad (vtk_hdf->attr_dtype_id, H5T_STR_NULLTERM);
  vtk_HDF_check_result (vtk_hdf, result);

  // Set the string value
  result = H5Tset_cset (vtk_hdf->attr_dtype_id, H5T_CSET_ASCII);
  vtk_HDF_check_result (vtk_hdf, result);

  // Create the attribute
  vtk_hdf->attr_id = H5Acreate2 (group_id,
                                 "Type",
                                 vtk_hdf->attr_dtype_id,
                                 vtk_hdf->attr_space_id,
                                 H5P_DEFAULT,
                                 H5P_DEFAULT);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->attr_id);

  // Write the attribute
  result = H5Awrite (vtk_hdf->attr_id, vtk_hdf->attr_dtype_id, type_name);
  vtk_HDF_check_object (vtk_hdf, result);

  // Close the open objects
  H5Aclose (vtk_hdf->attr_id);
  H5Tclose (vtk_hdf->attr_dtype_id);
  H5Sclose (vtk_hdf->attr_space_id);
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
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 *
 * @memberof vtkHDF
 */
void vtk_HDF_read_dataset (const char* dataset_name,
                           void** data,
                           const hid_t dtype_id,
                           const hid_t group_id,
                           const int rank,
                           hsize_t dims[],
                           hsize_t max_dims[],
                           vtkHDF* vtk_hdf) {
  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Open the dataset
  vtk_hdf->dset_id = H5Dopen2 (group_id, dataset_name, H5P_DEFAULT);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dset_id);

  // Dataspace
  vtk_hdf->dset_space_id = H5Dget_space (vtk_hdf->dset_id);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dset_space_id);

  // Rank
  int file_rank = H5Sget_simple_extent_ndims (vtk_hdf->dset_space_id);

  if (file_rank != rank) {
    perror ("vtk_HDF_read_dataset : Rank mismatch.\n");
    exit (1);
  }

  // Dims / max_dims
  result = H5Sget_simple_extent_dims (vtk_hdf->dset_space_id, dims, max_dims);

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
    H5Dread (vtk_hdf->dset_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, *data);
  vtk_HDF_check_result (vtk_hdf, result);

  // Close opened objects
  H5Dclose (vtk_hdf->dset_id);
  H5Sclose (vtk_hdf->dset_space_id);
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
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 *
 * @memberof vtkHDF
 */
void vtk_HDF_read_dataset_dims (const char* dataset_name, 
                                const hid_t dtype_id,
                                const hid_t group_id,
                                const int rank,
                                hsize_t dims[],
                                hsize_t max_dims[],
                                vtkHDF* vtk_hdf) {
  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Open the dataset
  vtk_hdf->dset_id = H5Dopen2 (group_id, dataset_name, H5P_DEFAULT);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dset_id);

  // Dataspace
  vtk_hdf->dset_space_id = H5Dget_space (vtk_hdf->dset_id);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dset_space_id);

  // Rank
  int file_rank = H5Sget_simple_extent_ndims (vtk_hdf->dset_space_id);

  if (file_rank != rank) {
    perror ("vtk_HDF_read_dataset_dims: Rank mismatch.\n");
    exit (1);
  }

  // Dims / max_dims
  result = H5Sget_simple_extent_dims (vtk_hdf->dset_space_id, dims, max_dims);

  // Close opened objects
  H5Dclose (vtk_hdf->dset_id);
  H5Sclose (vtk_hdf->dset_space_id);
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
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 *
 * @memberof vtkHDF
 */
void vtk_HDF_write_dataset (const char* dataset_name,
                            const void* data,
                            const hid_t dtype_id,
                            const hid_t group_id,
                            const int rank,
                            const hsize_t dims[],
                            vtkHDF* vtk_hdf) {
  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the data space with dimensions and maximum dimensions
  vtk_hdf->dset_space_id = H5Screate_simple (rank, dims, dims);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dset_space_id);

  // Create a property list for this dataset
  vtk_hdf->dcpl_id = H5Pcreate (H5P_DATASET_CREATE);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dcpl_id);

  // Set the datatype
  vtk_hdf->dset_dtype_id = H5Tcopy (dtype_id);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dset_dtype_id);

  // Create the dataset
  vtk_hdf->dset_id = H5Dcreate2 (group_id,
                                 dataset_name,
                                 vtk_hdf->dset_dtype_id,
                                 vtk_hdf->dset_space_id,
                                 H5P_DEFAULT,
                                 vtk_hdf->dcpl_id,
                                 H5P_DEFAULT);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dset_id);

  result = H5Dwrite (vtk_hdf->dset_id,
                     vtk_hdf->dset_dtype_id,
                     H5S_ALL,
                     H5S_ALL,
                     H5P_DEFAULT,
                     data);
  vtk_HDF_check_result (vtk_hdf, result);

  // Close opened objects
  H5Dclose (vtk_hdf->dset_id);
  H5Tclose (vtk_hdf->dset_dtype_id);
  H5Pclose (vtk_hdf->dcpl_id);
  H5Sclose (vtk_hdf->dset_space_id);
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
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 */
void vtk_HDF_write_chunked_dataset (const char* dataset_name,
                                    const void* data,
                                    const hid_t dtype_id,
                                    const hid_t group_id,
                                    const int rank,
                                    const hsize_t dims[],
                                    const hsize_t max_dims[],
                                    const hsize_t chunk_dims[],
                                    vtkHDF* vtk_hdf) {

  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the data space with dimensions and maximum dimensions
  vtk_hdf->dset_space_id = H5Screate_simple (rank, dims, max_dims);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dset_space_id);

  // Create a property list for this dataset
  vtk_hdf->dcpl_id = H5Pcreate (H5P_DATASET_CREATE);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dcpl_id);

  // Set the chunk size for the dataset
  result = H5Pset_chunk (vtk_hdf->dcpl_id, rank, chunk_dims);
  vtk_HDF_check_result (vtk_hdf, result);

  // Set the datatype
  vtk_hdf->dset_dtype_id = H5Tcopy (dtype_id);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dset_dtype_id);

  // Create the dataset
  vtk_hdf->dset_id = H5Dcreate2 (group_id,
                                 dataset_name,
                                 vtk_hdf->dset_dtype_id,
                                 vtk_hdf->dset_space_id,
                                 H5P_DEFAULT,
                                 vtk_hdf->dcpl_id,
                                 H5P_DEFAULT);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dset_id);

  result = H5Dwrite (vtk_hdf->dset_id,
                     vtk_hdf->dset_dtype_id,
                     H5S_ALL,
                     H5S_ALL,
                     H5P_DEFAULT,
                     data);
  vtk_HDF_check_result (vtk_hdf, result);

  // Close opened objects
  H5Dclose (vtk_hdf->dset_id);
  H5Tclose (vtk_hdf->dset_dtype_id);
  H5Pclose (vtk_hdf->dcpl_id);
  H5Sclose (vtk_hdf->dset_space_id);
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
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 */

void vtk_HDF_append_chunked_dataset (const char* dataset_name,
                                     const void* data,
                                     const hid_t dtype_id,
                                     const hid_t group_id,
                                     const int rank,
                                     const hsize_t dims[],
                                     vtkHDF* vtk_hdf) {
  herr_t result = 0;

  // Open existing dataset
  vtk_hdf->dset_id = H5Dopen2 (group_id, dataset_name, H5P_DEFAULT);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dset_id);

  // Get current dataspace and its dimensions
  vtk_hdf->dset_space_id = H5Dget_space (vtk_hdf->dset_id);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dset_space_id);

  hsize_t cur_dims[H5S_MAX_RANK];
  hsize_t cur_max_dims[H5S_MAX_RANK];
  int file_rank =
    H5Sget_simple_extent_dims (vtk_hdf->dset_space_id, cur_dims, cur_max_dims);
  vtk_HDF_check_result (vtk_hdf, file_rank);

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
  result = H5Dset_extent (vtk_hdf->dset_id, new_dims);
  vtk_HDF_check_result (vtk_hdf, result);

  // After extending, the old dataspace is invalid; close and re-open
  H5Sclose (vtk_hdf->dset_space_id);
  vtk_hdf->dset_space_id = H5Dget_space (vtk_hdf->dset_id);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dset_space_id);

  // 4. Select hyperslab corresponding to the newly appended region
  hsize_t start[H5S_MAX_RANK];
  hsize_t count[H5S_MAX_RANK];

  start[0] = cur_dims[0]; // start right after old end
  count[0] = dims[0];     // append this many along dim 0

  for (int i = 1; i < rank; ++i) {
    start[i] = 0;
    count[i] = dims[i]; // whole extent in other dims
  }

  result = H5Sselect_hyperslab (vtk_hdf->dset_space_id,
                                H5S_SELECT_SET,
                                start,
                                NULL, // stride
                                count,
                                NULL); // block

  vtk_HDF_check_result (vtk_hdf, result);

  // Create memspace for the data block we’re appending
  vtk_hdf->mem_space = H5Screate_simple (rank, dims, NULL);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->mem_space);

  // Write the new chunk into the extended portion
  result = H5Dwrite (vtk_hdf->dset_id,
                     dtype_id,               // in-memory type
                     vtk_hdf->mem_space,     // mem dataspace selection
                     vtk_hdf->dset_space_id, // file dataspace selection
                     H5P_DEFAULT,            // xfer plist
                     data);
  vtk_HDF_check_result (vtk_hdf, result);

  // Close everything
  H5Sclose (vtk_hdf->mem_space);
  H5Sclose (vtk_hdf->dset_space_id);
  H5Dclose (vtk_hdf->dset_id);
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
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 *
 * @memberof vtkHDF
 */
void vtk_HDF_write_compressed_dataset (const char* dataset_name,
                                       const void* data,
                                       const hid_t dtype_id,
                                       const hid_t group_id,
                                       const int rank,
                                       const hsize_t dims[],
                                       const hsize_t max_dims[],
                                       const hsize_t chunk_dims[],
                                       unsigned int compression_level,
                                       vtkHDF* vtk_hdf) {

  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the data space with dimensions and maximum dimensions
  vtk_hdf->dset_space_id = H5Screate_simple (rank, dims, max_dims);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dset_space_id);

  // Create a property list for this dataset
  vtk_hdf->dcpl_id = H5Pcreate (H5P_DATASET_CREATE);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dcpl_id);

  // Set the chunk size for the dataset
  result = H5Pset_chunk (vtk_hdf->dcpl_id, rank, chunk_dims);
  vtk_HDF_check_result (vtk_hdf, result);

  // Set compression on the dataset
  result = H5Pset_deflate (vtk_hdf->dcpl_id, compression_level);
  vtk_HDF_check_result (vtk_hdf, result);

  // Set the datatype
  vtk_hdf->dset_dtype_id = H5Tcopy (dtype_id);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dset_dtype_id);

  // Create the dataset
  vtk_hdf->dset_id = H5Dcreate2 (group_id,
                                 dataset_name,
                                 vtk_hdf->dset_dtype_id,
                                 vtk_hdf->dset_space_id,
                                 H5P_DEFAULT,
                                 vtk_hdf->dcpl_id,
                                 H5P_DEFAULT);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dset_id);

  // Write the dataset
  result = H5Dwrite (vtk_hdf->dset_id,
                     vtk_hdf->dset_dtype_id,
                     H5S_ALL,
                     H5S_ALL,
                     H5P_DEFAULT,
                     data);
  vtk_HDF_check_result (vtk_hdf, result);

  // Close opened objects
  H5Dclose (vtk_hdf->dset_id);
  H5Tclose (vtk_hdf->dset_dtype_id);
  H5Pclose (vtk_hdf->dcpl_id);
  H5Sclose (vtk_hdf->dset_space_id);
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
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 *
 * @memberof vtkHDF
 */
void vtk_HDF_collective_write_dataset (const char* dataset_name,
                                       const void* data,
                                       const hid_t dtype_id,
                                       const hid_t group_id,
                                       const int rank,
                                       const hsize_t dims[],
                                       const hsize_t local_size[],
                                       const hsize_t local_offset[],
                                       vtkHDF* vtk_hdf) {

  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the data space with dimensions and maximum dimensions
  vtk_hdf->dset_space_id = H5Screate_simple (rank, dims, dims);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dset_space_id);

  // Set the datatype
  vtk_hdf->dset_dtype_id = H5Tcopy (dtype_id);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dset_dtype_id);

  // Create the dataset
  vtk_hdf->dset_id = H5Dcreate2 (group_id,
                                 dataset_name,
                                 vtk_hdf->dset_dtype_id,
                                 vtk_hdf->dset_space_id,
                                 H5P_DEFAULT,
                                 H5P_DEFAULT,
                                 H5P_DEFAULT);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dset_id);

  // Create a property list for MPI-IO transfer
  vtk_hdf->xfer_plist = H5Pcreate (H5P_DATASET_XFER);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->xfer_plist);

  // Set MPIO_COLLECTIVE writing
  result = H5Pset_dxpl_mpio (vtk_hdf->xfer_plist, H5FD_MPIO_COLLECTIVE);
  vtk_HDF_check_result (vtk_hdf, result);

  // Get the dataset space
  vtk_hdf->file_space = H5Dget_space (vtk_hdf->dset_id);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->file_space);

  // Select a hyperslab for our process
  result = H5Sselect_hyperslab (
    vtk_hdf->file_space, H5S_SELECT_SET, local_offset, NULL, local_size, NULL);
  vtk_HDF_check_result (vtk_hdf, result);

  // Create a memory dataspace for our process
  vtk_hdf->mem_space = H5Screate_simple (rank, local_size, NULL);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->mem_space);

  // Actually do the collective write to the file
  result =
    H5Dwrite (vtk_hdf->dset_id,       /* dataset handle */
              vtk_hdf->dset_dtype_id, /* H5T_IEEE_F64LE */
              vtk_hdf->mem_space,     /* memory dataspace [local_nx] */
              vtk_hdf->file_space, /* file dataspace with hyperslab selected */
              vtk_hdf->xfer_plist, /* collective MPI‐IO transfer property */
              data                 /* pointer to local data */
    );
  vtk_HDF_check_result (vtk_hdf, result);

  // Close opened objects
  H5Pclose (vtk_hdf->xfer_plist);
  vtk_hdf->xfer_plist = H5I_INVALID_HID;
  H5Sclose (vtk_hdf->mem_space);
  vtk_hdf->mem_space = H5I_INVALID_HID;
  H5Sclose (vtk_hdf->file_space);
  vtk_hdf->file_space = H5I_INVALID_HID;
  H5Tclose (vtk_hdf->dset_dtype_id);
  H5Sclose (vtk_hdf->dset_space_id);
  H5Dclose (vtk_hdf->dset_id);
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
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 *
 * @memberof vtkHDF
 */
void vtk_HDF_collective_write_chunked_dataset (const char* dataset_name,
                                               const void* data,
                                               const hid_t dtype_id,
                                               const hid_t group_id,
                                               const int rank,
                                               const hsize_t dims[],
                                               const hsize_t max_dims[],
                                               const hsize_t chunk_dims[],
                                               const hsize_t local_size[],
                                               const hsize_t local_offset[],
                                               vtkHDF* vtk_hdf) {

  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the data space with dimensions and maximum dimensions
  vtk_hdf->dset_space_id = H5Screate_simple (rank, dims, max_dims);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dset_space_id);

  // Create a property list for this dataset
  vtk_hdf->dcpl_id = H5Pcreate (H5P_DATASET_CREATE);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dcpl_id);

  // Set the chunk size for the dataset
  result = H5Pset_chunk (vtk_hdf->dcpl_id, rank, chunk_dims);
  vtk_HDF_check_result (vtk_hdf, result);

  // Set the datatype
  vtk_hdf->dset_dtype_id = H5Tcopy (dtype_id);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dset_dtype_id);

  // Create the dataset
  vtk_hdf->dset_id = H5Dcreate2 (group_id,
                                 dataset_name,
                                 vtk_hdf->dset_dtype_id,
                                 vtk_hdf->dset_space_id,
                                 H5P_DEFAULT,
                                 vtk_hdf->dcpl_id,
                                 H5P_DEFAULT);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dset_id);

  // Create a property list for MPI-IO transfer
  vtk_hdf->xfer_plist = H5Pcreate (H5P_DATASET_XFER);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->xfer_plist);

  // Set MPIO_COLLECTIVE writing
  result = H5Pset_dxpl_mpio (vtk_hdf->xfer_plist, H5FD_MPIO_COLLECTIVE);
  vtk_HDF_check_result (vtk_hdf, result);

  // Get the dataset space
  vtk_hdf->file_space = H5Dget_space (vtk_hdf->dset_id);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->file_space);

  // Select a hyperslab for our process
  result = H5Sselect_hyperslab (
    vtk_hdf->file_space, H5S_SELECT_SET, local_offset, NULL, local_size, NULL);
  vtk_HDF_check_result (vtk_hdf, result);

  // Create a memory dataspace for our process
  vtk_hdf->mem_space = H5Screate_simple (rank, local_size, NULL);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->mem_space);

  // Actually do the collective write to the file
  result =
    H5Dwrite (vtk_hdf->dset_id,       /* dataset handle */
              vtk_hdf->dset_dtype_id, /* H5T_IEEE_F64LE */
              vtk_hdf->mem_space,     /* memory dataspace */
              vtk_hdf->file_space, /* file dataspace with hyperslab selected */
              vtk_hdf->xfer_plist, /* collective MPI‐IO transfer property */
              data                 /* pointer to local data */
    );
  vtk_HDF_check_result (vtk_hdf, result);

  // Close opened objects
  H5Pclose (vtk_hdf->xfer_plist);
  vtk_hdf->xfer_plist = H5I_INVALID_HID;
  H5Sclose (vtk_hdf->mem_space);
  vtk_hdf->mem_space = H5I_INVALID_HID;
  H5Sclose (vtk_hdf->file_space);
  vtk_hdf->file_space = H5I_INVALID_HID;
  H5Tclose (vtk_hdf->dset_dtype_id);
  H5Pclose (vtk_hdf->dcpl_id);
  H5Sclose (vtk_hdf->dset_space_id);
  H5Dclose (vtk_hdf->dset_id);
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
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 */
void vtk_HDF_collective_append_chunked_dataset (
  const char* dataset_name,
  const void* data,
  const hid_t dtype_id,
  const hid_t group_id,
  const int rank,
  const hsize_t dims[],         // global size of appended block
  const hsize_t local_size[],   // local tile within appended block
  const hsize_t local_offset[], // local offset within appended block
  vtkHDF* vtk_hdf) {
  herr_t result = 0;

  // 1. Open existing dataset collectively
  vtk_hdf->dset_id = H5Dopen2 (group_id, dataset_name, H5P_DEFAULT);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dset_id);

  // 2. Get current global dims
  vtk_hdf->dset_space_id = H5Dget_space (vtk_hdf->dset_id);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dset_space_id);

  hsize_t cur_dims[H5S_MAX_RANK];
  hsize_t cur_max_dims[H5S_MAX_RANK];
  int file_rank =
    H5Sget_simple_extent_dims (vtk_hdf->dset_space_id, cur_dims, cur_max_dims);
  if (file_rank < 0)
    vtk_HDF_error (vtk_hdf);
  if (file_rank != rank) {
    fprintf (stderr,
             "vtk_HDF_collective_append_chunked_dataset: rank mismatch "
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
  result = H5Dset_extent (vtk_hdf->dset_id, new_dims);
  vtk_HDF_check_result (vtk_hdf, result);

  // Old dataspace no longer valid
  H5Sclose (vtk_hdf->dset_space_id);
  vtk_hdf->dset_space_id = H5I_INVALID_HID;

  // 5. Get new file dataspace for selection
  vtk_hdf->file_space = H5Dget_space (vtk_hdf->dset_id);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->file_space);

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

  result = H5Sselect_hyperslab (vtk_hdf->file_space,
                                H5S_SELECT_SET,
                                start,
                                NULL, // stride
                                count,
                                NULL); // block
  vtk_HDF_check_result (vtk_hdf, result);

  // 7. Create memspace matching this rank's local buffer
  //    (shape of 'data' as seen locally)
  hsize_t mem_dims[H5S_MAX_RANK];
  for (int i = 0; i < rank; ++i)
    mem_dims[i] = local_size[i];

  vtk_hdf->mem_space = H5Screate_simple (rank, mem_dims, NULL);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->mem_space);

  // 8. Create transfer plist and enable collective MPI-IO
  vtk_hdf->xfer_plist = H5Pcreate (H5P_DATASET_XFER);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->xfer_plist);

  result = H5Pset_dxpl_mpio (vtk_hdf->xfer_plist, H5FD_MPIO_COLLECTIVE);
  vtk_HDF_check_result (vtk_hdf, result);

  // 9. Collective write into the appended region
  result =
    H5Dwrite (vtk_hdf->dset_id,
              dtype_id,            // in-memory type (e.g. H5T_NATIVE_DOUBLE)
              vtk_hdf->mem_space,  // local mem dataspace
              vtk_hdf->file_space, // file dataspace (with hyperslab selected)
              vtk_hdf->xfer_plist, // collective property list
              data);               // local buffer
  vtk_HDF_check_result (vtk_hdf, result);

  // 10. Close everything we opened
  H5Pclose (vtk_hdf->xfer_plist);
  vtk_hdf->xfer_plist = H5I_INVALID_HID;

  H5Sclose (vtk_hdf->mem_space);
  vtk_hdf->mem_space = H5I_INVALID_HID;

  H5Sclose (vtk_hdf->file_space);
  vtk_hdf->file_space = H5I_INVALID_HID;

  H5Dclose (vtk_hdf->dset_id);
  vtk_hdf->dset_id = H5I_INVALID_HID;
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
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 *
 * @memberof vtkHDF
 */

void vtk_HDF_collective_write_compressed_dataset (
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
  vtkHDF* vtk_hdf) {
  // Store the result of HDF5 functions that do not return an identifier
  herr_t result = 0;

  // Create the data space with dimensions and maximum dimensions
  vtk_hdf->dset_space_id = H5Screate_simple (rank, dims, max_dims);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dset_space_id);

  // Create a property list for this dataset
  vtk_hdf->dcpl_id = H5Pcreate (H5P_DATASET_CREATE);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dcpl_id);

  // Set the chunk size for the dataset
  result = H5Pset_chunk (vtk_hdf->dcpl_id, rank, chunk_dims);
  vtk_HDF_check_result (vtk_hdf, result);

  // Set compression on the dataset
  result = H5Pset_deflate (vtk_hdf->dcpl_id, compression_level);
  vtk_HDF_check_result (vtk_hdf, result);

  // Set the datatype
  vtk_hdf->dset_dtype_id = H5Tcopy (dtype_id);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dset_dtype_id);

  // Create the dataset
  vtk_hdf->dset_id = H5Dcreate2 (group_id,
                                 dataset_name,
                                 vtk_hdf->dset_dtype_id,
                                 vtk_hdf->dset_space_id,
                                 H5P_DEFAULT,
                                 vtk_hdf->dcpl_id,
                                 H5P_DEFAULT);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->dset_id);

  // Create a property list for MPI-IO transfer
  vtk_hdf->xfer_plist = H5Pcreate (H5P_DATASET_XFER);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->xfer_plist);

  // Set MPIO_COLLECTIVE writing
  result = H5Pset_dxpl_mpio (vtk_hdf->xfer_plist, H5FD_MPIO_COLLECTIVE);
  vtk_HDF_check_result (vtk_hdf, result);

  // Get the dataset space
  vtk_hdf->file_space = H5Dget_space (vtk_hdf->dset_id);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->file_space);

  // Select a hyperslab for our process
  result = H5Sselect_hyperslab (
    vtk_hdf->file_space, H5S_SELECT_SET, local_offset, NULL, local_size, NULL);
  vtk_HDF_check_result (vtk_hdf, result);

  // Create a memory dataspace for our process
  vtk_hdf->mem_space = H5Screate_simple (rank, local_size, NULL);
  vtk_HDF_check_object (vtk_hdf, vtk_hdf->mem_space);

  // Actually do the collective write to the file
  result =
    H5Dwrite (vtk_hdf->dset_id,       /* dataset handle */
              vtk_hdf->dset_dtype_id, /* H5T_IEEE_F64LE */
              vtk_hdf->mem_space,     /* memory dataspace [local_nx] */
              vtk_hdf->file_space, /* file dataspace with hyperslab selected */
              vtk_hdf->xfer_plist, /* collective MPI‐IO transfer property */
              data                 /* pointer to local data */
    );
  vtk_HDF_check_result (vtk_hdf, result);

  // Close opened objects
  H5Pclose (vtk_hdf->xfer_plist);
  vtk_hdf->xfer_plist = H5I_INVALID_HID;
  H5Sclose (vtk_hdf->mem_space);
  vtk_hdf->mem_space = H5I_INVALID_HID;
  H5Sclose (vtk_hdf->file_space);
  vtk_hdf->file_space = H5I_INVALID_HID;
  H5Tclose (vtk_hdf->dset_dtype_id);
  H5Pclose (vtk_hdf->dcpl_id);
  H5Sclose (vtk_hdf->dset_space_id);
  H5Dclose (vtk_hdf->dset_id);
}
