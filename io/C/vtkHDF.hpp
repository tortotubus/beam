#pragma once

#include "config/config.hpp"
#include "general/error.hpp"

#include <cstring>
#include <string>

#include <hdf5.h>
#include <mpi.h>

namespace ELFF {
namespace io {
namespace C {

/**
 * @brief vtkHDF file type
 * 
 * @memberof vtkHDF
 */
typedef enum
{
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
typedef struct
{
  int major;
  int minor;
} vtkHDFVersion;

/**
 * @brief Parent class with many helper functions for writing VTKHDF files
 */
typedef struct
{
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
void
vtk_HDF_close(vtkHDF* vtk_hdf);

/**
 * @brief Function called on any error which immediately closes the file to
 * avoid corruptions
 *
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 * 
 * @memberof vtkHDF
 */
void
vtk_HDF_error(vtkHDF* vtk_hdf);

/**
 * @brief Check for error of some operation that returns herr_t type
 *
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 * @param result Result of a HDF5 operation to check
 *
 * @memberof vtkHDF
 */
void
vtk_HDF_check_result(vtkHDF* vtk_hdf, herr_t result);

/**
 * @brief Check for error of some operation that returns an object id
 *
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 * @param object_id Object id to check
 *
 * @memberof vtkHDF
 */
void
vtk_HDF_check_object(vtkHDF* vtk_hdf, hid_t object_id);

/**
 * @brief Initialize the vtkHDF struct by writing a new file
 *
 * @memberof vtkHDF
 */
vtkHDF
vtk_HDF_init(const char* fname, bool overwrite);

/**
 * @brief Initialize the vtkHDF struct by opening an existing file
 *
 * @memberof vtkHDF
 */
vtkHDF
vtk_HDF_open(const char* fname);

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
void
vtk_HDF_write_attribute(const char* attribute_name,
                        const void* data,
                        const hid_t dtype_id,
                        const hid_t group_id,
                        const hsize_t dims[],
                        vtkHDF* vtk_hdf);

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

void
vtk_HDF_write_scalar_attribute(const char* attribute_name,
                               const void* data,
                               const hid_t dtype_id,
                               const hid_t group_id,
                               vtkHDF* vtk_hdf);

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
void
vtk_HDF_modify_scalar_attribute(const char* attribute_name,
                                const void* data,
                                const hid_t dtype_id,
                                const hid_t group_id,
                                vtkHDF* vtk_hdf);
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

void
vtk_HDF_read_scalar_attribute(const char* attribute_name,
                              void* data,
                              const hid_t group_id,
                              vtkHDF* vtk_hdf);

/**
 * @brief Write new attribute
 *
 * @param type_name The name of the attrbiute
 * @param group_id The HDF5 group to write the attribute in
 * @param vtk_hdf Pointer or reference to vtkHDF object/struct
 *
 * @memberof vtkHDF
 */
void
vtk_HDF_write_type_attribute(const char* type_name,
                             const hid_t group_id,
                             vtkHDF* vtk_hdf);

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
void
vtk_HDF_read_dataset(const char* dataset_name,
                     void** data,
                     const hid_t dtype_id,
                     const hid_t group_id,
                     const int rank,
                     hsize_t dims[],
                     hsize_t max_dims[],
                     vtkHDF* vtk_hdf);

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
void
vtk_HDF_write_dataset(const char* dataset_name,
                      const void* data,
                      const hid_t dtype_id,
                      const hid_t group_id,
                      const int rank,
                      const hsize_t dims[],
                      vtkHDF* vtk_hdf);

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
void
vtk_HDF_write_chunked_dataset(const char* dataset_name,
                              const void* data,
                              const hid_t dtype_id,
                              const hid_t group_id,
                              const int rank,
                              const hsize_t dims[],
                              const hsize_t max_dims[],
                              const hsize_t chunk_dims[],
                              vtkHDF* vtk_hdf);

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

void
vtk_HDF_append_chunked_dataset(const char* dataset_name,
                               const void* data,
                               const hid_t dtype_id,
                               const hid_t group_id,
                               const int rank,
                               const hsize_t dims[],
                               vtkHDF* vtk_hdf);

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
void
vtk_HDF_write_compressed_dataset(const char* dataset_name,
                                 const void* data,
                                 const hid_t dtype_id,
                                 const hid_t group_id,
                                 const int rank,
                                 const hsize_t dims[],
                                 const hsize_t max_dims[],
                                 const hsize_t chunk_dims[],
                                 unsigned int compression_level,
                                 vtkHDF* vtk_hdf);

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
void
vtk_HDF_collective_write_dataset(const char* dataset_name,
                                 const void* data,
                                 const hid_t dtype_id,
                                 const hid_t group_id,
                                 const int rank,
                                 const hsize_t dims[],
                                 const hsize_t local_size[],
                                 const hsize_t local_offset[],
                                 vtkHDF* vtk_hdf);
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
void
vtk_HDF_collective_write_chunked_dataset(const char* dataset_name,
                                         const void* data,
                                         const hid_t dtype_id,
                                         const hid_t group_id,
                                         const int rank,
                                         const hsize_t dims[],
                                         const hsize_t max_dims[],
                                         const hsize_t chunk_dims[],
                                         const hsize_t local_size[],
                                         const hsize_t local_offset[],
                                         vtkHDF* vtk_hdf);
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

void
vtk_HDF_collective_write_compressed_dataset(const char* dataset_name,
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
                                            vtkHDF* vtk_hdf);
} // namespace C
} // namespace io
} // namespace ELFF