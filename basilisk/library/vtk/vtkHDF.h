/** @defgroup basilisk Basilisk library
 * 
 *  @{
 */

#include <hdf5.h>

/**
 * @brief Struct for creating a vtkHDF file.
 *
 * This is intended to be "subclassed" for various vtkHDF files; very little is
 * neccesarily in common between different grid types in the VTKHDF schemas,
 * other than the requirement that we put our stuff into a VTKHDF group, which
 * we create here.
 *
 */
typedef struct
{
  hid_t file_id;       /**< File id */
  hid_t fapl;          /**< File access property list  */
  hid_t grp_vtkhdf_id; /**< VTKHDF group id */
} vtkHDF;

/**
 * @brief Close any still open objects
 *
 * @memberof vtkHDF
 */
void
vtk_HDF_close(vtkHDF* vtk_hdf)
{
  if (vtk_hdf->grp_vtkhdf_id >= 0)
    H5Gclose(vtk_hdf->grp_vtkhdf_id);
  if (vtk_hdf->fapl >= 0)
    H5Pclose(vtk_hdf->fapl);
  if (vtk_hdf->file_id >= 0)
    H5Fclose(vtk_hdf->file_id);
}

/**
 * @brief Function called on any error which immediately closes the file to
 * avoid corruptions
 *
 * @memberof vtkHDF
 */
void
vtk_HDF_error(vtkHDF* vtk_hdf)
{
  vtk_HDF_close(vtk_hdf);
}

/**
 * @brief Initialize the vtkHDF struct
 *
 * @memberof vtkHDF
 */
vtkHDF
vtk_HDF_init(const char* fname)
{
  /* File and group identifiers */
  vtkHDF vtk_hdf = { .file_id = H5I_INVALID_HID,
                     .fapl = H5I_INVALID_HID,
                     .grp_vtkhdf_id = H5I_INVALID_HID };

#if _MPI
  // Create a file access property list
  vtk_hdf.fapl = H5Pcreate(H5P_FILE_ACCESS);
  if (vtk_hdf.fapl < 0)
    vtk_HDF_error(&vtk_hdf);

#if MPI_SINGLE_FILE
  // If we use MPI_SINGLE_FILE, we set all procs to have access and write into
  // the same file
  if (H5Pset_fapl_mpio(vtk_hdf.fapl, MPI_COMM_WORLD, MPI_INFO_NULL) < 0)
    vtk_HDF_error(&vtk_hdf);
#else
  // Otherwise, each proccess gets its own file
  if (H5Pset_fapl_mpio(vtk_hdf.fapl, MPI_COMM_SELF, MPI_INFO_NULL) < 0)
    vtk_HDF_error(&vtk_hdf);
#endif

  // Create the acutal file on disk
  vtk_hdf.file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, vtk_hdf.fapl);
  if (vtk_hdf.file_id < 0)
    vtk_HDF_error(&vtk_hdf);

  // Create the VTKHDF group
  vtk_hdf.grp_vtkhdf_id = H5Gcreate2(
    vtk_hdf.file_id, "/VTKHDF", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (vtk_hdf.grp_vtkhdf_id < 0)
    vtk_HDF_error(&vtk_hdf);
#else

  // Create the actual file on disk
  vtk_hdf.file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (vtk_hdf.file_id < 0)
    vtk_HDF_error(&vtk_hdf);

  // Create the VTKHDF group
  vtk_hdf.grp_vtkhdf_id = H5Gcreate2(
    vtk_hdf.file_id, "/VTKHDF", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (vtk_hdf.grp_vtkhdf_id < 0)
    vtk_HDF_error(&vtk_hdf);
#endif

  return vtk_hdf;
}


/** @} */