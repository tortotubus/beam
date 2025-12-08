#pragma once

#include "vtkHDF.hpp"
#include "vtkPolyData.hpp"

namespace beam {
namespace io {
namespace C {

/**
 * @struct vtkHDFPolyData
 *
 * @brief Derived class for writing vtkHDF Polydata files, designed to consume
 * @ref vtkPolyData
 */
typedef struct
{
  /// @privatesection
  vtkHDF vtk_hdf;         /**< Parent @ref vtkHDF */
  hid_t grp_vertices_id;  /**< Group id of /VTKHDF/Vertices  */
  hid_t grp_lines_id;     /**< Group id of /VTKHDF/Lines  */
  hid_t grp_strips_id;    /**< Group id of /VTKHDF/Strips  */
  hid_t grp_polygons_id;  /**< Group id of /VTKHDF/Polygons  */
  hid_t grp_celldata_id;  /**< Group id of /VTKHDF/CellData  */
  hid_t grp_pointdata_id; /**< Group id of /VTKHDF/PointData  */
  hid_t grp_steps_id;     /**< Group id of /VTKHDF/Steps  */
  hid_t
    grp_celldata_offsets_id; /**< Group id of /VTKHDF/Steps/CellDataOffsets  */
  hid_t
    grp_pointdata_offsets_id; /**< Group id of /VTKHDF/Steps/PointDataOffsets */
} vtkHDFPolyData;

/**
 * @brief Close a vtkHDF Polydata file
 *
 * @memberof vtkHDFPolyData
 */
void
vtk_HDF_polydata_close(vtkHDFPolyData* vtk_hdf_pd);

/**
 * @brief Create a new time-independent vtkHDF Polydata file
 *
 * @memberof vtkHDFPolyData
 */
vtkHDFPolyData
vtk_HDF_polydata_init_static(const char* fname,
                             bool overwrite,
                             vtkPolyData* vtk_pd);

/**
 * @brief Create a new time-dependent vtkHDF Polydata file
 *
 * @memberof vtkHDFPolyData
 */
vtkHDFPolyData
vtk_HDF_polydata_init_transient(const char* fname,
                                bool overwrite,
                                vtkPolyData* vtk_pd,
                                float time);

/**
 * @brief Append an existing time-dependent vtkHDF Polydata file
 *
 * @memberof vtkHDFPolyData
 */
vtkHDFPolyData
vtk_HDF_polydata_append_transient(const char* fname,
                                  vtkPolyData* vtk_pd,
                                  float time);

} // namespace C
} // namespace io
} // namespace beam