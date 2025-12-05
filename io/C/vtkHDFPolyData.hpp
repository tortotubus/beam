#pragma once

#include "vtkHDF.hpp"
#include "vtkPolyData.hpp"

namespace beam {
namespace io {
namespace C {

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

void
vtk_HDF_polydata_close(vtkHDFPolyData* vtk_hdf_pd);

vtkHDFPolyData
vtk_HDF_polydata_init_static(const char* fname,
                             bool overwrite,
                             vtkPolyData* vtk_pd);
vtkHDFPolyData
vtk_HDF_polydata_init_transient(const char* fname, vtkPolyData* vtk_pd);

vtkHDFPolyData
vtk_HDF_polydata_append_transient(const char* fname,
                                  vtkPolyData* vtk_pd,
                                  float time);

} // namespace C
} // namespace io
} // namespace beam