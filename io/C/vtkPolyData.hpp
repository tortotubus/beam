#pragma once

#include <cstddef>
#include <cstdint>
#include <cstdlib>

namespace ELFF {
namespace IO {
namespace C {

/**
 * @brief
 *
 * @memberof vtkPolyData
 */
typedef enum
{
  BUILDING = 0,
  SEALED = 1
} vtkPolyDataState;

/**
 * @brief Data object for consumption by @ref vtkHDFPolyData class
 */
typedef struct
{
  vtkPolyDataState points_state;

  float* points;
  size_t n_points;
  size_t m_points;

  vtkPolyDataState connectivity_state;

  int64_t* vertices_connectivity;
  size_t m_vertices_connectivity;
  size_t n_vertices_connectivity;
  int64_t* vertices_offsets;
  size_t n_vertices_offsets;
  size_t m_vertices_offsets;

  int64_t* lines_connectivity;
  size_t n_lines_connectivity;
  size_t m_lines_connectivity;
  int64_t* lines_offsets;
  size_t n_lines_offsets;
  size_t m_lines_offsets;

  int64_t* strips_connectivity;
  size_t n_strips_connectivity;
  size_t m_strips_connectivity;
  int64_t* strips_offsets;
  size_t n_strips_offsets;
  size_t m_strips_offsets;

  int64_t* polygons_connectivity;
  size_t n_polygons_connectivity;
  size_t m_polygons_connectivity;
  int64_t* polygons_offsets;
  size_t n_polygons_offsets;
  size_t m_polygons_offsets;

  vtkPolyDataState fields_state;

  size_t n_pointdata;
  size_t m_pointdata;
  char** pointdata_names;
  size_t* pointdata_ncomp;
  double** pointdata_data;

  size_t n_celldata;
  size_t m_celldata;
  char** celldata_names;
  size_t* celldata_ncomp;
  double** celldata_data;
} vtkPolyData;

/* Function declarations */

/**
 * @memberof vtkPolyData
 */
vtkPolyData
vtk_polydata_init(size_t n_points,
                  size_t n_vertices,
                  size_t n_lines,
                  size_t n_strips,
                  size_t n_polygons,
                  size_t n_pointdata = 0,
                  size_t n_celldata = 0);

/**
 * @memberof vtkPolyData
 */
void
vtk_polydata_free(vtkPolyData* pd);

/**
 * @memberof vtkPolyData
 */
bool
vtk_polydata_points_is_sealed(vtkPolyData* pd);

/**
 * @memberof vtkPolyData
 */
bool
vtk_polydata_connectivity_is_sealed(vtkPolyData* pd);

/**
 * @memberof vtkPolyData
 */
bool
vtk_polydata_data_is_sealed(vtkPolyData* pd);

/**
 * @memberof vtkPolyData
 */
bool
vtk_polydata_point_exists(vtkPolyData* pd, int64_t point);

/**
 * @memberof vtkPolyData
 */
void
vtk_polydata_on_add_points(vtkPolyData* pd);

/**
 * @memberof vtkPolyData
 */
void
vtk_polydata_on_add_connectivity(vtkPolyData* pd);

/**
 * @memberof vtkPolyData
 */
void
vtk_polydata_on_add_data(vtkPolyData* pd);

/**
 * @memberof vtkPolyData
 */
void
vtk_polydata_validate_name(vtkPolyData* pd, const char* name);

/**
 * @memberof vtkPolyData
 */
size_t
vtk_polydata_number_of_points(vtkPolyData* pd);

/**
 * @memberof vtkPolyData
 */
size_t
vtk_polydata_number_of_vertices(vtkPolyData* pd);

/**
 * @memberof vtkPolyData
 */
size_t
vtk_polydata_number_of_lines(vtkPolyData* pd);

/**
 * @memberof vtkPolyData
 */
size_t
vtk_polydata_number_of_strips(vtkPolyData* pd);

/**
 * @memberof vtkPolyData
 */
size_t
vtk_polydata_number_of_polygons(vtkPolyData* pd);

/**
 * @memberof vtkPolyData
 */
size_t
vtk_polydata_number_of_cells(vtkPolyData* pd);

/**
 * @memberof vtkPolyData
 */
void
vtk_polydata_free_points(vtkPolyData* pd);

/**
 * @memberof vtkPolyData
 */
void
vtk_polydata_malloc_points(vtkPolyData* pd, size_t n);

/**
 * @memberof vtkPolyData
 */
void
vtk_polydata_free_vertices(vtkPolyData* pd);

/**
 * @memberof vtkPolyData
 */
void
vtk_polydata_malloc_vertices(vtkPolyData* pd, size_t n);

/**
 * @memberof vtkPolyData
 */
void
vtk_polydata_free_lines(vtkPolyData* pd);

/**
 * @memberof vtkPolyData
 */
void
vtk_polydata_malloc_lines(vtkPolyData* pd, size_t n);

/**
 * @memberof vtkPolyData
 */
void
vtk_polydata_free_strips(vtkPolyData* pd);

/**
 * @memberof vtkPolyData
 */
void
vtk_polydata_malloc_strips(vtkPolyData* pd, size_t n);

/**
 * @memberof vtkPolyData
 */
void
vtk_polydata_free_polygons(vtkPolyData* pd);

/**
 * @memberof vtkPolyData
 */
void
vtk_polydata_malloc_polygons(vtkPolyData* pd, size_t n);

/**
 * @memberof vtkPolyData
 */
int64_t
vtk_polydata_add_point(vtkPolyData* pd, float x, float y, float z);

/**
 * @memberof vtkPolyData
 */
int64_t
vtk_polydata_add_vertex(vtkPolyData* pd, int64_t vertex_point);

/**
 * @memberof vtkPolyData
 */
int64_t
vtk_polydata_add_line(vtkPolyData* pd,
                      int64_t line_point_1,
                      int64_t line_point_2);

/**
 * @memberof vtkPolyData
 *
 * @brief Append source polydata into destination polydata
 *
 * Point-data schema must match exactly (same number of fields, same names,
 * same component counts).
 */
void
vtk_polydata_append(vtkPolyData* dst, const vtkPolyData* src);

/**
 * @memberof vtkPolyData
 */
int64_t
vtk_polydata_add_pointdata_scalar(vtkPolyData* pd, const char* name);

/**
 * @memberof vtkPolyData
 */
int64_t
vtk_polydata_add_pointdata_vector(vtkPolyData* pd,
                                  const char* name,
                                  size_t ncomp);

/**
 * @memberof vtkPolyData
 */
double*
vtk_polydata_get_pointdata(vtkPolyData* pd, int64_t field);

/**
 * @memberof vtkPolyData
 */
double*
vtk_polydata_get_pointdata_data(vtkPolyData* pd, int64_t field);

/**
 * @memberof vtkPolyData
 */
int64_t
vtk_polydata_malloc_pointdata_scalar(vtkPolyData* pd);

/**
 * @memberof vtkPolyData
 */
int64_t
vtk_polydata_malloc_pointdata_vector(vtkPolyData* pd, size_t ncomp);

/**
 * @memberof vtkPolyData
 */
void
vtk_polydata_free_pointdata_field(vtkPolyData* pd, int64_t field);

/**
 * @memberof vtkPolyData
 */
void
vtk_polydata_malloc_pointdata(vtkPolyData* pd, size_t n);

/**
 * @memberof vtkPolyData
 */
void
vtk_polydata_free_pointdata(vtkPolyData* pd);

/**
 * @memberof vtkPolyData
 */
int64_t
vtk_polydata_add_celldata_scalar(vtkPolyData* pd, const char* name);

/**
 * @memberof vtkPolyData
 */
int64_t
vtk_polydata_add_celldata_vector(vtkPolyData* pd,
                                 const char* name,
                                 size_t ncomp);

/**
 * @memberof vtkPolyData
 */
double*
vtk_polydata_get_celldata(vtkPolyData* pd, int64_t field);

/**
 * @memberof vtkPolyData
 */
double*
vtk_polydata_get_celldata_data(vtkPolyData* pd, int64_t field);

/**
 * @memberof vtkPolyData
 */
int64_t
vtk_polydata_malloc_celldata_scalar(vtkPolyData* pd);

/**
 * @memberof vtkPolyData
 */
int64_t
vtk_polydata_malloc_celldata_vector(vtkPolyData* pd, size_t ncomp);

/**
 * @memberof vtkPolyData
 */
void
vtk_polydata_free_celldata_field(vtkPolyData* pd, int64_t field);

/**
 * @memberof vtkPolyData
 */
void
vtk_polydata_malloc_celldata(vtkPolyData* pd, size_t n);

/**
 * @memberof vtkPolyData
 */
void
vtk_polydata_free_celldata(vtkPolyData* pd);

} // namespace C
} // namespace IO
} // namespace ELFF
