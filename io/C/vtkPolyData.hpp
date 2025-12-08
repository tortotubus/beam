#pragma once

#include <cstddef>
#include <cstdint>
#include <cstdlib>

namespace beam {
namespace io {
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
                  size_t n_polygons);

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
} // namespace C
} // namespace io
} // namespace beam