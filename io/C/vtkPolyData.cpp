#include "vtkPolyData.hpp"

#include <cstdlib>

namespace beam {
namespace io {
namespace C {

bool
vtk_polydata_points_is_sealed(vtkPolyData* pd)
{
  return pd->points_state == SEALED;
}

bool
vtk_polydata_connectivity_is_sealed(vtkPolyData* pd)
{
  return pd->connectivity_state == SEALED;
}

bool
vtk_polydata_point_exists(vtkPolyData* pd, int64_t point)
{
  return 0 <= point && point < pd->n_points;
}

void
vtk_polydata_on_add_points(vtkPolyData* pd)
{
  if (vtk_polydata_points_is_sealed(pd)) {
    abort();
  }
}

void
vtk_polydata_on_add_connectivity(vtkPolyData* pd)
{
  pd->points_state = SEALED;

  if (vtk_polydata_connectivity_is_sealed(pd)) {
    abort();
  }
}

size_t
vtk_polydata_number_of_points(vtkPolyData* pd)
{
  return pd->n_points;
}

size_t
vtk_polydata_number_of_vertices(vtkPolyData* pd)
{
  return pd->n_vertices_offsets - 1;
}

size_t
vtk_polydata_number_of_lines(vtkPolyData* pd)
{
  return pd->n_lines_offsets - 1;
}

size_t
vtk_polydata_number_of_strips(vtkPolyData* pd)
{
  return pd->n_strips_offsets - 1;
}

size_t
vtk_polydata_number_of_polygons(vtkPolyData* pd)
{
  return pd->n_polygons_offsets - 1;
}

void
vtk_polydata_free_points(vtkPolyData* pd)
{
  if (pd->points != nullptr) {
    free(pd->points);
    pd->points = nullptr;
    pd->n_points = 0;
    pd->m_points = 0;
  }
}

void
vtk_polydata_malloc_points(vtkPolyData* pd, size_t n)
{
  vtk_polydata_free_points(pd);

  pd->m_points = n;
  pd->n_points = 0;
  pd->points = (float*)malloc(sizeof(float) * pd->m_points * 3);
}

void
vtk_polydata_free_vertices(vtkPolyData* pd)
{
  if (pd->vertices_connectivity != nullptr) {
    free(pd->vertices_connectivity);
    pd->vertices_connectivity = nullptr;
    pd->n_vertices_connectivity = 0;
    pd->m_vertices_connectivity = 0;

    free(pd->vertices_offsets);
    pd->vertices_offsets = nullptr;
    pd->n_vertices_offsets = 0;
    pd->m_vertices_offsets = 0;
  }
}

void
vtk_polydata_malloc_vertices(vtkPolyData* pd, size_t n)
{
  vtk_polydata_free_vertices(pd);

  pd->m_vertices_connectivity = n;
  pd->n_vertices_connectivity = 0;
  pd->vertices_connectivity =
    (int64_t*)malloc(sizeof(int64_t) * pd->m_vertices_connectivity);

  pd->m_vertices_offsets = n + 1;
  pd->n_vertices_offsets = 1;
  pd->vertices_offsets =
    (int64_t*)malloc(sizeof(int64_t) * pd->m_vertices_offsets);
  pd->vertices_offsets[0] = 0;
}

void
vtk_polydata_free_lines(vtkPolyData* pd)
{
  if (pd->lines_connectivity != nullptr) {
    free(pd->lines_connectivity);
    pd->lines_connectivity = nullptr;
    pd->m_lines_connectivity = 0;
    pd->n_lines_connectivity = 0;

    free(pd->lines_offsets);
    pd->lines_offsets = nullptr;
    pd->m_lines_offsets = 0;
    pd->n_lines_offsets = 0;
  }
}

void
vtk_polydata_malloc_lines(vtkPolyData* pd, size_t n)
{
  vtk_polydata_free_lines(pd);

  pd->m_lines_connectivity = n;
  pd->n_lines_connectivity = 0;
  pd->lines_connectivity =
    (int64_t*)malloc(sizeof(int64_t) * pd->m_lines_connectivity);

  pd->m_lines_offsets = n + 1;
  pd->n_lines_offsets = 1;
  pd->lines_offsets = (int64_t*)malloc(sizeof(int64_t) * pd->m_lines_offsets);
  pd->lines_offsets[0] = 0;
}

void
vtk_polydata_free_strips(vtkPolyData* pd)
{
  if (pd->strips_connectivity != nullptr) {
    free(pd->strips_connectivity);
    pd->strips_connectivity = nullptr;
    pd->m_strips_connectivity = 0;
    pd->n_strips_connectivity = 0;

    free(pd->strips_offsets);
    pd->strips_offsets = nullptr;
    pd->m_strips_offsets = 0;
    pd->n_strips_offsets = 0;
  }
}

void
vtk_polydata_malloc_strips(vtkPolyData* pd, size_t n)
{
  vtk_polydata_free_strips(pd);

  pd->m_strips_connectivity = n;
  pd->n_strips_connectivity = 0;
  pd->strips_connectivity =
    (int64_t*)malloc(sizeof(int64_t) * pd->m_strips_connectivity);

  pd->m_strips_offsets = n + 1;
  pd->n_strips_offsets = 1;
  pd->strips_offsets = (int64_t*)malloc(sizeof(int64_t) * pd->m_strips_offsets);
  pd->strips_offsets[0] = 0;
}

void
vtk_polydata_free_polygons(vtkPolyData* pd)
{
  if (pd->polygons_connectivity != nullptr) {
    free(pd->polygons_connectivity);
    pd->polygons_connectivity = nullptr;
    pd->m_polygons_connectivity = 0;
    pd->n_polygons_connectivity = 0;

    free(pd->polygons_offsets);
    pd->polygons_offsets = nullptr;
    pd->m_polygons_offsets = 0;
    pd->n_polygons_offsets = 0;
  }
}

void
vtk_polydata_malloc_polygons(vtkPolyData* pd, size_t n)
{
  vtk_polydata_free_polygons(pd);

  pd->m_polygons_connectivity = n;
  pd->n_polygons_connectivity = 0;
  pd->polygons_connectivity =
    (int64_t*)malloc(sizeof(int64_t) * pd->m_polygons_connectivity);

  pd->m_polygons_offsets = n + 1;
  pd->n_polygons_offsets = 1;
  pd->polygons_offsets =
    (int64_t*)malloc(sizeof(int64_t) * pd->m_polygons_offsets);
  pd->polygons_offsets[0] = 0;
}

int64_t
vtk_polydata_add_point(vtkPolyData* pd, float x, float y, float z)
{

  vtk_polydata_on_add_points(pd);

  pd->points[(pd->n_points * 3) + 0] = x;
  pd->points[(pd->n_points * 3) + 1] = y;
  pd->points[(pd->n_points * 3) + 2] = z;

  pd->n_points++;

  return pd->n_points - 1;
}

int64_t
vtk_polydata_add_vertex(vtkPolyData* pd, int64_t vertex_point)
{

  vtk_polydata_on_add_connectivity(pd);

  if (!vtk_polydata_point_exists(pd, vertex_point)) {
    abort();
  }

  pd->vertices_connectivity[pd->n_vertices_connectivity] = vertex_point;
  pd->n_vertices_connectivity++;

  pd->vertices_offsets[pd->n_vertices_offsets] = pd->n_vertices_connectivity;
  pd->n_vertices_offsets++;

  return vtk_polydata_number_of_vertices(pd) - 1;
}

int64_t
vtk_polydata_add_line(vtkPolyData* pd,
                      int64_t line_point_1,
                      int64_t line_point_2)
{

  vtk_polydata_on_add_connectivity(pd);

  if (!vtk_polydata_point_exists(pd, line_point_1) ||
      !vtk_polydata_point_exists(pd, line_point_2)) {
    abort();
  }

  pd->lines_connectivity[pd->n_lines_connectivity] = line_point_1;
  pd->n_lines_connectivity++;

  pd->lines_connectivity[pd->n_lines_connectivity] = line_point_2;
  pd->n_lines_connectivity++;

  pd->lines_offsets[pd->n_lines_offsets] = pd->n_lines_connectivity;
  pd->n_lines_offsets++;

  return vtk_polydata_number_of_lines(pd) - 1;
}

vtkPolyData
vtk_polydata_init(size_t n_points,
                  size_t n_vertices,
                  size_t n_lines,
                  size_t n_strips,
                  size_t n_polygons)
{
  vtkPolyData pd = { .points_state = BUILDING,
                     .points = nullptr,
                     .n_points = 0,
                     .m_points = 0,
                     .connectivity_state = BUILDING,
                     .vertices_connectivity = nullptr,
                     .m_vertices_connectivity = 0,
                     .n_vertices_connectivity = 0,
                     .vertices_offsets = nullptr,
                     .n_vertices_offsets = 0,
                     .m_vertices_offsets = 0,
                     .lines_connectivity = nullptr,
                     .n_lines_connectivity = 0,
                     .m_lines_connectivity = 0,
                     .lines_offsets = nullptr,
                     .n_lines_offsets = 0,
                     .m_lines_offsets = 0,
                     .strips_connectivity = nullptr,
                     .n_strips_connectivity = 0,
                     .m_strips_connectivity = 0,
                     .strips_offsets = nullptr,
                     .n_strips_offsets = 0,
                     .m_strips_offsets = 0,
                     .polygons_connectivity = nullptr,
                     .n_polygons_connectivity = 0,
                     .m_polygons_connectivity = 0,
                     .polygons_offsets = nullptr,
                     .n_polygons_offsets = 0,
                     .m_polygons_offsets = 0,
                     .fields_state = BUILDING };

  vtk_polydata_malloc_points(&pd, n_points);
  vtk_polydata_malloc_vertices(&pd, n_vertices);
  vtk_polydata_malloc_lines(&pd, n_lines);
  vtk_polydata_malloc_strips(&pd, n_strips);
  vtk_polydata_malloc_polygons(&pd, n_polygons);

  return pd;
}

void
vtk_polydata_free(vtkPolyData* pd)
{
  vtk_polydata_free_points(pd);
  vtk_polydata_free_vertices(pd);
  vtk_polydata_free_lines(pd);
  vtk_polydata_free_strips(pd);
  vtk_polydata_free_polygons(pd);
}

} // namespace C
} // namespace io
} // namespace beam