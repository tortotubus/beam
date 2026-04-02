#include "elff/io/C/vtkPolyData.hpp"

#include <cassert>
#include <cstdlib>
#include <cstring>

namespace ELFF {
namespace IO {
namespace C {

static void
vtk_polydata_require_point_capacity(vtkPolyData* pd, size_t required_points)
{
  if (required_points <= pd->m_points) {
    return;
  }

  pd->points = static_cast<float*>(realloc(pd->points, sizeof(float) * required_points * 3));
  assert(pd->points != nullptr);
  pd->m_points = required_points;
}

static void
vtk_polydata_require_connectivity_capacity(int64_t** connectivity,
                                           size_t* m_connectivity,
                                           size_t required_connectivity,
                                           int64_t** offsets,
                                           size_t* m_offsets,
                                           size_t required_offsets)
{
  if (required_connectivity > *m_connectivity) {
    *connectivity = static_cast<int64_t*>(
      realloc(*connectivity, sizeof(int64_t) * required_connectivity));
    assert(*connectivity != nullptr);
    *m_connectivity = required_connectivity;
  }

  if (required_offsets > *m_offsets) {
    *offsets =
      static_cast<int64_t*>(realloc(*offsets, sizeof(int64_t) * required_offsets));
    assert(*offsets != nullptr);
    *m_offsets = required_offsets;
  }
}

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
vtk_polydata_data_is_sealed(vtkPolyData* pd)
{
  return pd->fields_state == SEALED;
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

void
vtk_polydata_on_add_data(vtkPolyData* pd)
{
  pd->connectivity_state = SEALED;

  if (vtk_polydata_data_is_sealed(pd)) {
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

size_t
vtk_polydata_number_of_cells(vtkPolyData* pd)
{
  return vtk_polydata_number_of_vertices(pd) +
         vtk_polydata_number_of_lines(pd) +
         vtk_polydata_number_of_polygons(pd) +
         vtk_polydata_number_of_strips(pd);
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
  pd->points = static_cast<float*>(malloc(sizeof(float) * pd->m_points * 3));
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
  pd->vertices_connectivity = static_cast<int64_t*>(
    malloc(sizeof(int64_t) * pd->m_vertices_connectivity));

  pd->m_vertices_offsets = n + 1;
  pd->n_vertices_offsets = 1;
  pd->vertices_offsets =
    static_cast<int64_t*>(malloc(sizeof(int64_t) * pd->m_vertices_offsets));
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
    static_cast<int64_t*>(malloc(sizeof(int64_t) * pd->m_lines_connectivity));

  pd->m_lines_offsets = n + 1;
  pd->n_lines_offsets = 1;
  pd->lines_offsets =
    static_cast<int64_t*>(malloc(sizeof(int64_t) * pd->m_lines_offsets));
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
  pd->strips_connectivity = static_cast<int64_t*>(
    malloc(sizeof(int64_t) * pd->m_strips_connectivity));

  pd->m_strips_offsets = n + 1;
  pd->n_strips_offsets = 1;
  pd->strips_offsets =
    static_cast<int64_t*>(malloc(sizeof(int64_t) * pd->m_strips_offsets));
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
  pd->polygons_connectivity = static_cast<int64_t*>(
    malloc(sizeof(int64_t) * pd->m_polygons_connectivity));

  pd->m_polygons_offsets = n + 1;
  pd->n_polygons_offsets = 1;
  pd->polygons_offsets =
    static_cast<int64_t*>(malloc(sizeof(int64_t) * pd->m_polygons_offsets));
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

void
vtk_polydata_append(vtkPolyData* dst, const vtkPolyData* src)
{
  assert(dst != nullptr);
  assert(src != nullptr);
  assert(dst != src);

  assert(dst->n_pointdata == src->n_pointdata);
  for (size_t i = 0; i < dst->n_pointdata; ++i) {
    assert(dst->pointdata_ncomp[i] == src->pointdata_ncomp[i]);
    assert(strcmp(dst->pointdata_names[i], src->pointdata_names[i]) == 0);
  }
  assert(dst->n_celldata == src->n_celldata);
  for (size_t i = 0; i < dst->n_celldata; ++i) {
    assert(dst->celldata_ncomp[i] == src->celldata_ncomp[i]);
    assert(strcmp(dst->celldata_names[i], src->celldata_names[i]) == 0);
  }

  const size_t old_n_points = dst->n_points;
  const size_t src_n_points = src->n_points;
  const size_t old_n_cells = (dst->n_vertices_offsets - 1) + (dst->n_lines_offsets - 1) +
                             (dst->n_polygons_offsets - 1) + (dst->n_strips_offsets - 1);
  const size_t src_n_cells = (src->n_vertices_offsets - 1) + (src->n_lines_offsets - 1) +
                             (src->n_polygons_offsets - 1) + (src->n_strips_offsets - 1);
  const int64_t point_offset = static_cast<int64_t>(old_n_points);

  // Points
  vtk_polydata_require_point_capacity(dst, old_n_points + src_n_points);
  if (src_n_points > 0) {
    memcpy(dst->points + (old_n_points * 3), src->points, sizeof(float) * src_n_points * 3);
  }
  dst->n_points = old_n_points + src_n_points;

  auto append_connectivity = [point_offset](int64_t* dst_connectivity,
                                            size_t* dst_n_connectivity,
                                            int64_t* dst_offsets,
                                            size_t* dst_n_offsets,
                                            const int64_t* src_connectivity,
                                            size_t src_n_connectivity,
                                            const int64_t* src_offsets,
                                            size_t src_n_offsets) {
    const size_t old_n_connectivity = *dst_n_connectivity;

    for (size_t i = 0; i < src_n_connectivity; ++i) {
      dst_connectivity[*dst_n_connectivity + i] = src_connectivity[i] + point_offset;
    }
    *dst_n_connectivity += src_n_connectivity;

    for (size_t i = 1; i < src_n_offsets; ++i) {
      dst_offsets[*dst_n_offsets] = static_cast<int64_t>(old_n_connectivity) + src_offsets[i];
      *dst_n_offsets += 1;
    }
  };

  // Vertices
  vtk_polydata_require_connectivity_capacity(&dst->vertices_connectivity,
                                             &dst->m_vertices_connectivity,
                                             dst->n_vertices_connectivity +
                                               src->n_vertices_connectivity,
                                             &dst->vertices_offsets,
                                             &dst->m_vertices_offsets,
                                             dst->n_vertices_offsets + src->n_vertices_offsets - 1);
  append_connectivity(dst->vertices_connectivity,
                      &dst->n_vertices_connectivity,
                      dst->vertices_offsets,
                      &dst->n_vertices_offsets,
                      src->vertices_connectivity,
                      src->n_vertices_connectivity,
                      src->vertices_offsets,
                      src->n_vertices_offsets);

  // Lines
  vtk_polydata_require_connectivity_capacity(
    &dst->lines_connectivity,
    &dst->m_lines_connectivity,
    dst->n_lines_connectivity + src->n_lines_connectivity,
    &dst->lines_offsets,
    &dst->m_lines_offsets,
    dst->n_lines_offsets + src->n_lines_offsets - 1);
  append_connectivity(dst->lines_connectivity,
                      &dst->n_lines_connectivity,
                      dst->lines_offsets,
                      &dst->n_lines_offsets,
                      src->lines_connectivity,
                      src->n_lines_connectivity,
                      src->lines_offsets,
                      src->n_lines_offsets);

  // Strips
  vtk_polydata_require_connectivity_capacity(
    &dst->strips_connectivity,
    &dst->m_strips_connectivity,
    dst->n_strips_connectivity + src->n_strips_connectivity,
    &dst->strips_offsets,
    &dst->m_strips_offsets,
    dst->n_strips_offsets + src->n_strips_offsets - 1);
  append_connectivity(dst->strips_connectivity,
                      &dst->n_strips_connectivity,
                      dst->strips_offsets,
                      &dst->n_strips_offsets,
                      src->strips_connectivity,
                      src->n_strips_connectivity,
                      src->strips_offsets,
                      src->n_strips_offsets);

  // Polygons
  vtk_polydata_require_connectivity_capacity(
    &dst->polygons_connectivity,
    &dst->m_polygons_connectivity,
    dst->n_polygons_connectivity + src->n_polygons_connectivity,
    &dst->polygons_offsets,
    &dst->m_polygons_offsets,
    dst->n_polygons_offsets + src->n_polygons_offsets - 1);
  append_connectivity(dst->polygons_connectivity,
                      &dst->n_polygons_connectivity,
                      dst->polygons_offsets,
                      &dst->n_polygons_offsets,
                      src->polygons_connectivity,
                      src->n_polygons_connectivity,
                      src->polygons_offsets,
                      src->n_polygons_offsets);

  // Point-data
  for (size_t i = 0; i < dst->n_pointdata; ++i) {
    const size_t ncomp = dst->pointdata_ncomp[i];
    const size_t old_tuples = old_n_points * ncomp;
    const size_t src_tuples = src_n_points * ncomp;
    const size_t new_tuples = old_tuples + src_tuples;

    dst->pointdata_data[i] =
      static_cast<double*>(realloc(dst->pointdata_data[i], sizeof(double) * new_tuples));
    assert(dst->pointdata_data[i] != nullptr);

    if (src_tuples > 0) {
      memcpy(dst->pointdata_data[i] + old_tuples,
             src->pointdata_data[i],
             sizeof(double) * src_tuples);
    }
  }
  for (size_t i = 0; i < dst->n_celldata; ++i) {
    const size_t ncomp = dst->celldata_ncomp[i];
    const size_t old_tuples = old_n_cells * ncomp;
    const size_t src_tuples = src_n_cells * ncomp;
    const size_t new_tuples = old_tuples + src_tuples;

    dst->celldata_data[i] =
      static_cast<double*>(realloc(dst->celldata_data[i], sizeof(double) * new_tuples));
    assert(dst->celldata_data[i] != nullptr);

    if (src_tuples > 0) {
      memcpy(dst->celldata_data[i] + old_tuples,
             src->celldata_data[i],
             sizeof(double) * src_tuples);
    }
  }

  if (vtk_polydata_number_of_cells(dst) > 0) {
    dst->points_state = SEALED;
  }
  if (dst->n_pointdata > 0 || dst->n_celldata > 0) {
    dst->connectivity_state = SEALED;
  }
}

int64_t
vtk_polydata_add_pointdata_scalar(vtkPolyData* pd, const char* name)
{
  vtk_polydata_on_add_data(pd);

  int64_t id = vtk_polydata_malloc_pointdata_scalar(pd);

  vtk_polydata_validate_name(pd, name);

  pd->pointdata_names[id] = static_cast<char*>(malloc(sizeof(char) * strlen(name) + 1));
  strcpy(pd->pointdata_names[id], name);

  return id;
}

double*
vtk_polydata_get_pointdata(vtkPolyData* pd, int64_t field)
{
  return vtk_polydata_get_pointdata_data(pd, field);
}

double*
vtk_polydata_get_pointdata_data(vtkPolyData* pd, int64_t field)
{
  size_t id = static_cast<size_t>(field);

  assert(id < pd->n_pointdata);
  assert(pd->pointdata_data[id] != nullptr);

  return pd->pointdata_data[id];
}

int64_t
vtk_polydata_add_pointdata_vector(vtkPolyData* pd,
                                  const char* name,
                                  size_t ncomp)
{
  vtk_polydata_on_add_data(pd);

  int64_t id = vtk_polydata_malloc_pointdata_vector(pd, ncomp);

  vtk_polydata_validate_name(pd, name);

  pd->pointdata_names[id] = static_cast<char*>(malloc(sizeof(char) * strlen(name) + 1));
  strcpy(pd->pointdata_names[id], name);

  return id;
}

void
vtk_polydata_malloc_pointdata(vtkPolyData* pd, size_t n)
{
  vtk_polydata_free_pointdata(pd);

  pd->m_pointdata = n;
  pd->n_pointdata = 0;

  pd->pointdata_data = static_cast<double**>(calloc(pd->m_pointdata, sizeof(double*)));
  pd->pointdata_names = static_cast<char**>(calloc(pd->m_pointdata, sizeof(char*)));
  pd->pointdata_ncomp = static_cast<size_t*>(calloc(pd->m_pointdata, sizeof(size_t)));
}

void
vtk_polydata_free_pointdata(vtkPolyData* pd)
{
  for (size_t i = 0; i < pd->n_pointdata; i++) {
    vtk_polydata_free_pointdata_field(pd, static_cast<int64_t>(i));
  }

  free(pd->pointdata_data);
  pd->pointdata_data = nullptr;
  free(pd->pointdata_names);
  pd->pointdata_names = nullptr;
  free(pd->pointdata_ncomp);
  pd->pointdata_ncomp = nullptr;

  pd->m_pointdata = 0;
  pd->n_pointdata = 0;
}

int64_t
vtk_polydata_malloc_pointdata_scalar(vtkPolyData* pd)
{
  assert(pd->n_pointdata < pd->m_pointdata);

  size_t id = pd->n_pointdata;
  pd->n_pointdata++;

  pd->pointdata_ncomp[id] = 1;
  pd->pointdata_data[id] =
    static_cast<double*>(calloc(pd->n_points, sizeof(double)));

  return static_cast<int64_t>(id);
}

int64_t
vtk_polydata_malloc_pointdata_vector(vtkPolyData* pd, size_t ncomp)
{
  assert(pd->n_pointdata < pd->m_pointdata);
  assert(ncomp > 0);

  size_t id = pd->n_pointdata;
  pd->n_pointdata++;

  pd->pointdata_ncomp[id] = ncomp;
  pd->pointdata_data[id] = static_cast<double*>(
    calloc(pd->n_points * pd->pointdata_ncomp[id], sizeof(double)));

  return static_cast<int64_t>(id);
}

void
vtk_polydata_free_pointdata_field(vtkPolyData* pd, int64_t field)
{
  size_t id = static_cast<size_t>(field);
  assert(id < pd->n_pointdata);

  pd->pointdata_ncomp[id] = 0;

  free(pd->pointdata_names[id]);
  pd->pointdata_names[id] = nullptr;

  free(pd->pointdata_data[id]);
  pd->pointdata_data[id] = nullptr;
}

void
vtk_polydata_validate_name(vtkPolyData* pd, const char* name)
{
  (void)pd;
  (void)name;

  return;
}

vtkPolyData
vtk_polydata_init(size_t n_points,
                  size_t n_vertices,
                  size_t n_lines,
                  size_t n_strips,
                  size_t n_polygons,
                  size_t n_pointdata,
                  size_t n_celldata)
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
                     .fields_state = BUILDING,
                     .n_pointdata = 0,
                     .m_pointdata = 0,
                     .pointdata_names = nullptr,
                     .pointdata_ncomp = nullptr,
                     .pointdata_data = nullptr,
                     .n_celldata = 0,
                     .m_celldata = 0,
                     .celldata_names = nullptr,
                     .celldata_ncomp = nullptr,
                     .celldata_data = nullptr };

  vtk_polydata_malloc_points(&pd, n_points);
  vtk_polydata_malloc_vertices(&pd, n_vertices);
  vtk_polydata_malloc_lines(&pd, n_lines);
  vtk_polydata_malloc_strips(&pd, n_strips);
  vtk_polydata_malloc_polygons(&pd, n_polygons);
  vtk_polydata_malloc_pointdata(&pd, n_pointdata);
  vtk_polydata_malloc_celldata(&pd, n_celldata);

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
  vtk_polydata_free_pointdata(pd);
  vtk_polydata_free_celldata(pd);
}

int64_t
vtk_polydata_add_celldata_scalar(vtkPolyData* pd, const char* name)
{
  vtk_polydata_on_add_data(pd);

  int64_t id = vtk_polydata_malloc_celldata_scalar(pd);

  vtk_polydata_validate_name(pd, name);

  pd->celldata_names[id] =
    static_cast<char*>(malloc(sizeof(char) * strlen(name) + 1));
  strcpy(pd->celldata_names[id], name);

  return id;
}

double*
vtk_polydata_get_celldata(vtkPolyData* pd, int64_t field)
{
  return vtk_polydata_get_celldata_data(pd, field);
}

double*
vtk_polydata_get_celldata_data(vtkPolyData* pd, int64_t field)
{
  size_t id = static_cast<size_t>(field);

  assert(id < pd->n_celldata);
  assert(pd->celldata_data[id] != nullptr);

  return pd->celldata_data[id];
}

int64_t
vtk_polydata_add_celldata_vector(vtkPolyData* pd,
                                 const char* name,
                                 size_t ncomp)
{
  vtk_polydata_on_add_data(pd);

  int64_t id = vtk_polydata_malloc_celldata_vector(pd, ncomp);

  vtk_polydata_validate_name(pd, name);

  pd->celldata_names[id] =
    static_cast<char*>(malloc(sizeof(char) * strlen(name) + 1));
  strcpy(pd->celldata_names[id], name);

  return id;
}

void
vtk_polydata_malloc_celldata(vtkPolyData* pd, size_t n)
{
  vtk_polydata_free_celldata(pd);

  pd->m_celldata = n;
  pd->n_celldata = 0;

  pd->celldata_data =
    static_cast<double**>(calloc(pd->m_celldata, sizeof(double*)));
  pd->celldata_names =
    static_cast<char**>(calloc(pd->m_celldata, sizeof(char*)));
  pd->celldata_ncomp =
    static_cast<size_t*>(calloc(pd->m_celldata, sizeof(size_t)));
}

void
vtk_polydata_free_celldata(vtkPolyData* pd)
{
  for (size_t i = 0; i < pd->n_celldata; i++) {
    vtk_polydata_free_celldata_field(pd, static_cast<int64_t>(i));
  }

  free(pd->celldata_data);
  pd->celldata_data = nullptr;
  free(pd->celldata_names);
  pd->celldata_names = nullptr;
  free(pd->celldata_ncomp);
  pd->celldata_ncomp = nullptr;

  pd->m_celldata = 0;
  pd->n_celldata = 0;
}

int64_t
vtk_polydata_malloc_celldata_scalar(vtkPolyData* pd)
{
  assert(pd->n_celldata < pd->m_celldata);

  size_t id = pd->n_celldata;
  pd->n_celldata++;

  pd->celldata_ncomp[id] = 1;
  pd->celldata_data[id] = static_cast<double*>(
    calloc(vtk_polydata_number_of_cells(pd), sizeof(double)));

  return static_cast<int64_t>(id);
}

int64_t
vtk_polydata_malloc_celldata_vector(vtkPolyData* pd, size_t ncomp)
{
  assert(pd->n_celldata < pd->m_celldata);
  assert(ncomp > 0);

  size_t id = pd->n_celldata;
  pd->n_celldata++;

  pd->celldata_ncomp[id] = ncomp;
  pd->celldata_data[id] = static_cast<double*>(
    calloc(vtk_polydata_number_of_cells(pd) * pd->celldata_ncomp[id],
           sizeof(double)));

  return static_cast<int64_t>(id);
}

void
vtk_polydata_free_celldata_field(vtkPolyData* pd, int64_t field)
{
  size_t id = static_cast<size_t>(field);
  assert(id < pd->n_celldata);

  pd->celldata_ncomp[id] = 0;

  free(pd->celldata_names[id]);
  pd->celldata_names[id] = nullptr;

  free(pd->celldata_data[id]);
  pd->celldata_data[id] = nullptr;
}

} // namespace C
} // namespace IO
} // namespace ELFF
