#pragma once

/**
 * @brief
 *
 * @memberof vtkPolyData
 */
typedef enum { BUILDING = 0, SEALED = 1 } vtkPolyDataState;

/**
 * @brief Data object for consumption by @ref vtkHDFPolyData class
 */
typedef struct {
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

  vtkPolyDataState data_state;

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

//
//
// Function declarations
//
//

vtkPolyData vtk_polydata_init (size_t n_points,
                               size_t n_vertices,
                               size_t n_lines,
                               size_t n_strips,
                               size_t n_polygons,
                               size_t n_pointdata);

void vtk_polydata_free (vtkPolyData* pd);

bool vtk_polydata_points_is_sealed (vtkPolyData* pd);
bool vtk_polydata_connectivity_is_sealed (vtkPolyData* pd);
bool vtk_polydata_data_is_sealed (vtkPolyData* pd);

void vtk_polydata_on_add_points (vtkPolyData* pd);
void vtk_polydata_on_add_connectivity (vtkPolyData* pd);
void vtk_polydata_on_add_data (vtkPolyData* pd);

void vtk_polydata_validate_name (vtkPolyData* pd, const char* name);

bool vtk_polydata_point_exists (vtkPolyData* pd, int64_t point);

int64_t vtk_polydata_add_point (vtkPolyData* pd, float x, float y, float z);
size_t vtk_polydata_number_of_points (vtkPolyData* pd);
void vtk_polydata_free_points (vtkPolyData* pd);
void vtk_polydata_malloc_points (vtkPolyData* pd, size_t n);

int64_t vtk_polydata_add_vertex (vtkPolyData* pd, int64_t vertex_point);
size_t vtk_polydata_number_of_vertices (vtkPolyData* pd);
void vtk_polydata_free_vertices (vtkPolyData* pd);
void vtk_polydata_malloc_vertices (vtkPolyData* pd, size_t n);

int64_t vtk_polydata_add_line (vtkPolyData* pd,
                               int64_t line_point_1,
                               int64_t line_point_2);
size_t vtk_polydata_number_of_lines (vtkPolyData* pd);
void vtk_polydata_free_lines (vtkPolyData* pd);
void vtk_polydata_malloc_lines (vtkPolyData* pd, size_t n);

// int64_t vtk_polydata_add_strip (vtkPolyData* pd);
size_t vtk_polydata_number_of_strips (vtkPolyData* pd);
void vtk_polydata_free_strips (vtkPolyData* pd);
void vtk_polydata_malloc_strips (vtkPolyData* pd, size_t n);

// int64_t vtk_polydata_add_polygon (vtkPolyData* pd);
size_t vtk_polydata_number_of_polygons (vtkPolyData* pd);
void vtk_polydata_free_polygons (vtkPolyData* pd);
void vtk_polydata_malloc_polygons (vtkPolyData* pd, size_t n);

int64_t vtk_polydata_add_pointdata_scalar (vtkPolyData* pd, const char* name);
int64_t vtk_polydata_add_pointdata_vector (vtkPolyData* pd,
                                           const char* name,
                                           size_t ncomp);
double* vtk_polydata_get_pointdata (vtkPolyData* pd, int64_t field);

int64_t vtk_polydata_malloc_pointdata_scalar (vtkPolyData* pd);
int64_t vtk_polydata_malloc_pointdata_vector (vtkPolyData* pd, size_t ncomp);
void vtk_polydata_free_pointdata_field (vtkPolyData* pd, int64_t field);
void vtk_polydata_malloc_pointdata (vtkPolyData* pd, size_t n);
void vtk_polydata_free_pointdata (vtkPolyData* pd);

//
//
// Function definitions
//
//

/**
 * @memberof vtkPolyData
 */
vtkPolyData vtk_polydata_init (size_t n_points,
                               size_t n_vertices,
                               size_t n_lines,
                               size_t n_strips,
                               size_t n_polygons,
                               size_t n_pointdata) {
  // vtkPolyData pd = {.points_state = BUILDING,
  //                   .points = NULL,
  //                   .n_points = 0,
  //                   .m_points = 0,
  //                   .connectivity_state = BUILDING,
  //                   .vertices_connectivity = NULL,
  //                   .m_vertices_connectivity = 0,
  //                   .n_vertices_connectivity = 0,
  //                   .vertices_offsets = NULL,
  //                   .n_vertices_offsets = 0,
  //                   .m_vertices_offsets = 0,
  //                   .lines_connectivity = NULL,
  //                   .n_lines_connectivity = 0,
  //                   .m_lines_connectivity = 0,
  //                   .lines_offsets = NULL,
  //                   .n_lines_offsets = 0,
  //                   .m_lines_offsets = 0,
  //                   .strips_connectivity = NULL,
  //                   .n_strips_connectivity = 0,
  //                   .m_strips_connectivity = 0,
  //                   .strips_offsets = NULL,
  //                   .n_strips_offsets = 0,
  //                   .m_strips_offsets = 0,
  //                   .polygons_connectivity = NULL,
  //                   .n_polygons_connectivity = 0,
  //                   .m_polygons_connectivity = 0,
  //                   .polygons_offsets = NULL,
  //                   .n_polygons_offsets = 0,
  //                   .m_polygons_offsets = 0,
  //                   // .data_state = BUILDING,
  //                   .data_state = 0,
  //                   .m_pointdata = 0,
  //                   .n_pointdata = 0,
  //                   .pointdata_names = NULL,
  //                   .pointdata_data = NULL,
  //                   .pointdata_ncomp = NULL};

  vtkPolyData pd = {0};

  vtk_polydata_malloc_points (&pd, n_points);
  vtk_polydata_malloc_vertices (&pd, n_vertices);
  vtk_polydata_malloc_lines (&pd, n_lines);
  vtk_polydata_malloc_strips (&pd, n_strips);
  vtk_polydata_malloc_polygons (&pd, n_polygons);
  vtk_polydata_malloc_pointdata (&pd, n_pointdata);

  return pd;
}

/**
 * @memberof vtkPolyData
 */
void vtk_polydata_free (vtkPolyData* pd) {
  vtk_polydata_free_points (pd);
  vtk_polydata_free_vertices (pd);
  vtk_polydata_free_lines (pd);
  vtk_polydata_free_strips (pd);
  vtk_polydata_free_polygons (pd);
  vtk_polydata_free_pointdata (pd);
}

/**
 * @memberof vtkPolyData
 */
bool vtk_polydata_points_is_sealed (vtkPolyData* pd) {
  return pd->points_state == 1;
}

/**
 * @memberof vtkPolyData
 */
bool vtk_polydata_connectivity_is_sealed (vtkPolyData* pd) {
  return pd->connectivity_state == 1; // SEALED = 1
}

/**
 * @memberof vtkPolyData
 */
bool vtk_polydata_data_is_sealed (vtkPolyData* pd) {
  return pd->data_state == 1; // SEALED = 1
}

/**
 * @memberof vtkPolyData
 */
bool vtk_polydata_point_exists (vtkPolyData* pd, int64_t point) {
  return 0 <= point && point < pd->n_points;
}

/**
 * @memberof vtkPolyData
 */
void vtk_polydata_on_add_points (vtkPolyData* pd)

{
  if (vtk_polydata_points_is_sealed (pd)) {
    abort ();
  }
}

/**
 * @memberof vtkPolyData
 */
void vtk_polydata_on_add_connectivity (vtkPolyData* pd) {
  pd->points_state = SEALED;

  if (vtk_polydata_connectivity_is_sealed (pd)) {
    abort ();
  }
}

/**
 * @memberof vtkPolyData
 */
void vtk_polydata_on_add_data (vtkPolyData* pd) {
  pd->connectivity_state = SEALED;

  if (vtk_polydata_data_is_sealed (pd)) {
    abort ();
  }
}

/**
 * @memberof vtkPolyData
 */
size_t vtk_polydata_number_of_points (vtkPolyData* pd) {
  return pd->n_points;
}

/**
 * @memberof vtkPolyData
 */
size_t vtk_polydata_number_of_vertices (vtkPolyData* pd) {
  return pd->n_vertices_offsets - 1;
}

/**
 * @memberof vtkPolyData
 */
size_t vtk_polydata_number_of_lines (vtkPolyData* pd) {
  return pd->n_lines_offsets - 1;
}

/**
 * @memberof vtkPolyData
 */
size_t vtk_polydata_number_of_strips (vtkPolyData* pd) {
  return pd->n_strips_offsets - 1;
}

/**
 * @memberof vtkPolyData
 */
size_t vtk_polydata_number_of_polygons (vtkPolyData* pd) {
  return pd->n_polygons_offsets - 1;
}

/**
 * @memberof vtkPolyData
 */
void vtk_polydata_free_points (vtkPolyData* pd) {
  if (pd->points != nullptr) {
    free (pd->points);
    pd->points = nullptr;
    pd->n_points = 0;
    pd->m_points = 0;
  }
}

/**
 * @memberof vtkPolyData
 */
void vtk_polydata_malloc_points (vtkPolyData* pd, size_t n) {
  vtk_polydata_free_points (pd);

  pd->m_points = n;
  pd->n_points = 0;
  pd->points = (float*) malloc (sizeof (float) * pd->m_points * 3);
}

/**
 * @memberof vtkPolyData
 */
void vtk_polydata_free_vertices (vtkPolyData* pd) {
  if (pd->vertices_connectivity != nullptr) {
    free (pd->vertices_connectivity);
    pd->vertices_connectivity = nullptr;
    pd->n_vertices_connectivity = 0;
    pd->m_vertices_connectivity = 0;

    free (pd->vertices_offsets);
    pd->vertices_offsets = nullptr;
    pd->n_vertices_offsets = 0;
    pd->m_vertices_offsets = 0;
  }
}

/**
 * @memberof vtkPolyData
 */
void vtk_polydata_malloc_vertices (vtkPolyData* pd, size_t n) {
  vtk_polydata_free_vertices (pd);

  pd->m_vertices_connectivity = n;
  pd->n_vertices_connectivity = 0;
  pd->vertices_connectivity =
    (int64_t*) malloc (sizeof (int64_t) * pd->m_vertices_connectivity);

  pd->m_vertices_offsets = n + 1;
  pd->n_vertices_offsets = 1;
  pd->vertices_offsets =
    (int64_t*) malloc (sizeof (int64_t) * pd->m_vertices_offsets);
  pd->vertices_offsets[0] = 0;
}

/**
 * @memberof vtkPolyData
 */
void vtk_polydata_free_lines (vtkPolyData* pd) {
  if (pd->lines_connectivity != nullptr) {
    free (pd->lines_connectivity);
    pd->lines_connectivity = nullptr;
    pd->m_lines_connectivity = 0;
    pd->n_lines_connectivity = 0;

    free (pd->lines_offsets);
    pd->lines_offsets = nullptr;
    pd->m_lines_offsets = 0;
    pd->n_lines_offsets = 0;
  }
}
/**
 * @memberof vtkPolyData
 */
void vtk_polydata_malloc_lines (vtkPolyData* pd, size_t n) {
  vtk_polydata_free_lines (pd);

  pd->m_lines_connectivity = n;
  pd->n_lines_connectivity = 0;
  pd->lines_connectivity =
    (int64_t*) malloc (sizeof (int64_t) * pd->m_lines_connectivity);

  pd->m_lines_offsets = n + 1;
  pd->n_lines_offsets = 1;
  pd->lines_offsets =
    (int64_t*) malloc (sizeof (int64_t) * pd->m_lines_offsets);
  pd->lines_offsets[0] = 0;
}

/**
 * @memberof vtkPolyData
 */
void vtk_polydata_free_strips (vtkPolyData* pd) {
  if (pd->strips_connectivity != nullptr) {
    free (pd->strips_connectivity);
    pd->strips_connectivity = nullptr;
    pd->m_strips_connectivity = 0;
    pd->n_strips_connectivity = 0;

    free (pd->strips_offsets);
    pd->strips_offsets = nullptr;
    pd->m_strips_offsets = 0;
    pd->n_strips_offsets = 0;
  }
}

/**
 * @memberof vtkPolyData
 */
void vtk_polydata_malloc_strips (vtkPolyData* pd, size_t n) {
  vtk_polydata_free_strips (pd);

  pd->m_strips_connectivity = n;
  pd->n_strips_connectivity = 0;
  pd->strips_connectivity =
    (int64_t*) malloc (sizeof (int64_t) * pd->m_strips_connectivity);

  pd->m_strips_offsets = n + 1;
  pd->n_strips_offsets = 1;
  pd->strips_offsets =
    (int64_t*) malloc (sizeof (int64_t) * pd->m_strips_offsets);
  pd->strips_offsets[0] = 0;
}

/**
 * @memberof vtkPolyData
 */
void vtk_polydata_free_polygons (vtkPolyData* pd) {
  if (pd->polygons_connectivity != nullptr) {
    free (pd->polygons_connectivity);
    pd->polygons_connectivity = nullptr;
    pd->m_polygons_connectivity = 0;
    pd->n_polygons_connectivity = 0;

    free (pd->polygons_offsets);
    pd->polygons_offsets = nullptr;
    pd->m_polygons_offsets = 0;
    pd->n_polygons_offsets = 0;
  }
}

/**
 * @memberof vtkPolyData
 */
void vtk_polydata_malloc_polygons (vtkPolyData* pd, size_t n) {
  vtk_polydata_free_polygons (pd);

  pd->m_polygons_connectivity = n;
  pd->n_polygons_connectivity = 0;
  pd->polygons_connectivity =
    (int64_t*) malloc (sizeof (int64_t) * pd->m_polygons_connectivity);

  pd->m_polygons_offsets = n + 1;
  pd->n_polygons_offsets = 1;
  pd->polygons_offsets =
    (int64_t*) malloc (sizeof (int64_t) * pd->m_polygons_offsets);
  pd->polygons_offsets[0] = 0;
}

/**
 * @memberof vtkPolyData
 */
int64_t vtk_polydata_add_point (vtkPolyData* pd, float x, float y, float z) {

  vtk_polydata_on_add_points (pd);

  pd->points[(pd->n_points * 3) + 0] = x;
  pd->points[(pd->n_points * 3) + 1] = y;
  pd->points[(pd->n_points * 3) + 2] = z;

  pd->n_points++;

  return pd->n_points - 1;
}

/**
 * @memberof vtkPolyData
 */
int64_t vtk_polydata_add_vertex (vtkPolyData* pd, int64_t vertex_point) {

  vtk_polydata_on_add_connectivity (pd);

  if (!vtk_polydata_point_exists (pd, vertex_point)) {
    abort ();
  }

  pd->vertices_connectivity[pd->n_vertices_connectivity] = vertex_point;
  pd->n_vertices_connectivity++;

  pd->vertices_offsets[pd->n_vertices_offsets] = pd->n_vertices_connectivity;
  pd->n_vertices_offsets++;

  return vtk_polydata_number_of_vertices (pd) - 1;
}

/**
 * @memberof vtkPolyData
 */
int64_t vtk_polydata_add_line (vtkPolyData* pd,
                               int64_t line_point_1,
                               int64_t line_point_2) {

  vtk_polydata_on_add_connectivity (pd);

  if (!vtk_polydata_point_exists (pd, line_point_1) ||
      !vtk_polydata_point_exists (pd, line_point_2)) {
    abort ();
  }

  pd->lines_connectivity[pd->n_lines_connectivity] = line_point_1;
  pd->n_lines_connectivity++;

  pd->lines_connectivity[pd->n_lines_connectivity] = line_point_2;
  pd->n_lines_connectivity++;

  pd->lines_offsets[pd->n_lines_offsets] = pd->n_lines_connectivity;
  pd->n_lines_offsets++;

  return vtk_polydata_number_of_lines (pd) - 1;
}

// int64_t vtk_polydata_add_pointdata (vtkPolyData *pd,
//   const char * pontdata_name,
//   size_t pointdata_ncomp
// ) {

// }

/**
 * @memberof vtkPolyData
 */
int64_t vtk_polydata_add_pointdata_scalar (vtkPolyData* pd, const char* name) {

  vtk_polydata_on_add_data (pd);

  // malloc, get index
  int64_t id = vtk_polydata_malloc_pointdata_scalar (pd);

  // validate name
  vtk_polydata_validate_name (pd, name);

  // set name
  pd->pointdata_names[id] = (char*) malloc (sizeof (char) * strlen (name) + 1);
  strcpy (pd->pointdata_names[id], name);

  return id;
}

/**
 * @memberof vtkPolyData
 */
double* vtk_polydata_get_pointdata_data (vtkPolyData* pd, int64_t field) {
  size_t id = (size_t) field;

  assert (id < pd->n_pointdata);
  assert (pd->pointdata_data[id]);

  return pd->pointdata_data[id];
}

/**
 * @memberof vtkPolyData
 */
int64_t vtk_polydata_add_pointdata_vector (vtkPolyData* pd,
                                           const char* name,
                                           size_t ncomp) {
  vtk_polydata_on_add_data (pd);

  // malloc, get index
  int64_t id = vtk_polydata_malloc_pointdata_vector (pd, ncomp);

  // validate name
  vtk_polydata_validate_name (pd, name);

  // set name
  pd->pointdata_names[id] = (char*) malloc (sizeof (char) * strlen (name) + 1);
  strcpy (pd->pointdata_names[id], name);

  return id;
}

/**
 * @brief
 *
 * @memberof vtkPolyData
 */
void vtk_polydata_malloc_pointdata (vtkPolyData* pd, size_t n) {
  vtk_polydata_free_pointdata (pd);

  pd->m_pointdata = n;
  pd->n_pointdata = 0;

  pd->pointdata_data = (double**) calloc (pd->m_pointdata, sizeof (double*));
  pd->pointdata_names = (char**) calloc (pd->m_pointdata, sizeof (char*));
  pd->pointdata_ncomp = (size_t*) calloc (pd->m_pointdata, sizeof (size_t));
}

/**
 * @brief
 *
 * @memberof vtkPolyData
 */
void vtk_polydata_free_pointdata (vtkPolyData* pd) {
  for (size_t i = 0; i < pd->n_pointdata; i++) {
    vtk_polydata_free_pointdata_field (pd, i);
  }

  free (pd->pointdata_data);
  pd->pointdata_data = NULL;
  free (pd->pointdata_names);
  pd->pointdata_names = NULL;
  free (pd->pointdata_ncomp);
  pd->pointdata_ncomp = NULL;

  pd->m_pointdata = 0;
  pd->n_pointdata = 0;
}

/**
 * @memberof vtkPolyData
 */
int64_t vtk_polydata_malloc_pointdata_scalar (vtkPolyData* pd) {

  assert (pd->n_pointdata < pd->m_pointdata);

  size_t id = pd->n_pointdata;
  pd->n_pointdata++;

  pd->pointdata_ncomp[id] = 1;
  // pd->pointdata_data[id] = (double*) malloc (sizeof (double) * pd->n_points);
  pd->pointdata_data[id] = (double*) calloc (pd->n_points, sizeof (double));

  return id;
}

/**
 * @memberof vtkPolyData
 */
int64_t vtk_polydata_malloc_pointdata_vector (vtkPolyData* pd, size_t ncomp) {
  assert (pd->n_pointdata < pd->m_pointdata);
  assert (ncomp > 0);

  size_t id = pd->n_pointdata;
  pd->n_pointdata++;

  pd->pointdata_ncomp[id] = ncomp;
  pd->pointdata_data[id] =
    (double*) calloc (pd->n_points * pd->pointdata_ncomp[id], sizeof (double));

  return id;
}

/**
 * @memberof vtkPolyData
 */
void vtk_polydata_free_pointdata_field (vtkPolyData* pd, int64_t field) {

  size_t id = (size_t) field;
  assert (id < pd->n_pointdata);

  pd->pointdata_ncomp[id] = 0;

  free (pd->pointdata_names[id]);
  pd->pointdata_names[id] = NULL;

  free (pd->pointdata_data[id]);
  pd->pointdata_data[id] = NULL;
}

/**
 * @memberof vtkPolyData
 */
void vtk_polydata_validate_name (vtkPolyData* pd, const char* name) {

  return;

  // Check if name is already used
  for (size_t i = 0; i < pd->n_pointdata; i++) {
    const char* other_name = pd->pointdata_names[i];
    bool matches = !strcmp (other_name, name);
    assert (other_name && !matches);
  }

  // Check if name follows naming conventions
  for (const char* p = name; *p; ++p) {
    assert (*p != '/');
    assert (*p != '.');
  }

  return;
}