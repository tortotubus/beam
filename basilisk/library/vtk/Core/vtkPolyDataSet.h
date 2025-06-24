#include <stdlib.h>
#include <string.h>

/* ── 1. A list of N points, each with D components (usually D=3) ────────── */
typedef struct {
    int64_t n_points;      /* number of points */
    int64_t n_components;  /* floats per point (e.g. 3 for XYZ) */
    float  *data;         /* length = n_points * n_components */
} vtkPoints_t;

/* ── 2. A generic “cell array”:
      - offsets[i] is the cumulative length up to cell i
      - connectivity is all point‐indices concatenated                      ── */
typedef struct {
    int64_t   n_cells;       /* how many cells of this type */
    int64_t  *offsets;       /* length = n_cells */
    int64_t  *connectivity;  /* length = offsets[n_cells-1] */
} vtkCellArray_t;

/* ── 3. The full PolyData, with one vtkCellArray for each VTK cell type ─── */
typedef struct {
    vtkPoints_t    points;           /* the X,Y,Z (or D) coordinates */
    vtkCellArray_t vertices;         /* VTK_VERTEX cells  (code=1) */
    // vtkCellArray_t poly_vertices;    /* VTK_POLY_VERTEX   (code=2) */
    // vtkCellArray_t lines;            /* VTK_LINE          (code=3) */
    // vtkCellArray_t poly_lines;       /* VTK_POLY_LINE     (code=4) */
    // vtkCellArray_t triangles;        /* VTK_TRIANGLE      (code=5) */
    // vtkCellArray_t triangle_strips;  /* VTK_TRIANGLE_STRIP(code=6) */
    // vtkCellArray_t polygons;         /* VTK_POLYGON       (code=7) */
} vtkPolyData_t;

/* ── 4. Initialization / cleanup API, all as static inline in header ────── */

/** Initialize points container.
 *  @return 0 on success, -1 on malloc failure.
 */
static inline int vtkPoints_init(vtkPoints_t *pts,
                                 int64_t       n_points,
                                 int64_t       n_components)
{
    pts->n_points     = n_points;
    pts->n_components = n_components;
    int64_t total      = n_points * n_components;
    pts->data         = (float*)malloc(total * sizeof(float));
    if (!pts->data) return -1;
    memset(pts->data, 0, total * sizeof(float));
    return 0;
}

/** Free points container. Safe to call on zero‐inited struct. */
static inline void vtkPoints_free(vtkPoints_t *pts)
{
    free(pts->data);
    pts->data = NULL;
    pts->n_points = pts->n_components = 0;
}

/** Initialize a cell‐array given the number of cells and
 *  the number of point‐indices per cell.
 *  cell_sizes[i] = # indices in the i-th cell.
 *  @return 0 on success, -1 on malloc failure.
 */
static inline int vtkCellArray_init(vtkCellArray_t *ca,
                                    int64_t           n_cells,
                                    const int64_t    *cell_sizes)
{
    ca->n_cells = n_cells;
    ca->offsets = (int64_t*)malloc(n_cells * sizeof(int64_t));
    if (!ca->offsets) return -1;
    int64_t total = 0;
    for (int64_t i = 0; i < n_cells; i++) {
        total += cell_sizes[i];
        ca->offsets[i] = total;
    }
    ca->connectivity = (int64_t*)malloc(total * sizeof(int64_t));
    if (!ca->connectivity) {
        free(ca->offsets);
        return -1;
    }
    return 0;
}

/** Free cell‐array. Safe to call on zero‐inited struct. */
static inline void vtkCellArray_free(vtkCellArray_t *ca)
{
    free(ca->offsets);
    free(ca->connectivity);
    ca->offsets      = NULL;
    ca->connectivity = NULL;
    ca->n_cells      = 0;
}

/** Zero‐initialize an empty PolyData. */
static inline void vtkPolyData_init(vtkPolyData_t *pd)
{
    memset(pd, 0, sizeof(*pd));
}

/** Free all internal buffers of a PolyData. */
static inline void vtkPolyData_free(vtkPolyData_t *pd)
{
    vtkPoints_free    (&pd->points);
    vtkCellArray_free (&pd->vertices);
    // vtkCellArray_free (&pd->poly_vertices);
    // vtkCellArray_free (&pd->lines);
    // vtkCellArray_free (&pd->poly_lines);
    // vtkCellArray_free (&pd->triangles);
    // vtkCellArray_free (&pd->triangle_strips);
    // vtkCellArray_free (&pd->polygons);
}

// #endif /* VTKPOLYDATA_H */
