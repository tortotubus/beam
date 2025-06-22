#include "../Utility/foreach_cell_bfs.h"
#include "vtkDataArray.h"
#include "vtkHyperTreeGridData.h"
#include "vtkBase64.h"

typedef struct {
  void *ctx;
  void (*ctx_free)(void *ctx);
  vtkHyperTreeGridData *data;
  vtkBase64Writer *b64_writer;
  vtkType type;
} vtkHyperTreeGridWriterCtx;

void vtk_hypertree_grid_writer_ctx_free(void *ctx) {
  vtkHyperTreeGridWriterCtx *c = (vtkHyperTreeGridWriterCtx *)ctx;
  if (c->ctx != NULL && c->ctx_free != NULL) {
    c->ctx_free(c->ctx);
    free(c->ctx);
  } else if (c->ctx != NULL) {
    free(c->ctx);
  }
}

void vtk_hypertreegrid_x_writer(void *ctx, vtkFormat format,
#if _MPI
                                MPI_Offset base, MPI_File fp
#else
                                FILE *fp
#endif
) {
  vtkHyperTreeGridWriterCtx *c = (vtkHyperTreeGridWriterCtx *)ctx;
  vtkHyperTreeGridData *data = c->data;
#if _MPI
  if (pid() == 0) {
    size_t n_x = data->n_x;
    double *x = data->x;
    MPI_File_write_at(fp, base, x, n_x, MPI_DOUBLE, MPI_STATUS_IGNORE);
  }
#else
  vtkBase64Writer *b64_writer = c->b64_writer;
  NOT_UNUSED(b64_writer);
  size_t n_x = data->n_x;
  double *x = data->x;

  switch (format) {
  case 0: { // ascii
    for (size_t i = 0; i < n_x; i++) {
      fprintf(fp, "%f ", x[i]);
    }
    // fputc('\n', fp);
    break;
  }
  case 1: { // raw
    fwrite(x, sizeof(double), n_x, fp);
    break;
  }
  case 2: { // base64
    vtkBase64WriterWrite(b64_writer, x, n_x * sizeof(double));
    break;
  }
  }
#endif
}

void vtk_hypertreegrid_y_writer(void *ctx, vtkFormat format,
#if _MPI
                                MPI_Offset base, MPI_File fp
#else
                                FILE *fp
#endif
) {
  vtkHyperTreeGridWriterCtx *c = (vtkHyperTreeGridWriterCtx *)ctx;
  vtkHyperTreeGridData *data = c->data;

#if _MPI
  if (pid() == 0) {
    size_t n_y = data->n_y;
    double *y = data->y;
    MPI_File_write_at(fp, base, y, n_y, MPI_DOUBLE, MPI_STATUS_IGNORE);
  }
#else
  vtkBase64Writer *b64_writer = c->b64_writer;
  NOT_UNUSED(b64_writer);
  size_t n_y = data->n_y;
  double *y = data->y;

  switch (format) {
  case 0: { // ascii
    for (size_t i = 0; i < n_y; i++) {
      fprintf(fp, "%f ", y[i]);
    }
    // fputc('\n', fp);
    break;
  }
  case 1: { // raw
    fwrite(y, sizeof(double), n_y, fp);
    break;
  }
  case 2: { // base64
    vtkBase64WriterWrite(b64_writer, y, n_y * sizeof(double));
    break;
  }
  }
#endif
}

void vtk_hypertreegrid_z_writer(void *ctx, vtkFormat format,
#if _MPI
                                MPI_Offset base, MPI_File fp
#else
                                FILE *fp
#endif
) {
  vtkHyperTreeGridWriterCtx *c = (vtkHyperTreeGridWriterCtx *)ctx;
  vtkHyperTreeGridData *data = c->data;

#if _MPI
  if (pid() == 0) {
    size_t n_z = data->n_z;
    double *z = data->z;
    MPI_File_write_at(fp, base, z, n_z, MPI_DOUBLE, MPI_STATUS_IGNORE);
  }
#else
  vtkBase64Writer *b64_writer = c->b64_writer;
  NOT_UNUSED(b64_writer);
  size_t n_z = data->n_z;
  double *z = data->z;

  switch (format) {
  case 0: { // ascii
    for (size_t i = 0; i < n_z; i++) {
      fprintf(fp, "%f ", z[i]);
    }
    break;
  }
  case 1: { // raw
    fwrite(z, sizeof(double), n_z, fp);
    break;
  }
  case 2: { // base64
    vtkBase64WriterWrite(b64_writer, z, n_z * sizeof(double));
    break;
  }
  }
#endif
}

void vtk_hypertreegrid_descriptors_writer(void *ctx, vtkFormat format,
#if _MPI
                                          MPI_Offset base, MPI_File fp
#else
                                          FILE *fp
#endif
) {
  vtkHyperTreeGridWriterCtx *c = (vtkHyperTreeGridWriterCtx *)ctx;
  vtkHyperTreeGridData *data = c->data;

#if _MPI
  if (pid() == 0) {
    Bit_t *descriptors = data->descriptors;
    size_t descriptors_size = data->descriptors_size;
    MPI_File_write_at(fp, base, descriptors, descriptors_size, MPI_BYTE,
                      MPI_STATUS_IGNORE);
  }
#else
  vtkBase64Writer *b64_writer = c->b64_writer;
  NOT_UNUSED(b64_writer);
  switch (format) {
  case 0: { // ascii
    size_t n_descriptors = data->n_descriptors;
    size_t iter = 0;
    foreach_cell_BFS() {
      if (iter < n_descriptors) {
        fprintf(fp, "%i ", is_leaf(cell) ? 0 : 1);
      }
      iter++;
    }
    break;
  }
  case 1: { // raw
    Bit_t *descriptors = data->descriptors;
    size_t descriptors_size = data->descriptors_size;
    fwrite(descriptors, sizeof(Bit_t), descriptors_size, fp);
    break;
  }
  case 2: { // base64
    Bit_t *descriptors = data->descriptors;
    size_t descriptors_size = data->descriptors_size;
    vtkBase64WriterWrite(b64_writer, descriptors,
                         descriptors_size * sizeof(Bit_t));
    break;
  }
  }
#endif
}

void vtk_hypertreegrid_number_of_vertices_per_depth_writer(void *ctx,
                                                           vtkFormat format,
#if _MPI
                                                           MPI_Offset base,
                                                           MPI_File fp
#else
                                                           FILE *fp
#endif
) {
  vtkHyperTreeGridWriterCtx *c = (vtkHyperTreeGridWriterCtx *)ctx;
  vtkHyperTreeGridData *data = c->data;

#if _MPI
  if (pid() == 0) {
    int64_t *number_of_vertices_per_depth = data->number_of_vertices_per_depth;
    size_t max_depth = data->max_depth;
    MPI_File_write_at(fp, base, number_of_vertices_per_depth, max_depth,
                      MPI_INT64_T, MPI_STATUS_IGNORE);
  }
#else
  vtkBase64Writer *b64_writer = c->b64_writer;
  NOT_UNUSED(b64_writer);
  int64_t *number_of_vertices_per_depth = data->number_of_vertices_per_depth;
  size_t max_depth = data->max_depth;
  switch (format) {
  case 0: { // ascii
    for (size_t i = 0; i < max_depth; i++) {
      fprintf(fp, "%li ", number_of_vertices_per_depth[i]);
    }
    break;
  }
  case 1: { // raw
    fwrite(number_of_vertices_per_depth, sizeof(int64_t), max_depth, fp);
    break;
  }
  case 2: { // base64
    vtkBase64WriterWrite(b64_writer, number_of_vertices_per_depth,
                         max_depth * sizeof(int64_t));
    break;
  }
  }
#endif
}

void vtk_hypertreegrid_tree_ids_writer(void *ctx, vtkFormat format,
#if _MPI
                                       MPI_Offset base, MPI_File fp
#else
                                       FILE *fp
#endif
) {
  // vtkHyperTreeGridWriterCtx *c = (vtkHyperTreeGridWriterCtx *)ctx;
  // vtkHyperTreeGridData *data = c->data;

#if _MPI
  switch (format) {
  case 1: { // raw
    if (pid() == 0) {
      int64_t tree_id = 0;
      MPI_File_write_at(fp, base, &tree_id, 1, MPI_INT64_T, MPI_STATUS_IGNORE);
    }
    break;
  }
  default: {
    fprintf(stderr, "Format option \"%s\" is not availible with MPI.\n",
            vtkFormat_to_string(format));
    exit(EXIT_FAILURE);
    break;
  }
  }
#else
  vtkBase64Writer *b64_writer = c->b64_writer;
  NOT_UNUSED(b64_writer);
  // fseek(fp, start, SEEK_SET);
  switch (format) {
  case 0: { // ascii
    fputs("0 ", fp);
    break;
  }
  case 1: { // raw
    int64_t tree_id = 0;
    fwrite(&tree_id, sizeof(int64_t), 1, fp);
    break;
  }
  case 2: { // base64
    int64_t tree_id = 0;
    vtkBase64WriterWrite(b64_writer, &tree_id, sizeof(int64_t));
    break;
  }
  }
#endif
}

void vtk_hypertreegrid_depth_per_tree_writer(void *ctx, vtkFormat format,
#if _MPI
                                             MPI_Offset base, MPI_File fp
#else
                                             FILE *fp
#endif
) {
  vtkHyperTreeGridWriterCtx *c = (vtkHyperTreeGridWriterCtx *)ctx;
  vtkHyperTreeGridData *data = c->data;

#if _MPI
  if (pid() == 0) {
    uint32_t max_depth = (uint32_t)data->max_depth;
    MPI_File_write_at(fp, base, &max_depth, 1, MPI_UINT32_T, MPI_STATUS_IGNORE);
  }
#else
  vtkBase64Writer *b64_writer = c->b64_writer;
  NOT_UNUSED(b64_writer);
  switch (format) {
  case 0: { // ascii
    fprintf(fp, "%lu ", data->max_depth);
    break;
  }
  case 1: { // raw
    uint32_t max_depth = (uint32_t)data->max_depth;
    fwrite(&max_depth, sizeof(uint32_t), 1, fp);
    break;
  }
  case 2: { // base64
    uint32_t max_depth = (uint32_t)data->max_depth;
    vtkBase64WriterWrite(b64_writer, &max_depth, sizeof(uint32_t));
    break;
  }
  }
#endif
}

typedef struct {
  void *ctx;
  scalar s;
} vtkScalarFieldWriterCtx;

void vtk_hypertreegrid_scalar_writer(void *ctx, vtkFormat format,
#if _MPI
                                     MPI_Offset base, MPI_File fp
#else
                                     FILE *fp
#endif
) {
  vtkHyperTreeGridWriterCtx *c = (vtkHyperTreeGridWriterCtx *)ctx;
  vtkType type = c->type;

  vtkScalarFieldWriterCtx *s_ctx = (vtkScalarFieldWriterCtx *)c->ctx;
  scalar s = s_ctx->s;

#if _MPI
  range_t *ranges = c->data->ranges;
  size_t *offsets = c->data->offsets;
  size_t n_ranges = c->data->n_ranges;

  switch (type) {
  case 8: { // Float32_e
    for (size_t r = 0; r < n_ranges; r++) {
      range_t range = ranges[r];
      if (pid() == range.pid) {
        size_t size = (range.end - range.begin + 1);
        Float32_t *vals = malloc(size * sizeof(Float32_t));
        size_t i = 0;
        for (uint64_t mc = range.begin; mc <= range.end; mc++) {
          Point point = morton_to_point(mc);
          {
            vals[i] = (Float32_t)val(s);
            i++;
          }
        }
        MPI_Offset pos = base + ((MPI_Offset)(offsets[r] * sizeof(Float32_t)));
        MPI_File_write_at(fp, pos, vals, size, MPI_FLOAT, MPI_STATUS_IGNORE);
        free(vals);
      }
    }
    break;
  }
  case 9: { // Float64_e
    for (size_t r = 0; r < n_ranges; r++) {
      range_t range = ranges[r];
      if (pid() == range.pid) {
        size_t size = (range.end - range.begin + 1);
        Float64_t *vals = malloc(size * sizeof(Float64_t));
        size_t i = 0;
        for (uint64_t mc = range.begin; mc <= range.end; mc++) {
          Point point = morton_to_point(mc);
          {
            vals[i] = (Float64_t)val(s);
            i++;
          }
        }
        MPI_Offset pos = base + ((MPI_Offset)(offsets[r] * sizeof(Float64_t)));
        MPI_File_write_at(fp, pos, vals, size, MPI_DOUBLE, MPI_STATUS_IGNORE);
        free(vals);
      }
    }
    break;
  }
  default: {
    fprintf(stderr, "Type option \"%s\" is not availible.\n",
            vtkType_to_string(type));
    exit(EXIT_FAILURE);
    break;
  }
  }

#else
  vtkBase64Writer *b64_writer = c->b64_writer;
  NOT_UNUSED(b64_writer);
  switch (format) {
  case 0: { // ascii
    switch (type) {
    case 8: // Float32_e
      foreach_cell_BFS() {
        Float32_t val = (Float32_t)val(s);
        fprintf(fp, "%f ", val);
      }
      break;

    case 9: { // Float64_e
      foreach_cell_BFS() {
        Float64_t val = (Float64_t)val(s);
        fprintf(fp, "%f ", val);
      }
      break;
    }
    default: {
      fprintf(stderr, "Type option \"%s\" is not availible.\n",
              vtkFormat_to_string(type));
      exit(EXIT_FAILURE);
      break;
    }
    }
    break;
  }
  case 1: { // raw
    switch (type) {
    case 8: { // Float32_e
      foreach_cell_BFS() {
        Float32_t val = (Float32_t)val(s);
        fwrite(&val, sizeof(Float32_t), 1, fp);
      }
      break;
    }
    case 9: { // Float64_e
      foreach_cell_BFS() {
        Float64_t val = (Float64_t)val(s);
        fwrite(&val, sizeof(Float64_t), 1, fp);
      }
      break;
    }
    default: {
      fprintf(stderr, "Type option \"%s\" is not availible.\n",
              vtkFormat_to_string(type));
      exit(EXIT_FAILURE);
      break;
    }
    }
    break;
  }
  case 2: { // base64
    switch (type) {
    case 8: { // Float32_e
      foreach_cell_BFS() {
        Float32_t val = (Float32_t)val(s);
        vtkBase64WriterWrite(b64_writer, &val, sizeof(Float32_t));
      }
      break;
    }
    case 9: { // Float64_e
      foreach_cell_BFS() {
        Float64_t val = (Float64_t)val(s);
        vtkBase64WriterWrite(b64_writer, &val, sizeof(Float64_t));
      }
      break;
    }
    default: {
      fprintf(stderr, "Type option \"%s\" is not availible.\n",
              vtkFormat_to_string(type));
      exit(EXIT_FAILURE);
      break;
    }
    }
    break;
  }
  }
#endif
}

typedef struct {
  void *ctx;
  vector v;
} vtkVectorFieldWriterCtx;

void vtk_hypertreegrid_vector_writer(void *ctx, vtkFormat format,
#if _MPI
                                     MPI_Offset base, MPI_File fp
#else
                                     FILE *fp
#endif
) {
  vtkHyperTreeGridWriterCtx *c = (vtkHyperTreeGridWriterCtx *)ctx;
  vtkType type = c->type;

  vtkVectorFieldWriterCtx *v_ctx = (vtkVectorFieldWriterCtx *)c->ctx;
  vector v = v_ctx->v;

#if _MPI
  range_t *ranges = c->data->ranges;
  size_t *offsets = c->data->offsets;
  size_t n_ranges = c->data->n_ranges;

  switch (type) {
  case 8: { // Float32_e
    for (size_t r = 0; r < n_ranges; r++) {
      range_t range = ranges[r];
      if (pid() == range.pid) {
        size_t size = (range.end - range.begin + 1) * dimension;
        Float32_t *vals = malloc(size * sizeof(Float32_t));
        size_t i = 0;
        for (uint64_t mc = range.begin; mc <= range.end; mc++) {
          Point point = morton_to_point(mc);
          {
#if dimension == 1
            vals[i] = (Float32_t)val(v.x);
#elif dimension == 2
            vals[i] = (Float32_t)val(v.x);
            vals[i + 1] = (Float32_t)val(v.y);
#else // dimension == 3
            vals[i] = (Float32_t)val(v.x);
            vals[i + 1] = (Float32_t)val(v.y);
            vals[i + 2] = (Float32_t)val(v.z);
#endif
            i += dimension;
          }
        }

        MPI_Offset pos =
            base + ((MPI_Offset)(offsets[r] * sizeof(Float32_t) * dimension));
        MPI_File_write_at(fp, pos, vals, size, MPI_FLOAT, MPI_STATUS_IGNORE);
        free(vals);
      }
    }
    break;
  }
  case 9: { // Float64_e
    for (size_t r = 0; r < n_ranges; r++) {
      range_t range = ranges[r];
      if (pid() == range.pid) {
        size_t size = (range.end - range.begin + 1) * dimension;
        Float64_t *vals = malloc(size * sizeof(Float64_t));
        size_t i = 0;
        for (uint64_t mc = range.begin; mc <= range.end; mc++) {
          Point point = morton_to_point(mc);
          {
#if dimension == 1
            vals[i] = (Float64_t)val(v.x);
#elif dimension == 2
            vals[i] = (Float64_t)val(v.x);
            vals[i + 1] = (Float64_t)val(v.y);
#else // dimension == 3
            vals[i] = (Float64_t)val(v.x);
            vals[i + 1] = (Float64_t)val(v.y);
            vals[i + 2] = (Float64_t)val(v.z);
#endif
            i += dimension;
          }
        }

        MPI_Offset pos =
            base + ((MPI_Offset)(offsets[r] * sizeof(Float64_t) * dimension));
        MPI_File_write_at(fp, pos, vals, size, MPI_DOUBLE, MPI_STATUS_IGNORE);
        free(vals);
      }
    }
    break;
  }
  default: {
    fprintf(stderr, "Type option \"%s\" is not availible.\n",
            vtkType_to_string(type));
    exit(EXIT_FAILURE);
    break;
  }
  }

#else
  vtkBase64Writer *b64_writer = c->b64_writer;
  NOT_UNUSED(b64_writer);
  switch (format) {
  case 0: { // ascii
    switch (type) {
    case 8: { // Float32_e
      foreach_cell_BFS() {
#if dimension == 1
        Float32_t val_x = (Float32_t)val(v.x);
        fprintf(fp, "%f ", val_x);
#elif dimension == 2
        Float32_t val_x = (Float32_t)val(v.x);
        Float32_t val_y = (Float32_t)val(v.y);
        fprintf(fp, "%f ", val_x);
        fprintf(fp, "%f ", val_y);
#else // dimension == 3
        Float32_t val_x = (Float32_t)val(v.x);
        Float32_t val_y = (Float32_t)val(v.y);
        Float32_t val_z = (Float32_t)val(v.z);
        fprintf(fp, "%f ", val_x);
        fprintf(fp, "%f ", val_y);
        fprintf(fp, "%f ", val_z);
#endif
      }
      break;
    }
    case 9: { // Float64_e
      foreach_cell_BFS() {
#if dimension == 1
        Float64_t val_x = (Float64_t)val(v.x);
        fprintf(fp, "%f ", val_x);
#elif dimension == 2
        Float64_t val_x = (Float64_t)val(v.x);
        Float64_t val_y = (Float64_t)val(v.y);
        fprintf(fp, "%f ", val_x);
        fprintf(fp, "%f ", val_y);
#else // dimension == 3
        Float64_t val_x = (Float64_t)val(v.x);
        Float64_t val_y = (Float64_t)val(v.y);
        Float64_t val_z = (Float64_t)val(v.z);
        fprintf(fp, "%f ", val_x);
        fprintf(fp, "%f ", val_y);
        fprintf(fp, "%f ", val_z);
#endif
      }
      break;
    }
    default: {
      fprintf(stderr, "Type option \"%s\" is not availible.\n",
              vtkType_to_string(type));
      exit(EXIT_FAILURE);
      break;
    }
    }
    break;
  }
  case 1: { // raw
    switch (type) {
    case 8: { // Float32_e
      foreach_cell_BFS() {
#if dimension == 1
        Float32_t val_x = (Float32_t)val(v.x);
        fwrite(&val_x, sizeof(Float32_t), 1, fp);
#elif dimension == 2
        Float32_t val_x = (Float32_t)val(v.x);
        Float32_t val_y = (Float32_t)val(v.y);
        fwrite(&val_x, sizeof(Float32_t), 1, fp);
        fwrite(&val_y, sizeof(Float32_t), 1, fp);
#else // dimension == 3
        Float32_t val_x = (Float32_t)val(v.x);
        Float32_t val_y = (Float32_t)val(v.y);
        Float32_t val_z = (Float32_t)val(v.z);
        fwrite(&val_x, sizeof(Float32_t), 1, fp);
        fwrite(&val_y, sizeof(Float32_t), 1, fp);
        fwrite(&val_z, sizeof(Float32_t), 1, fp);
#endif
      }
      break;
    }
    case 9: { // Float64_e
      foreach_cell_BFS() {
#if dimension == 1
        Float64_t val_x = (Float64_t)val(v.x);
        fwrite(&val_x, sizeof(Float64_t), 1, fp);
#elif dimension == 2
        Float64_t val_x = (Float64_t)val(v.x);
        Float64_t val_y = (Float64_t)val(v.y);
        fwrite(&val_x, sizeof(Float64_t), 1, fp);
        fwrite(&val_y, sizeof(Float64_t), 1, fp);
#else // dimension == 3
        Float64_t val_x = (Float64_t)val(v.x);
        Float64_t val_y = (Float64_t)val(v.y);
        Float64_t val_z = (Float64_t)val(v.z);
        fwrite(&val_x, sizeof(Float64_t), 1, fp);
        fwrite(&val_y, sizeof(Float64_t), 1, fp);
        fwrite(&val_z, sizeof(Float64_t), 1, fp);
#endif
      }
      break;
    }
    default: {
      fprintf(stderr, "Type option \"%s\" is not availible.\n",
              vtkType_to_string(type));
      exit(EXIT_FAILURE);
      break;
    }
    }
    break;
  }
  case 2: { // base64
    switch (type) {
    case 8: { // Float32_e
      foreach_cell_BFS() {
#if dimension == 1
        Float32_t val_x = (Float32_t)val(v.x);
        vtkBase64WriterWrite(b64_writer, &val_x, sizeof(Float32_t));
#elif dimension == 2
        Float32_t val_x = (Float32_t)val(v.x);
        Float32_t val_y = (Float32_t)val(v.y);
        vtkBase64WriterWrite(b64_writer, &val_x, sizeof(Float32_t));
        vtkBase64WriterWrite(b64_writer, &val_y, sizeof(Float32_t));
#else // dimension == 3
        Float32_t val_x = (Float32_t)val(v.x);
        Float32_t val_y = (Float32_t)val(v.y);
        Float32_t val_z = (Float32_t)val(v.z);
        vtkBase64WriterWrite(b64_writer, &val_x, sizeof(Float32_t));
        vtkBase64WriterWrite(b64_writer, &val_y, sizeof(Float32_t));
        vtkBase64WriterWrite(b64_writer, &val_z, sizeof(Float32_t));
#endif
      }
      break;
    }
    case 9: { // Float64_e
      foreach_cell_BFS() {
#if dimension == 1
        Float64_t val_x = (Float64_t)val(v.x);
        vtkBase64WriterWrite(b64_writer, &val_x, sizeof(Float64_t));
#elif dimension == 2
        Float64_t val_x = (Float64_t)val(v.x);
        Float64_t val_y = (Float64_t)val(v.y);
        vtkBase64WriterWrite(b64_writer, &val_x, sizeof(Float64_t));
        vtkBase64WriterWrite(b64_writer, &val_y, sizeof(Float64_t));
#else // dimension == 3
        Float64_t val_x = (Float64_t)val(v.x);
        Float64_t val_y = (Float64_t)val(v.y);
        Float64_t val_z = (Float64_t)val(v.z);
        vtkBase64WriterWrite(b64_writer, &val_x, sizeof(Float64_t));
        vtkBase64WriterWrite(b64_writer, &val_y, sizeof(Float64_t));
        vtkBase64WriterWrite(b64_writer, &val_z, sizeof(Float64_t));
#endif
      }
      break;
    }
    default: {
      fprintf(stderr, "Type option \"%s\" is not availible.\n",
              vtkType_to_string(type));
      exit(EXIT_FAILURE);
      break;
    }
    }
    break;
  }
  }
#endif
}

void vtk_hypertreegrid_time_value_writer(void *ctx, vtkFormat format,
#if _MPI
                                         MPI_Offset base, MPI_File fp
#else
                                         FILE *fp
#endif
) {
  vtkHyperTreeGridWriterCtx *c = (vtkHyperTreeGridWriterCtx *)ctx;
  vtkHyperTreeGridData *data = c->data;

#if _MPI
  if (pid() == 0) {
    MPI_File_write_at(fp, base, &data->time, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
  }
#else
  switch (format) {
  case 0: { // ascii
    fprintf(fp, "%f", data->time);
  } break;
  case 1: { // raw
    fwrite(&data->time, sizeof(double), 1, fp);
    break;
  }
  case 2: { // base64
    vtkBase64WriterWrite(c->b64_writer, &data->time, sizeof(double));
  }
  }
#endif
}

typedef struct {
  vtkDataArray x, y, z;
  vtkDataArray descriptors;
  vtkDataArray number_of_vertices_per_depth;
  vtkDataArray tree_ids;
  vtkDataArray depth_per_tree;
  vtkDataArray time_value;
  vtkDataArrayList scalars, vectors;
  vtkHyperTreeGridData *data;
} vtkHyperTreeGrid;

void vtk_hypertreegrid_free(vtkHyperTreeGrid *vtk_hypertreegrid) {
  if (vtk_hypertreegrid) {
    vtk_data_array_free(&vtk_hypertreegrid->x);
    vtk_data_array_free(&vtk_hypertreegrid->y);
    vtk_data_array_free(&vtk_hypertreegrid->z);
    vtk_data_array_free(&vtk_hypertreegrid->descriptors);
    vtk_data_array_free(&vtk_hypertreegrid->number_of_vertices_per_depth);
    vtk_data_array_free(&vtk_hypertreegrid->tree_ids);
    vtk_data_array_free(&vtk_hypertreegrid->depth_per_tree);
    vtk_data_array_free(&vtk_hypertreegrid->time_value);
    vtk_data_array_list_clear(&vtk_hypertreegrid->scalars);
    vtk_data_array_list_clear(&vtk_hypertreegrid->vectors);
    vtk_hypertreegrid_data_free(vtk_hypertreegrid->data);
    free(vtk_hypertreegrid);
  }
}

vtkHyperTreeGrid *vtk_hypertreegrid_init(void) {
  vtkHyperTreeGrid *vtk_hypertreegrid = calloc(1, sizeof(vtkHyperTreeGrid));
  vtkHyperTreeGridWriterCtx *c = NULL;

  vtk_hypertreegrid->scalars = vtk_data_array_list_init();
  vtk_hypertreegrid->vectors = vtk_data_array_list_init();
  vtk_hypertreegrid->data = vtk_hypertreegrid_data_init();

  vtk_hypertreegrid->x.type = 9; // Float64
  vtk_hypertreegrid->x.name = strdup("XCoordinates");
  vtk_hypertreegrid->x.number_of_tuples = (dimension >= 1) ? 2 : 1;
  vtk_hypertreegrid->x.number_of_components = 0;
  vtk_hypertreegrid->x.writer = vtk_hypertreegrid_x_writer;

  c = malloc(sizeof(vtkHyperTreeGridWriterCtx));
  c->data = vtk_hypertreegrid->data;
  c->type = 9;
  c->ctx_free = NULL;
  c->ctx = NULL;

  vtk_hypertreegrid->x.writer_ctx = c;
  vtk_hypertreegrid->x.writer_ctx_free = vtk_hypertree_grid_writer_ctx_free;
  c = NULL;

  vtk_hypertreegrid->y.type = 9; // Float64
  vtk_hypertreegrid->y.name = strdup("YCoordinates");
  vtk_hypertreegrid->y.number_of_tuples = (dimension >= 2) ? 2 : 1;
  vtk_hypertreegrid->y.number_of_components = 0;
  vtk_hypertreegrid->y.writer = vtk_hypertreegrid_y_writer;

  c = malloc(sizeof(vtkHyperTreeGridWriterCtx));
  c->data = vtk_hypertreegrid->data;
  c->type = 9;
  c->ctx_free = NULL;
  c->ctx = NULL;

  vtk_hypertreegrid->y.writer_ctx = c;
  vtk_hypertreegrid->y.writer_ctx_free = vtk_hypertree_grid_writer_ctx_free;
  c = NULL;

  vtk_hypertreegrid->z.type = 9; // Float64
  vtk_hypertreegrid->z.name =
      strdup("ZCoordinates"); // replaced direct assignment with strdup
  vtk_hypertreegrid->z.number_of_tuples = (dimension >= 3) ? 2 : 1;
  vtk_hypertreegrid->z.number_of_components = 0;
  vtk_hypertreegrid->z.writer = vtk_hypertreegrid_z_writer;

  c = malloc(sizeof(vtkHyperTreeGridWriterCtx));
  c->data = vtk_hypertreegrid->data;
  c->type = 9;
  c->ctx_free = NULL;
  c->ctx = NULL;

  vtk_hypertreegrid->z.writer_ctx = c;
  vtk_hypertreegrid->z.writer_ctx_free = vtk_hypertree_grid_writer_ctx_free;
  c = NULL;

  vtk_hypertreegrid->descriptors.type = 10;                    // Bit
  vtk_hypertreegrid->descriptors.name = strdup("Descriptors"); // replaced
  vtk_hypertreegrid->descriptors.number_of_tuples =
      vtk_hypertreegrid->data->n_descriptors;
  vtk_hypertreegrid->descriptors.number_of_components = 0;
  vtk_hypertreegrid->descriptors.writer = vtk_hypertreegrid_descriptors_writer;

  c = malloc(sizeof(vtkHyperTreeGridWriterCtx));
  c->data = vtk_hypertreegrid->data;
  c->type = 10;
  c->ctx_free = NULL;
  c->ctx = NULL;

  vtk_hypertreegrid->descriptors.writer_ctx = c;
  vtk_hypertreegrid->descriptors.writer_ctx_free =
      vtk_hypertree_grid_writer_ctx_free;

  c = NULL;

  vtk_hypertreegrid->number_of_vertices_per_depth.type = 7; // Int64
  vtk_hypertreegrid->number_of_vertices_per_depth.name =
      strdup("NumberOfVerticesPerDepth"); // replaced
  vtk_hypertreegrid->number_of_vertices_per_depth.number_of_tuples =
      vtk_hypertreegrid->data->max_depth;
  vtk_hypertreegrid->number_of_vertices_per_depth.number_of_components = 0;
  vtk_hypertreegrid->number_of_vertices_per_depth.writer =
      vtk_hypertreegrid_number_of_vertices_per_depth_writer;

  c = malloc(sizeof(vtkHyperTreeGridWriterCtx));
  c->data = vtk_hypertreegrid->data;
  c->type = 7;
  c->ctx_free = NULL;
  c->ctx = NULL;

  vtk_hypertreegrid->number_of_vertices_per_depth.writer_ctx = c;
  vtk_hypertreegrid->number_of_vertices_per_depth.writer_ctx_free =
      vtk_hypertree_grid_writer_ctx_free;

  c = NULL;

  vtk_hypertreegrid->tree_ids.type = 7;                 // Int64
  vtk_hypertreegrid->tree_ids.name = strdup("TreeIds"); // replaced
  vtk_hypertreegrid->tree_ids.number_of_tuples = 1;
  vtk_hypertreegrid->tree_ids.number_of_components = 0;
  vtk_hypertreegrid->tree_ids.writer = vtk_hypertreegrid_tree_ids_writer;

  c = malloc(sizeof(vtkHyperTreeGridWriterCtx));
  c->data = vtk_hypertreegrid->data;
  c->type = 7;
  c->ctx_free = NULL;
  c->ctx = NULL;

  vtk_hypertreegrid->tree_ids.writer_ctx = c;
  vtk_hypertreegrid->tree_ids.writer_ctx_free =
      vtk_hypertree_grid_writer_ctx_free;
  c = NULL;

  vtk_hypertreegrid->depth_per_tree.type = 2;                      // UInt32
  vtk_hypertreegrid->depth_per_tree.name = strdup("DepthPerTree"); // replaced
  vtk_hypertreegrid->depth_per_tree.number_of_tuples = 1;
  vtk_hypertreegrid->depth_per_tree.number_of_components = 0;
  vtk_hypertreegrid->depth_per_tree.writer =
      vtk_hypertreegrid_depth_per_tree_writer;

  c = malloc(sizeof(vtkHyperTreeGridWriterCtx));
  c->data = vtk_hypertreegrid->data;
  c->type = 2;
  c->ctx_free = NULL;
  c->ctx = NULL;

  vtk_hypertreegrid->depth_per_tree.writer_ctx = c;
  vtk_hypertreegrid->depth_per_tree.writer_ctx_free =
      vtk_hypertree_grid_writer_ctx_free;
  c = NULL;

  vtk_hypertreegrid->time_value.type = 9;                   // Float64
  vtk_hypertreegrid->time_value.name = strdup("TimeValue"); // replaced
  vtk_hypertreegrid->time_value.number_of_tuples = 1;
  vtk_hypertreegrid->time_value.number_of_components = 0;
  vtk_hypertreegrid->time_value.writer = vtk_hypertreegrid_time_value_writer;

  c = malloc(sizeof(vtkHyperTreeGridWriterCtx));
  c->data = vtk_hypertreegrid->data;
  c->type = 9;
  c->ctx_free = NULL;
  c->ctx = NULL;

  vtk_hypertreegrid->time_value.writer_ctx = c;
  vtk_hypertreegrid->time_value.writer_ctx_free =
      vtk_hypertree_grid_writer_ctx_free;
  c = NULL;

  return vtk_hypertreegrid;
}

vtkDataArray *vtk_hypertreegrid_add_scalar(vtkHyperTreeGrid *vtk_hypertreegrid,
                                           scalar s, vtkType type) {

  vtkHyperTreeGridWriterCtx *c = malloc(sizeof(vtkHyperTreeGridWriterCtx));

  c->data = vtk_hypertreegrid->data;
  c->type = type;
  c->b64_writer = NULL;

  vtkScalarFieldWriterCtx *s_ctx = malloc(sizeof(vtkScalarFieldWriterCtx));
  s_ctx->s = s;
  c->ctx = s_ctx;
  c->ctx_free = NULL;

  vtkDataArray *vtk_data_array_scalar = malloc(sizeof(vtkDataArray));
  vtk_data_array_scalar->type = type;
  vtk_data_array_scalar->name = strdup(s.name);
  vtk_data_array_scalar->writer_ctx = c;
  vtk_data_array_scalar->writer_ctx_free = vtk_hypertree_grid_writer_ctx_free;
  vtk_data_array_scalar->writer = vtk_hypertreegrid_scalar_writer;

  vtk_data_array_list_add(&vtk_hypertreegrid->scalars, vtk_data_array_scalar,
                          true);

  return vtk_data_array_scalar;
}

vtkDataArray *vtk_hypertreegrid_add_vector(vtkHyperTreeGrid *vtk_hypertreegrid,
                                           vector v, vtkType type) {

  vtkHyperTreeGridWriterCtx *c = malloc(sizeof(vtkHyperTreeGridWriterCtx));

  c->data = vtk_hypertreegrid->data;
  c->type = type;
  c->b64_writer = NULL;

  vtkVectorFieldWriterCtx *v_ctx = malloc(sizeof(vtkVectorFieldWriterCtx));
  v_ctx->v = v;
  c->ctx = v_ctx;
  c->ctx_free = NULL;

  vtkDataArray *vtk_data_array_vector = malloc(sizeof(vtkDataArray));
  vtk_data_array_vector->type = type;
  vtk_data_array_vector->writer_ctx = c;
  vtk_data_array_vector->writer_ctx_free = vtk_hypertree_grid_writer_ctx_free;
  vtk_data_array_vector->writer = vtk_hypertreegrid_vector_writer;

  size_t trunc_len = (size_t)(strlen(v.x.name) - 2);
  vtk_data_array_vector->name = malloc((trunc_len + 1) * sizeof(char));
  strncpy(vtk_data_array_vector->name, v.x.name, trunc_len);
  vtk_data_array_vector->name[trunc_len] = '\0';

  vtk_data_array_list_add(&vtk_hypertreegrid->scalars, vtk_data_array_vector,
                          true);

  return vtk_data_array_vector;
}