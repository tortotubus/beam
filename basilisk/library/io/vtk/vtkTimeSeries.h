#include "library/io/output-common.h"

// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>
// #include <unistd.h>

static inline void
vtk_series_index_filename (char* dst,
                           size_t n,
                           const char* series_filename)
{
  snprintf (dst, n, "%s.idx", series_filename);
}

trace void output_vtk_series (const char* entry,
                               const char* series_filename,
                               int iter = i,
                               double time = t)
{
  if (pid () != 0)
    return;

  char idx_filename[1024];
  vtk_series_index_filename (
    idx_filename, sizeof (idx_filename), series_filename);

  FILE* fp = fopen (idx_filename, "ab+");
  if (!fp) {
    perror ("output_vtk_series(): fopen idx");
    return;
  }

  long keep_size = 0;

  if (fseek (fp, 0, SEEK_END) != 0) {
    perror ("output_vtk_series(): fseek idx end");
    fclose (fp);
    return;
  }

  long pos = ftell (fp);
  if (pos < 0) {
    perror ("output_vtk_series(): ftell idx");
    fclose (fp);
    return;
  }

  pos--;

  // Walk backward through the ledger until we find the last entry with
  // idx_iter < iter. Any later entries are discarded before we append.
  while (pos >= 0) {
    int ch = 0;

    if (fseek (fp, pos, SEEK_SET) != 0) {
      perror ("output_vtk_series(): fseek idx");
      fclose (fp);
      return;
    }

    ch = fgetc (fp);
    if (ch == '\n') {
      pos--;
      continue;
    }

    long line_end = pos + 1;

    while (pos >= 0) {
      if (fseek (fp, pos, SEEK_SET) != 0) {
        perror ("output_vtk_series(): fseek idx");
        fclose (fp);
        return;
      }

      ch = fgetc (fp);
      if (ch == '\n')
        break;
      pos--;
    }

    long line_start = pos + 1;
    long line_len = line_end - line_start;

    if (line_len <= 0) {
      pos--;
      continue;
    }

    char* line = (char*) malloc ((size_t) line_len + 1);
    if (!line) {
      perror ("output_vtk_series(): malloc line");
      fclose (fp);
      return;
    }

    if (fseek (fp, line_start, SEEK_SET) != 0) {
      perror ("output_vtk_series(): fseek line");
      free (line);
      fclose (fp);
      return;
    }

    if (fread (line, 1, (size_t) line_len, fp) != (size_t) line_len) {
      perror ("output_vtk_series(): fread line");
      free (line);
      fclose (fp);
      return;
    }

    line[line_len] = '\0';

    int idx_iter = -1;
    if (sscanf (line, "%d", &idx_iter) == 1 && idx_iter < iter) {
      keep_size = line_end;
      free (line);
      break;
    }

    free (line);
    pos = line_start - 2;
  }

  if (fflush (fp) != 0) {
    perror ("output_vtk_series(): fflush idx");
    fclose (fp);
    return;
  }

  if (ftruncate (fileno (fp), keep_size) != 0) {
    perror ("output_vtk_series(): ftruncate idx");
    fclose (fp);
    return;
  }

  if (fseek (fp, 0, SEEK_END) != 0) {
    perror ("output_vtk_series(): fseek append");
    fclose (fp);
    return;
  }

  if (keep_size > 0)
    fputc ('\n', fp);

  fprintf (fp, "%d\t%.17g\t%s\n", iter, time, entry);
  fclose (fp);

  char series_filename_incompl[1024];
  snprintf (
    series_filename_incompl, sizeof (series_filename_incompl), "%s~", series_filename);

  FILE* fp_idx = fopen (idx_filename, "rb");
  if (!fp_idx) {
    perror ("output_vtk_series(): reopen idx");
    return;
  }

  FILE* fp_series = fopen (series_filename_incompl, "wb");
  if (!fp_series) {
    perror ("output_vtk_series(): fopen series");
    fclose (fp_idx);
    return;
  }

  fprintf (fp_series,
           "{\n"
           "\t\"file-series-version\" : \"1.0\",\n"
           "\t\"files\" : [");

  bool first = true;
  char line[4096];
  while (fgets (line, sizeof (line), fp_idx)) {
    int idx_iter = -1;
    double idx_time = 0.;
    char idx_entry[3072];

    if (sscanf (line, "%d\t%lf\t%3071s", &idx_iter, &idx_time, idx_entry) != 3)
      continue;

    fprintf (fp_series,
             "%s\n\t\t{ \"iter\" : %d, \"time\" : %.17g, \"name\" : \"%s\" }",
             first ? "" : ",",
             idx_iter,
             idx_time,
             idx_entry);
    first = false;
  }

  fprintf (fp_series, "\n\t]\n}\n");

  fclose (fp_idx);
  fclose (fp_series);

  if (rename (series_filename_incompl, series_filename) != 0) {
    perror ("output_vtk_series(): rename series");
    return;
  }
}
