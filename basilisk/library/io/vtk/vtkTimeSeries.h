#include "library/io/output-common.h"
#include <string.h>

trace void output_vtk_series (const char* entry,
                              const char* series_filename,
                              int iter = i,
                              double time = t) {
  if (pid () == 0) {
    static const char series_tail[] = "\n\t]\n}";
    const long series_tail_len = (long) strlen (series_tail);

    char* series_filename_incompl =
      (char*) malloc ((strlen (series_filename) + 2) * sizeof (char));
    strcpy (series_filename_incompl, series_filename);
    strcat (series_filename_incompl, "~");

    bool series_exists = file_exists (series_filename);
    bool series_incompl_exists = file_exists (series_filename_incompl);

    if (!series_incompl_exists && !series_exists) {
      FILE* fp_incompl = fopen (series_filename_incompl, "w");
      if (!fp_incompl) {
        perror ("output_vtk_series(): fopen incomplete");
        free (series_filename_incompl);
        return;
      }
      fprintf (fp_incompl,
               "{\n"
               "\t\"file-series-version\" : \"1.0\",\n"
               "\t\"files\" : [\n"
               "\t\t{ \"iter\" : %d, \"time\" : %f, \"name\" : \"%s\" }",
               iter,
               time,
               entry);
      fclose (fp_incompl);
    } else {
      if (!series_incompl_exists && series_exists) {
        FILE* fp = fopen (series_filename, "rb");
        FILE* fp_incompl = fopen (series_filename_incompl, "wb");
        if (!fp || !fp_incompl) {
          if (fp)
            fclose (fp);
          if (fp_incompl)
            fclose (fp_incompl);
          perror ("output_vtk_series(): reopen series");
          free (series_filename_incompl);
          return;
        }

        if (fseek (fp, 0, SEEK_END) != 0) {
          perror ("output_vtk_series(): fseek");
          fclose (fp);
          fclose (fp_incompl);
          free (series_filename_incompl);
          return;
        }

        long file_size = ftell (fp);
        if (file_size < 0) {
          perror ("output_vtk_series(): ftell");
          fclose (fp);
          fclose (fp_incompl);
          free (series_filename_incompl);
          return;
        }

        if (fseek (fp, 0, SEEK_SET) != 0) {
          perror ("output_vtk_series(): fseek");
          fclose (fp);
          fclose (fp_incompl);
          free (series_filename_incompl);
          return;
        }

        long keep_size = file_size - series_tail_len;
        if (keep_size < 0)
          keep_size = file_size;

        char buf[1 << 12];
        while (keep_size > 0) {
          size_t chunk =
            (size_t) (keep_size < (long) sizeof (buf) ? keep_size : (long) sizeof (buf));
          if (fread (buf, 1, chunk, fp) != chunk) {
            perror ("output_vtk_series(): fread");
            fclose (fp);
            fclose (fp_incompl);
            free (series_filename_incompl);
            return;
          }
          if (fwrite (buf, 1, chunk, fp_incompl) != chunk) {
            perror ("output_vtk_series(): fwrite");
            fclose (fp);
            fclose (fp_incompl);
            free (series_filename_incompl);
            return;
          }
          keep_size -= (long) chunk;
        }

        fclose (fp);
        fclose (fp_incompl);
      }

      FILE* fp_incompl = fopen (series_filename_incompl, "a");
      if (!fp_incompl) {
        perror ("output_vtk_series(): fopen incomplete");
        free (series_filename_incompl);
        return;
      }
      fprintf (fp_incompl,
               ",\n\t\t{ \"iter\" : %d, \"time\" : %f, \"name\" : \"%s\" }",
               iter,
               time,
               entry);
      fclose (fp_incompl);
    }

    if (copy_file (series_filename_incompl, series_filename) != 0) {
      perror ("output_vtk_series(): copy_file");
      free (series_filename_incompl);
      return;
    }

    FILE* fp = fopen (series_filename, "a");
    if (!fp) {
      perror ("output_vtk_series(): fopen series");
      free (series_filename_incompl);
      return;
    }
    fprintf (fp, "%s", series_tail);
    fclose (fp);

    free (series_filename_incompl);
  }
}
