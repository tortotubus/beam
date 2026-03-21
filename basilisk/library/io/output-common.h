#include <time.h>
// #include <unistd.h>
#include <errno.h>
#include <stdio.h>
#include <sys/stat.h> // mkdir
#include <sys/types.h>

#if _MPI
#include <mpi.h>
#endif

bool file_exists (const char* fname) {
  FILE* fp = fopen (fname, "rb");
  if (!fp) {
    return false;
  } else {
    fclose (fp);
    return true;
  }
}

int copy_file (const char* src, const char* dst) {
  FILE* in = fopen (src, "rb");
  if (!in) {
    perror ("fopen src");
    return -1;
  }

  FILE* out = fopen (dst, "wb");
  if (!out) {
    perror ("fopen dst");
    fclose (in);
    return -1;
  }

  char buf[1 << 16]; // 64 KB
  size_t n;

  while ((n = fread (buf, 1, sizeof (buf), in)) > 0) {
    if (fwrite (buf, 1, n, out) != n) {
      perror ("fwrite");
      fclose (in);
      fclose (out);
      return -1;
    }
  }

  if (ferror (in)) {
    perror ("fread");
    fclose (in);
    fclose (out);
    return -1;
  }

  fclose (in);
  fclose (out);
  return 0;
}

int erase_file (const char* fname) {
  if (remove (fname) == 0)
    return 0;
  if (errno == ENOENT)
    return 0;
  return -1;
}

int create_path (const char* path) {
  if (mkdir (path, 0777) == 0)
    return 0;
  if (errno == EEXIST)
    return 0;
  return -1;
}

static char executable_name[4096] = "";

const char* get_executable_name (void) {
  if (executable_name[0] == '\0') {
    ssize_t n = readlink (
      "/proc/self/exe", executable_name, sizeof (executable_name) - 1);
    if (n < 0) {
      // fallback
      return "unknown";
    }
    executable_name[n] = '\0';
  }

  const char* basename = executable_name;
  for (const char* p = executable_name; *p; ++p)
    if (*p == '/')
      basename = p + 1;

  return basename;
}

#if _MPI
const char* get_datetime_string (void) {
  static char datetime[20] = "";
  int rank = 0;

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    time_t now = time (NULL);
    struct tm* tm_now = localtime (&now);

    if (tm_now == NULL)
      datetime[0] = '\0';
    else if (strftime (datetime, sizeof (datetime), "%Y-%m-%d_%H-%M-%S", tm_now) == 0)
      datetime[0] = '\0';
  }

  MPI_Bcast (datetime, sizeof (datetime), MPI_CHAR, 0, MPI_COMM_WORLD);

  return datetime;
}
#else
const char* get_datetime_string (void) {
  static char datetime[20] = "";

  time_t now = time (NULL);
  struct tm* tm_now = localtime (&now);

  if (tm_now == NULL)
    return "";

  if (strftime (datetime, sizeof (datetime), "%Y-%m-%d_%H-%M-%S", tm_now) == 0)
    return "";

  return datetime;
}
#endif

static char unique_basename[4096] = "";

const char* get_unique_basename (void) {
  if (unique_basename[0] == '\0')
    snprintf (unique_basename, sizeof (unique_basename), "%s-%s",
              get_executable_name(), get_datetime_string());

  return unique_basename;
}
